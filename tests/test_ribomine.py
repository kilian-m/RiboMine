"""Tests for the machinery that could corrupt a run quietly: the config guard, the
stage graph, the query's text filter, the hand-off to fqdissect, and the filters.
(The read-structure caller itself is tested in fqdissect.)
"""
from __future__ import annotations

import json
import os
import shutil

import pytest

from ribomine import config as cfgmod
from ribomine.config import ConfigError, stages_to_run
from ribomine.sra import query


# --- config ---------------------------------------------------------------
def test_unknown_key_is_an_error(tmp_path):
    """A misspelt config key raises instead of being silently ignored."""
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"qc": {"periodic_minn": 0.9}}))
    with pytest.raises(ConfigError, match="qc.periodic_minn"):
        cfgmod.load(str(p))


def test_defaults_are_complete_and_merge():
    cfg = cfgmod.load(None, {"qc": {"periodic_min": 0.9}})
    assert cfg["qc.periodic_min"] == 0.9          # override applied
    assert cfg["qc.tvd_min"] == 0.10              # sibling default survives
    assert cfg["process.umi_dedup"] is False      # dedup is off by default


def test_the_default_dedup_tool_can_run_the_default_dedup_method():
    """The default dedup tool implements the default method (not every backend
    has every method)."""
    from ribomine.process.dedup import TOOLS, UMICOLLAPSE_ALGO

    cfg = cfgmod.load(None)
    tool, method = cfg["process.umi_dedup_tool"], cfg["process.umi_dedup_method"]
    assert tool in TOOLS
    if tool == "umicollapse":
        assert method in UMICOLLAPSE_ALGO


def test_architecture_thresholds_reach_fqdissect():
    """The `architecture` config block mirrors fqdissect's thresholds, and an
    override reaches the call."""
    from dataclasses import fields

    from fqdissect.infer import Thresholds

    from ribomine import architecture
    names = {f.name for f in fields(Thresholds)}
    assert names | {"process_undetermined"} == set(cfgmod.DEFAULTS["architecture"])

    cfg = cfgmod.load(None, {"architecture": {"min_reads": 7}})
    call = architecture.call({"label": "x", "n_used": 5}, cfg)
    assert call["thresholds"]["min_reads"] == 7
    # too few alignments: undetermined, with a reason
    assert call["status"] == "undetermined" and "5" in call["reason"]
    assert call["structure"].startswith("UNDETERMINED")


@pytest.mark.parametrize("start,end,expect", [
    ("query", "bam", ["query", "qc", "architecture", "bam"]),
    ("query", "qc", ["query", "qc"]),
    ("accessions", "architecture", ["qc", "architecture"]),
    ("fastq", "bam", ["qc", "architecture", "bam"]),
    ("fastq", "qc", ["qc"]),
])
def test_stage_graph(start, end, expect):
    assert stages_to_run(start, end) == expect


# --- the query's text filter ----------------------------------------------
def _row(**kw):
    base = {"study_title": "", "experiment_title": "", "sample_title": "",
            "library_construction_protocol": "", "description": ""}
    return {**base, **kw}


@pytest.mark.parametrize("title", [
    "Ribosome profiling of HeLa cells",
    "Ribo-seq in primary neurons",
    "RiboSeq_HEK293_rep1",
    "ribosome footprinting after CNOT1 depletion",
    "Ribosome-protected fragments (RPFs) from liver",
    "QTI-seq maps start codons",
])
def test_strong_hits(title):
    assert query.classify(_row(study_title=title)) == "strong"


def test_rrna_depletion_never_vetoes_a_real_riboseq_study():
    """Ribo-seq protocols also deplete rRNA, so the 'rRNA depletion' exclude term
    must not veto a strong hit."""
    r = _row(study_title="Ribosome profiling with CNOT1 depletion",
             library_construction_protocol="RNA was rRNA-depleted using Ribo-Zero, then ...")
    assert query.classify(r) == "strong"


def test_exclude_still_vetoes_a_weak_hit():
    """A plain rRNA-depleted RNA-seq run mentioning 'translatome' is not Ribo-seq."""
    r = _row(study_title="Translatome-adjacent expression atlas",
             library_construction_protocol="Ribo-Zero rRNA depletion, TruSeq stranded")
    assert query.classify(r) == "excluded"


def test_plain_rnaseq_is_not_a_candidate():
    assert query.classify(_row(study_title="RNA-seq of airway smooth muscle")) == "none"


def test_ena_rejects_multiword_terms():
    """ENA's text index is tokenised, so a multi-word term matches nothing; it
    must raise."""
    cfg = cfgmod.load(None, {"query": {"terms": ["ribosome profiling"]}})
    with pytest.raises(ValueError, match="single words"):
        query._ena_query(cfg)


# --- trimming (fqdissect + cutadapt) ----------------------------------------
needs_cutadapt = pytest.mark.skipif(not shutil.which("cutadapt"), reason="cutadapt not on PATH")

FOOT = "ACGTACGTACGTACGTACGTACGTACGTAC"          # 30 nt "footprint"
ADAP = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
CALL = {
    "status": "ok", "footprint_len_mode": 30, "polyA_tail": "none",
    "p5_layout": [{"role": "umi5", "offset": 0, "len": 2}],
    "p3_layout": [{"role": "umi3", "len": 5},
                  {"role": "barcode3", "len": 5, "seq": "AGCTA"}],
    "adapter3_name": "illumina_truseq", "adapter3_seq": ADAP,
    "functional": {"trim_5p": 2, "dedup_umi_len": 7},
}


def _fastq(tmp_path, reads):
    p = tmp_path / "in.fastq"
    p.write_text("".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads)))
    return str(p)


def _read(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


@needs_cutadapt
def test_trimming_leaves_the_footprint_and_moves_the_umi_to_the_read_name(tmp_path):
    from ribomine import architecture

    cfg = cfgmod.load(None, {"project": {"threads": 1}})
    with_adapter = "GG" + FOOT + "TTTTT" + "AGCTA" + ADAP    # umi5=GG, umi3=TTTTT
    no_adapter = "GG" + FOOT + "CCCCCCCCCCCCCCCC"            # insert ran off the read
    out = str(tmp_path / "o.trimmed.fastq")
    st = architecture.trim_fastq(CALL, _fastq(tmp_path, [with_adapter, no_adapter]), out, cfg)

    (name, seq), = _read(out)             # the read without the adapter is dropped
    assert seq == FOOT
    assert name.endswith("_GGTTTTT")      # 5' UMI + 3' UMI: the form the deduplicators read
    # the columns mapping_summary.tsv reads from the trim record
    assert (st["n_reads_in"], st["n_reads_out"]) == (2, 1)
    assert st["frac_no_adapter"] == 0.5
    assert st["mean_len_in"] > st["mean_len_out"] == 30.0


@needs_cutadapt
def test_a_barcode_between_two_umi_blocks_stays_out_of_the_umi(tmp_path):
    from ribomine import architecture

    cfg = cfgmod.load(None, {"project": {"threads": 1}})
    call = dict(CALL, p5_layout=[
        {"role": "umi5", "offset": 0, "len": 2},
        {"role": "barcode5", "offset": 2, "len": 3, "seq": "GAT"},
        {"role": "umi5", "offset": 5, "len": 2},
    ], functional={"trim_5p": 7, "dedup_umi_len": 9})
    read = "AC" + "GAT" + "TG" + FOOT + "TTTTT" + "AGCTA" + ADAP
    out = str(tmp_path / "o.trimmed.fastq")
    architecture.trim_fastq(call, _fastq(tmp_path, [read]), out, cfg)
    (name, seq), = _read(out)
    assert seq == FOOT
    umi = name.split("_")[-1]             # all three UMI blocks, and not the barcode GAT
    assert sorted(umi) == sorted("AC" + "TG" + "TTTTT")


# --- shipped reference data ------------------------------------------------
def test_contaminant_fasta_is_bundled():
    """The contaminant reference ships with the package and is the default, so
    the filter is never skipped for lack of one."""
    from ribomine import data
    from ribomine.qc import contaminants

    p = data.human_contaminants()
    with open(p) as fh:
        heads = [ln for ln in fh if ln.startswith(">")]
    assert len(heads) > 3000
    blob = "".join(heads)
    for kind in ("tRNA", "snRNA", "snoRNA", "rRNA"):
        assert kind in blob, f"{kind} missing from the bundled contaminant reference"
    # The pre-rRNA and the rDNA repeat carry the transcribed spacers (ITS1/2, 5'/3'ETS),
    # which no mature rRNA sequence contains; without them those fragments reach the
    # aligner and come back as multimappers.
    for acc in ("NR_046235.3", "U13369.1"):
        assert acc in blob, f"{acc} (pre-rRNA/rDNA) missing from the contaminant reference"

    # and an unconfigured run resolves to it rather than skipping the filter
    cfg = cfgmod.load(None)
    assert contaminants.resolve_fasta(cfg) == p


def test_explicit_contaminant_fasta_wins(tmp_path):
    fa = tmp_path / "mine.fa"
    fa.write_text(">x\nACGT\n")
    from ribomine.qc import contaminants
    cfg = cfgmod.load(None, {"reference": {"contaminant_fasta": str(fa)}})
    assert contaminants.resolve_fasta(cfg) == str(fa)


def test_missing_contaminant_fasta_is_an_error():
    from ribomine.qc import contaminants
    cfg = cfgmod.load(None, {"reference": {"contaminant_fasta": "/nope/absent.fa"}})
    with pytest.raises(ValueError, match="not found"):
        contaminants.resolve_fasta(cfg)


def test_no_option_can_leave_a_bam_unindexed():
    """Every BAM left on disk is sorted and indexed; no config key turns that off."""
    assert "sort_index_bam" not in cfgmod.DEFAULTS["process"]


def test_null_means_is_documented_for_every_nullable_key():
    """For these keys `null` means 'determine it automatically', so every nullable
    default needs an explanation in the generated config."""
    def nullable(d, prefix=""):
        out = []
        for k, v in d.items():
            if isinstance(v, dict):
                out += nullable(v, f"{prefix}{k}.")
            elif v is None:
                out.append(f"{prefix}{k}")
        return out

    undocumented = [k for k in nullable(cfgmod.DEFAULTS) if k not in cfgmod.NULL_MEANS]
    assert not undocumented, f"nullable keys with no explanation of what null means: {undocumented}"


def test_generated_config_explains_itself_and_still_loads(tmp_path):
    p = tmp_path / "c.json"
    cfgmod.write_example(str(p))
    raw = json.loads(p.read_text())
    assert "bundled" in raw["_null_means"]["reference.contaminant_fasta"]
    # the comment keys must not trip the unknown-key guard
    cfg = cfgmod.load(str(p))
    assert cfg["qc.periodic_min"] == cfgmod.DEFAULTS["qc"]["periodic_min"]


def test_query_does_not_demand_the_accession_list_it_may_produce():
    """`ribomine query` does not validate pipeline.start's inputs: it may be what
    produces the accession list. `run` still does."""
    cfg = cfgmod.load(None, {"pipeline": {"start": "accessions",
                                          "accession_list": "does/not/exist.txt"}})
    cfg.validate(need_reference=False, need_inputs=False)      # must not raise
    with pytest.raises(ConfigError, match="accession list not found"):
        cfg.validate(need_reference=False)                     # but `run` still checks


# --- gene counts / the read-count matrix ------------------------------------
READS_PER_GENE = "\n".join([
    # STAR's four summary rows, then the genes.
    # columns: gene_id, unstranded, sense, antisense
    "N_unmapped\t100\t100\t100",
    "N_multimapping\t50\t50\t50",
    "N_noFeature\t30\t35\t900",
    "N_ambiguous\t10\t8\t4",
    "ENSG01\t210\t200\t10",
    "ENSG02\t0\t0\t0",
    "ENSG03\t95\t90\t5",
]) + "\n"

GENE_INFO = "\n".join([        # STAR's own gene table: a count, then id/name/biotype
    "3",
    "ENSG01\tAAA\tprotein_coding",
    "ENSG02\tBBB\tlincRNA",
    "ENSG03\tCCC\tprotein_coding",
]) + "\n"


def _star_index(tmp_path):
    d = tmp_path / "index"
    d.mkdir()
    (d / "geneInfo.tab").write_text(GENE_INFO)
    return str(d)


def _counts_tab(tmp_path, name: str, text: str = READS_PER_GENE) -> str:
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_gene_counts_are_read_off_the_sense_strand(tmp_path):
    """Counts come from the sense column (a footprint maps to the transcript's
    own strand), with the summary statistics alongside."""
    from ribomine.process import counts

    c, stats = counts.read_counts(_counts_tab(tmp_path, "r.tab"))
    assert c == {"ENSG01": 200, "ENSG02": 0, "ENSG03": 90}   # sense column, not 210/95
    assert stats["n_in_genes"] == 290
    assert stats["n_genes_detected"] == 2                    # the zero gene is not "detected"
    # reads outside the counts, which complete the denominator
    assert stats["n_no_feature"] == 35 and stats["n_ambiguous"] == 8
    assert stats["frac_in_genes"] == round(290 / (290 + 35 + 8), 4)
    # STAR's N_multimapping row is 0 when multimap_nmax=1 (multimappers are dropped
    # before counting), so it is not reported; Log.final.out has the rate
    assert "n_multimapping" not in stats
    assert stats["sense_over_antisense"] == round(290 / 15, 1)


def test_a_library_that_is_not_sense_stranded_is_called_out(tmp_path, caplog):
    """A library with similar sense and antisense counts triggers a warning: its
    sense column undercounts it."""
    from ribomine.process import counts

    flipped = "\n".join([
        "N_noFeature\t0\t0\t0",
        "N_ambiguous\t0\t0\t0",
        "ENSG01\t2000\t1000\t1000",     # 1:1 -- unstranded or reversed
    ]) + "\n"
    with caplog.at_level("WARNING"):
        _, stats = counts.read_counts(_counts_tab(tmp_path, "f.tab", flipped), label="SRRX")
    assert stats["sense_over_antisense"] == 1.0
    assert "does not look sense-stranded" in caplog.text


def test_the_matrix_has_a_row_for_every_gene_including_the_zero_ones(tmp_path):
    """Rows come from the annotation, not from the data, so every matrix has the
    same row set."""
    from ribomine.process import counts
    from ribomine.utils import read_tsv

    out = str(tmp_path / "m.tsv")
    counts.matrix(out, _star_index(tmp_path),
                  [("SRR1", _counts_tab(tmp_path, "a.tab")),
                   ("SRR2", _counts_tab(tmp_path, "b.tab"))])
    rows = read_tsv(out)
    assert [r["gene_id"] for r in rows] == ["ENSG01", "ENSG02", "ENSG03"]
    assert rows[0]["gene_name"] == "AAA"          # names come from the STAR index
    assert rows[0]["SRR1"] == "200" and rows[0]["SRR2"] == "200"
    assert rows[1]["SRR1"] == "0"                 # a gene with no reads is a 0, not a gap


def test_a_run_with_no_counts_is_left_out_rather_than_left_blank(tmp_path):
    """A run without counts is omitted as a column, not written as empty cells."""
    from ribomine.process import counts
    from ribomine.utils import read_tsv

    out = str(tmp_path / "m.tsv")
    counts.matrix(out, _star_index(tmp_path),
                  [("SRR1", _counts_tab(tmp_path, "a.tab")),
                   ("SRR_MISSING", str(tmp_path / "nope.tab"))])
    rows = read_tsv(out)
    assert "SRR_MISSING" not in rows[0]
    assert rows[0]["SRR1"] == "200"


# --- what a run leaves behind ------------------------------------------------
def test_by_default_only_the_bam_survives():
    """By default only the BAM is kept; every other file is an intermediate that
    can be rebuilt from the accession."""
    cfg = cfgmod.load(None)
    keep = cfg["keep"]
    assert keep["bam"] is True
    assert [k for k, v in keep.items() if v] == ["bam"], f"kept by default: {keep}"


def test_resume_does_not_re_download_a_run_whose_bam_was_deleted_on_purpose(tmp_path):
    """With keep.bam off, a missing BAM does not make `resume` re-process a
    finished sample."""
    from ribomine.pipeline import is_processed
    from ribomine.utils import Sample, write_json

    s = Sample("SRR1", str(tmp_path))
    kept = cfgmod.load(None, {"project": {"workdir": str(tmp_path)}})
    dropped = cfgmod.load(None, {"project": {"workdir": str(tmp_path)},
                                 "keep": {"bam": False}})

    assert not is_processed(dropped, s), "nothing has run yet"
    write_json(s.process_json, {"run_accession": "SRR1", "bam_bytes": 123})

    assert is_processed(dropped, s), "a finished sample with no BAM is still done"
    assert not is_processed(kept, s), "but a BAM that was meant to be kept and is not there is a re-run"


# --- the streamed read sample ------------------------------------------------
def _fake_fastq_bytes(n: int) -> bytes:
    return b"".join(b"@r%d\nACGTACGTAC\n+\nIIIIIIIIII\n" % i for i in range(n))


class _BrokenStream:
    """What a dropped HTTPS connection looks like from inside gzip."""

    def readline(self):
        raise EOFError("Compressed file ended before the end-of-stream marker was reached")

    def close(self):
        pass


def test_a_dropped_read_sample_stream_is_retried_not_lost(tmp_path, monkeypatch):
    """A truncated stream is retried by re-opening it: curl writes to a pipe, so
    its own retry would corrupt the gzip stream."""
    import io

    from ribomine.sra import download

    monkeypatch.setattr(download.metadata, "fastq_urls",
                        lambda acc: ["https://ena.example/x.fastq.gz"])
    opened = []

    def fake_open(src):
        opened.append(src)
        if len(opened) == 1:
            return _BrokenStream(), None          # the connection drops
        return io.BytesIO(_fake_fastq_bytes(50)), None

    monkeypatch.setattr(download, "_open_stream", fake_open)
    out = str(tmp_path / "s.fastq")
    st = download.sample_reads("SRRFAKE", out, n=10, scan=1000, seed=1, backoff_s=0)

    assert len(opened) == 2, "the stream must be re-OPENED, not resumed"
    assert st["n_sampled"] == 10
    assert os.path.getsize(out) > 0


def test_a_corrupt_local_fastq_is_not_retried_four_times(tmp_path, monkeypatch):
    """A corrupt local file fails the same way every time, so it is not retried."""
    from ribomine.sra import download

    local = tmp_path / "reads.fastq"
    local.write_bytes(_fake_fastq_bytes(5))
    opened = []

    def fake_open(src):
        opened.append(src)
        return _BrokenStream(), None

    monkeypatch.setattr(download, "_open_stream", fake_open)
    with pytest.raises(RuntimeError, match="could not stream a read sample"):
        download.sample_reads(str(local), str(tmp_path / "o.fastq"), n=10, backoff_s=0)
    assert len(opened) == 1, "a local file is not a flaky network"


# --- bowtie2 thread cap ------------------------------------------------------
def test_bowtie2_is_capped_at_eight_threads_because_un_corrupts(monkeypatch, tmp_path):
    """bowtie2 is capped at 8 threads: with more, the records of its `--un` file
    interleave and the FASTQ is corrupt without any error."""
    from ribomine.qc import contaminants

    seen = {}

    class _Proc:
        stderr = "100 reads; of these:\n  40 (40.00%) aligned exactly 1 time\n"

    def fake_run(cmd, **kw):
        seen["cmd"] = [str(c) for c in cmd]
        return _Proc()

    monkeypatch.setattr(contaminants, "run", fake_run)
    monkeypatch.setattr(contaminants, "_bowtie2_bin", lambda: "bowtie2")

    contaminants._run_bowtie2("in.fq", "idx", str(tmp_path / "un.fq"), threads=24)
    cmd = seen["cmd"]
    assert "--un" in cmd, "the cap exists BECAUSE of --un; if it goes, revisit the cap"
    assert cmd[cmd.index("-p") + 1] == "8", f"24 threads must be capped to 8: {cmd}"

    # ... and a smaller request is honoured, not inflated to the cap
    contaminants._run_bowtie2("in.fq", "idx", str(tmp_path / "un.fq"), threads=4)
    cmd = seen["cmd"]
    assert cmd[cmd.index("-p") + 1] == "4"


def test_a_short_un_file_fails_the_sample_instead_of_being_mapped(monkeypatch, tmp_path):
    """bowtie2 read N and aligned M, so `--un` must hold N-M records; a mismatch
    means a corrupt file and fails the sample."""
    from ribomine.qc import contaminants

    idx = tmp_path / "idx.1.bt2"
    idx.write_bytes(b"x")
    (tmp_path / "in.fastq").write_text("@r\nACGT\n+\nIIII\n")

    cfg = cfgmod.load(None, {"project": {"workdir": str(tmp_path)},
                             "reference": {"contaminant_index": str(tmp_path / "idx")}})

    # bowtie2 claims 100 reads, 40 aligned -> --un must hold 60 ...
    monkeypatch.setattr(contaminants, "_run_bowtie2", lambda *a, **k: (100, 40))
    # ... but the screen only finds 58 records in it: two were lost in the interleave
    monkeypatch.setattr(contaminants, "_screen",
                        lambda *a, **k: {"n_in": 58, "n_low_complexity": 0, "n_kept": 58})

    with pytest.raises(RuntimeError, match="--un output is corrupt"):
        contaminants.filter_fastq(str(tmp_path / "in.fastq"), str(tmp_path / "o.fastq"),
                                  cfg, threads=8)

    # and matching counts pass
    monkeypatch.setattr(contaminants, "_screen",
                        lambda *a, **k: {"n_in": 60, "n_low_complexity": 5, "n_kept": 55})
    st = contaminants.filter_fastq(str(tmp_path / "in.fastq"), str(tmp_path / "o.fastq"),
                                   cfg, threads=8)
    assert st["n_input"] == 100 and st["n_kept"] == 55


def test_screen_drops_overlength_and_malformed_reads_that_crash_star(tmp_path):
    """_screen drops over-length and malformed (len(qual) != len(seq)) records,
    which make STAR fail or segfault, but still counts them in n_in so that
    filter_fastq's `--un` integrity check balances."""
    from ribomine.qc import contaminants

    in_fq = tmp_path / "in.fastq"
    out_fq = tmp_path / "out.fastq"
    long_seq = "ACGT" * 200                          # 800 nt -- STAR would overflow on it
    in_fq.write_text(
        "@good1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
        f"@toolong\n{long_seq}\n+\n{'I' * len(long_seq)}\n"   # well-formed but > MAX_READ_LEN
        "@malformed\nACGTACGTACGT\n+\nIIIIIIIIIII\n"          # 12 nt seq, 11 nt qual
        "@good2\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n"
    )
    s = contaminants._screen(str(in_fq), str(out_fq), min_entropy=1.1, max_base_frac=0.85)

    assert s["n_overlong"] == 1
    assert s["n_malformed"] == 1
    assert s["n_in"] == 4            # every record is read, so the --un check still holds
    assert s["n_kept"] == 2
    written = out_fq.read_text()
    assert "@toolong" not in written and "@malformed" not in written
    assert "@good1" in written and "@good2" in written


def test_bowtie2_abort_on_solid_colourspace_is_retried_on_base_space_reads(monkeypatch, tmp_path):
    """bowtie2 aborts (SIGABRT, "more quality values than read characters") on
    SOLiD colour-space reads; filter_fastq catches that abort, drops the reads
    that are not base-space, and retries."""
    from ribomine.qc import contaminants
    from ribomine.utils import ToolError

    (tmp_path / "idx.1.bt2").write_bytes(b"x")
    in_fq = tmp_path / "in.fastq"
    in_fq.write_text(
        "@good\nACGTACGT\n+\nIIIIIIII\n"
        "@solid\nT2220001\n+\nIIIIIIII\n"       # colour-space: equal length, but digit bases
    )
    cfg = cfgmod.load(None, {"project": {"workdir": str(tmp_path)},
                             "reference": {"contaminant_index": str(tmp_path / "idx")}})

    calls = {"n": 0}
    def fake_run_bowtie2(fastq, index, noncontam_fq, *, threads, log=""):
        calls["n"] += 1
        if calls["n"] == 1:
            raise ToolError(["bowtie2"], 134,
                            "Error: Read x has more quality values than read characters.\n"
                            "(ERR): bowtie2-align died with signal 6 (ABRT) (core dumped)")
        # retry: the input must be the base-space copy with the colour-space read gone
        txt = open(fastq).read()
        assert "@solid" not in txt and "@good" in txt
        with open(noncontam_fq, "w") as fh:
            fh.write("@good\nACGTACGT\n+\nIIIIIIII\n")   # survivor to --un, none aligned
        return (1, 0)

    monkeypatch.setattr(contaminants, "_run_bowtie2", fake_run_bowtie2)

    out_fq = tmp_path / "out.fastq"
    st = contaminants.filter_fastq(str(in_fq), str(out_fq), cfg, threads=8)

    assert calls["n"] == 2                       # aborted once, then retried
    assert st["n_malformed"] == 1                # the colour-space read is accounted for
    assert st["n_input"] == 2                    # and added back into the input total
    assert st["n_kept"] == 1
    assert "@good" in out_fq.read_text() and "@solid" not in out_fq.read_text()


def test_bowtie2_error_that_is_not_a_qual_mismatch_is_not_swallowed(monkeypatch, tmp_path):
    """Only the quality-length abort is retried; any other bowtie2 failure
    propagates."""
    from ribomine.qc import contaminants
    from ribomine.utils import ToolError

    (tmp_path / "idx.1.bt2").write_bytes(b"x")
    (tmp_path / "in.fastq").write_text("@r\nACGT\n+\nIIII\n")
    cfg = cfgmod.load(None, {"project": {"workdir": str(tmp_path)},
                             "reference": {"contaminant_index": str(tmp_path / "idx")}})

    def fake_run_bowtie2(*a, **k):
        raise ToolError(["bowtie2"], 1, "Error: could not open index files")
    monkeypatch.setattr(contaminants, "_run_bowtie2", fake_run_bowtie2)

    with pytest.raises(ToolError, match="could not open index"):
        contaminants.filter_fastq(str(tmp_path / "in.fastq"), str(tmp_path / "out.fastq"),
                                  cfg, threads=8)


# --- the pile-up filter's length-concentration cut ---------------------------
def _bam(tmp_path, name, reads):
    """reads = [(chrom, pos, length, n)] -> a tiny coordinate-sorted, indexed BAM."""
    import pysam

    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": "1", "LN": 1_000_000}]}
    path = str(tmp_path / name)
    recs = []
    for chrom, pos, ln, n in reads:
        for i in range(n):
            a = pysam.AlignedSegment()
            a.query_name = f"{chrom}_{pos}_{ln}_{i}"
            a.reference_id = 0
            a.reference_start = pos
            a.query_sequence = "A" * ln
            a.query_qualities = pysam.qualitystring_to_array("I" * ln)
            a.cigarstring = f"{ln}M"
            a.flag = 0
            recs.append(a)
    recs.sort(key=lambda a: a.reference_start)
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        for a in recs:
            out.write(a)
    pysam.index(path)
    return path


def test_a_miRNA_pile_is_removed_but_a_translated_codon_is_not(tmp_path):
    """A miRNA pile (one 5' base, ~22 nt, length concentration ~0.78) is removed;
    a translated codon, whose footprints spread over ~26-34 nt, is kept. The
    length-concentration cut has to sit between the two."""
    from ribomine.qc import pileups

    # a miRNA-like pile: 800 reads on one base, 78% of them exactly 22nt (lc = 0.78)
    mirna = [("1", 1000, 22, 624), ("1", 1000, 21, 100), ("1", 1000, 23, 76)]
    # a translated codon: just as abundant, but its footprints spread across lengths
    codon = [("1", 5000, 28, 240), ("1", 5000, 29, 200), ("1", 5000, 30, 180),
             ("1", 5000, 27, 120), ("1", 5000, 31, 60)]
    # background, so neither pile is a >=10% "dominant" position (a different rule)
    bg = [("1", 20000 + i * 7, 29, 30) for i in range(300)]

    inp = _bam(tmp_path, "in.bam", mirna + codon + bg)
    cfg = cfgmod.load(None, {"project": {"workdir": str(tmp_path)}})
    st = pileups.filter_bam(inp, str(tmp_path / "out.bam"), cfg, label="t")

    pos = {p["pos"] for p in st["top_pileups"]}
    assert 1000 in pos, "the miRNA pile (length concentration 0.78) must be removed"
    assert 5000 not in pos, "a translated codon's footprints spread over lengths -- keep it"
    assert st["n_reads_removed"] == 800


def test_every_mapping_summary_column_is_actually_produced(tmp_path):
    """PROCESS_COLUMNS and the keys of a row match exactly: write_tsv leaves an
    unfilled column blank and drops an unlisted key."""
    from ribomine import reports
    from ribomine.utils import Sample, write_json

    cfg = cfgmod.load(None, {"project": {"workdir": str(tmp_path)}})
    s = Sample("SRR1", cfg.workdir)
    # a process record with every block the row builder reads
    write_json(s.process_json, {
        "run_accession": "SRR1",
        "download": {"route": "ena_https", "mb_per_s": 40.0, "bytes": 10},
        "trim": {"n_reads_in": 100, "n_reads_out": 90, "mean_len_in": 50.0,
                 "mean_len_out": 30.0, "frac_no_adapter": 0.1},
        "contaminants": {"n_input": 90, "n_kept": 40, "n_contaminant_rRNA_tRNA_etc": 50},
        "mapping": {"n_input": 40, "n_unique": 30, "frac_unique": 0.75,
                    "frac_multimapping": 0.1, "frac_unmapped": 0.15,
                    "avg_input_len": 31.0, "avg_mapped_len": 29.4},
        "periodicity": {"n_reads_in_bam": 30, "n_reads_scored": 30, "mean_mapped_len": 29.4,
                        "mean_mapped_softclip": 0.3, "read_len_mode": 30,
                        "periodicity_inframe_frac": 0.6,
                        "periodicity_tvd_uniform": 0.4, "n_cds_reads": 20,
                        "cds_frac_of_genic": 0.8},
        "counts": {"n_in_genes": 25, "frac_in_genes": 0.83, "n_genes_detected": 9,
                   "n_ambiguous": 1, "n_no_feature": 4, "sense_over_antisense": 20.0},
        "umi_dedup": False, "bam_bytes": 123, "keep": {"bam": True},
    })
    row = reports._process_row(cfg, "SRR1")

    missing = [c for c in reports.PROCESS_COLUMNS if c not in row]
    extra = [k for k in row if k not in reports.PROCESS_COLUMNS]
    assert not missing, f"columns in the header that no row fills: {missing}"
    assert not extra, f"row keys write_tsv would silently drop: {extra}"

    # the leading columns, in order
    assert reports.PROCESS_COLUMNS[:8] == [
        "run_accession", "verdict", "architecture", "n_mapped", "mean_footprint_len",
        "mean_mapped_len", "mean_mapped_softclip", "periodicity_tvd"]
    assert row["n_mapped"] == 30
    assert row["mean_mapped_softclip"] == 0.3
    # the footprint length is that of STAR's input (trimmed, contaminant-free), not
    # mean_len_after_trim, which still includes the contaminants
    assert row["mean_footprint_len"] == 31.0 != row["mean_len_after_trim"]
    assert row["mean_mapped_len"] == 29.4


class _DeadCurl:
    """A curl that exited non-zero without writing anything (e.g. ENA answering
    403), which leaves an empty stream rather than a broken one."""

    returncode = 22

    def __init__(self):
        import io
        self.stdout = io.BytesIO(b"")
        self.stderr = io.BytesIO(b"curl: (22) The requested URL returned error: 403\n")

    def poll(self):
        return 22

    def terminate(self):
        pass

    def wait(self, timeout=None):
        return 22


def test_a_failed_transfer_is_retried_and_never_read_as_an_empty_run(tmp_path, monkeypatch):
    """An empty stream from a failed curl (e.g. HTTP 403) is retried: the transfer
    is judged by curl's exit status, not by whether bytes arrived."""
    import io

    from ribomine.sra import download

    monkeypatch.setattr(download.metadata, "fastq_urls",
                        lambda acc: ["https://ena.example/x.fastq.gz"])
    opened = []

    def fake_open(src):
        opened.append(src)
        if len(opened) == 1:
            p = _DeadCurl()
            return p.stdout, p                  # empty stream + a curl that failed
        return io.BytesIO(_fake_fastq_bytes(40)), None

    monkeypatch.setattr(download, "_open_stream", fake_open)
    st = download.sample_reads("SRRX", str(tmp_path / "s.fastq"), n=10, scan=100,
                               seed=1, backoff_s=0)
    assert len(opened) == 2, "a 403 must be retried, not read as a run with no reads"
    assert st["n_sampled"] == 10


def test_a_transfer_that_dies_MID_stream_is_not_kept_as_a_short_sample(tmp_path, monkeypatch):
    """If curl dies partway, the partial (and therefore biased) sample is
    discarded and the transfer retried."""
    import io

    from ribomine.sra import download

    class _HalfCurl(_DeadCurl):
        returncode = 56

        def __init__(self):
            self.stdout = io.BytesIO(_fake_fastq_bytes(3))   # a few reads, then death
            self.stderr = io.BytesIO(b"curl: (56) OpenSSL SSL_read: error\n")

        def poll(self):
            return 56

    monkeypatch.setattr(download.metadata, "fastq_urls",
                        lambda acc: ["https://ena.example/x.fastq.gz"])
    opened = []

    def fake_open(src):
        opened.append(src)
        if len(opened) == 1:
            p = _HalfCurl()
            return p.stdout, p
        return io.BytesIO(_fake_fastq_bytes(40)), None

    monkeypatch.setattr(download, "_open_stream", fake_open)
    st = download.sample_reads("SRRX", str(tmp_path / "s.fastq"), n=10, scan=100,
                               seed=1, backoff_s=0)
    assert len(opened) == 2, "the truncated draw must be thrown away, not sampled from"
    assert st["n_sampled"] == 10


# --- the QC figure -----------------------------------------------------------
def _qc_dict(reason: str) -> dict:
    """The smallest QC record plot_qc will draw."""
    return {
        "label": "SRRX", "verdict": "RIBO-SEQ", "is_riboseq": True, "verdict_reason": "",
        "reasons": ["footprint length ok (mode 29 nt)", reason],
        "n_reads_scored": 1000, "read_len_mode": 29, "read_len_peak_frac": 0.8,
        "periodicity_inframe_frac": 0.6, "periodicity_tvd_uniform": 0.4,
        "cds_frac_of_genic": 0.8, "start_codon_ratio": 5.0, "top5p_locus_frac": 0.01,
        "n_cds_reads": 500, "read_len_hist": {28: 300, 29: 500, 30: 200},
        "frame_by_len": {29: 0.6}, "frame_by_len_n": {29: 500},
        "region_frac": {"CDS": 0.8, "intron": 0.1, "intergenic": 0.1},
        "metagene_start": {-12: 100, 0: 50}, "metagene_stop": {-24: 40},
        "pooled_frame_counts": [300, 100, 100],
        "mapping": {"n_input": 1000, "frac_unique": 0.6, "frac_multimapping": 0.3,
                    "frac_unmapped": 0.1},
        "thresholds": {"footprint_len_lo": 25, "footprint_len_hi": 36,
                       "periodic_min": 0.4, "periodic_strong": 0.5},
    }


def test_a_long_reason_wraps_instead_of_stretching_the_whole_figure(tmp_path):
    """A long verdict reason wraps. Unwrapped, it widens the saved figure
    (savefig crops to the artists' bounding box) and squeezes the panels."""
    from PIL import Image

    from ribomine.qc import plot

    long_reason = (
        "MITOCHONDRIAL-DOMINATED (42% of reads) -- this looks like mitoribosome profiling. "
        "The verdict above was decided on the 1,338 NUCLEAR CDS reads; the 586 MT-CDS reads "
        "are 54% in-frame (TVD 0.21) and are not scored")
    assert len(long_reason) > 200

    a = str(tmp_path / "short.png")
    b = str(tmp_path / "long.png")
    plot.plot_qc(_qc_dict("unique mapping ok (18%)"), a, dpi=110)
    plot.plot_qc(_qc_dict(long_reason), b, dpi=110)

    wa = Image.open(a).size[0]
    wb = Image.open(b).size[0]
    assert wb <= wa + 20, (
        f"a long reason widened the figure from {wa}px to {wb}px -- it must wrap, "
        f"or the panels get squeezed and their titles overlap")


# --- strandedness / TI-seq ---------------------------------------------------
def test_a_read_antisense_to_a_CDS_is_counted_as_antisense_not_just_intron():
    """A read on a CDS but on the wrong strand is classified 'intron', yet `on_cds`
    is set, so the caller can recognise a reverse-complemented library."""
    import numpy as np
    from ncls import NCLS

    from ribomine.qc.verdict import _classify

    # one + strand CDS, 100..200, inside a gene
    idx = {"cds": {"1": {"start": np.array([100]), "end": np.array([200]),
                         "strand": np.array([1], dtype=np.int8),
                         "frame": np.array([0], dtype=np.int8)}}}
    one = NCLS(np.array([100], dtype=np.int64), np.array([200], dtype=np.int64),
               np.array([0], dtype=np.int64))
    ncls = {"cds": {"1": one}, "utr5": {}, "utr3": {}, "exon_nc": {},
            "gene": {"1": one}}

    region, frame, on_cds = _classify(idx, ncls, "1", 150, 1)      # sense read
    assert region == "CDS" and frame is not None and on_cds

    region, frame, on_cds = _classify(idx, ncls, "1", 150, -1)     # antisense read
    assert region == "intron", "an antisense CDS read still falls through to 'intron'"
    assert on_cds, "... but the caller must be able to SEE that a CDS was there"


def test_the_tiseq_cut_labels_an_initiation_dominated_run():
    """The TI-seq cut is 30: in a 100-run test cohort the initiation-dominated
    runs scored 34 and 37, and the highest elongating run 14.9."""
    cfg = cfgmod.load(None)
    cut = cfg["qc.tiseq_ratio_min"]
    assert cut == 30
    assert 34.1 >= cut, "SRR12790151 is TI-seq and must be called as such"
    assert 14.9 < cut, "... while the best ELONGATING run in that cohort must not be"

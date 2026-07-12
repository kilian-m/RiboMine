"""Tests for the parts where a silent wrong answer is possible.

The scientific core (profile / verdict / infer) is verified against its reference
implementation on real BAMs, which no unit test can substitute for. What is tested
here is the machinery that would corrupt a run *quietly*: the config typo guard,
the stage graph, the query's tiered text filter, and the trimming/UMI invariants.
"""
from __future__ import annotations

import json

import pytest

from ribomine import config as cfgmod
from ribomine.arch import infer, trim
from ribomine.config import ConfigError, stages_to_run
from ribomine.sra import query


# --- config ---------------------------------------------------------------
def test_unknown_key_is_an_error(tmp_path):
    """A typo'd threshold that is silently ignored is worse than a crash: the run
    completes, looks fine, and used the wrong number."""
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"qc": {"periodic_minn": 0.9}}))
    with pytest.raises(ConfigError, match="qc.periodic_minn"):
        cfgmod.load(str(p))


def test_defaults_are_complete_and_merge():
    cfg = cfgmod.load(None, {"qc": {"periodic_min": 0.9}})
    assert cfg["qc.periodic_min"] == 0.9          # override applied
    assert cfg["qc.tvd_min"] == 0.10              # sibling default survives
    assert cfg["process.umi_dedup"] is False      # dedup is OFF by default


def test_every_architecture_threshold_is_wired():
    """A config key nothing reads is a lie to the user. Thresholds must round-trip."""
    from dataclasses import fields
    names = {f.name for f in fields(infer.Thresholds)}
    assert names == set(cfgmod.DEFAULTS["architecture"]), (
        "architecture config keys and Thresholds fields have drifted apart")


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
    """The bug this guards: every Ribo-seq protocol ALSO depletes rRNA, so an
    exclude list that vetoes on 'rRNA depletion' throws away genuine Ribo-seq.
    Measured, a naive exclude silently dropped 21 real runs."""
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
    """ENA's text index is tokenised: a multi-word wildcard silently matches NOTHING.
    Fail loudly rather than return an empty archive."""
    cfg = cfgmod.load(None, {"query": {"terms": ["ribosome profiling"]}})
    with pytest.raises(ValueError, match="single words"):
        query._ena_query(cfg)


# --- trimming and the UMI invariant ---------------------------------------
CALL = {
    "status": "ok",
    "umi5_len": 2, "umi3_len": 5, "nt3_len": 0,
    "barcode3_seq": "AGCTA", "adapter3_name": "illumina_truseq",
    "adapter3_seq": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "polyA_tail": "none", "footprint_len_mode": 30, "p5_layout": [],
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


FOOT = "ACGTACGTACGTACGTACGTACGTACGTAC"          # 30 nt "footprint"
ADAP = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def test_umi_goes_to_the_header_and_the_rt_base_stays(tmp_path):
    read = "GG" + FOOT + "TTTTT" + "AGCTA" + ADAP    # umi5=GG, umi3=TTTTT, bc=AGCTA
    out = str(tmp_path / "o.fastq")
    st = trim.trim_fastq(_fastq(tmp_path, [read]), CALL, out, min_len=20)
    (name, seq), = _read(out)
    assert seq == FOOT                    # exactly the footprint, nothing else
    assert name.endswith("_GGTTTTT")      # 5' UMI + 3' UMI, in that order
    assert st["n_reads_out"] == 1


def test_every_umi_has_the_same_length(tmp_path):
    """umi_tools ABORTS on a variable-length UMI ('not all umis are the same
    length'). A read without the adapter never sequenced its 3' UMI, so it would
    carry a 2 nt UMI where every other read carries 7 -- the exact failure this
    pipeline hit on SRR12285169. Both policies must keep the length fixed."""
    with_adapter = "GG" + FOOT + "TTTTT" + "AGCTA" + ADAP
    no_adapter = "GG" + FOOT + "CCCCCCCCCCCCCCCC"        # insert ran off the read

    # default: drop the read that has no adapter
    out = str(tmp_path / "drop.fastq")
    st = trim.trim_fastq(_fastq(tmp_path, [with_adapter, no_adapter]), CALL, out)
    assert st["n_reads_out"] == 1 and st["n_dropped_untrimmed"] == 1
    assert {len(n.split("_")[-1]) for n, _ in _read(out)} == {7}

    # keep them: the unsequenced UMI bases are written as N, so the length holds
    out2 = str(tmp_path / "keep.fastq")
    st2 = trim.trim_fastq(_fastq(tmp_path, [with_adapter, no_adapter]), CALL, out2,
                          discard_untrimmed=False)
    assert st2["n_reads_out"] == 2 and st2["n_umi_padded"] == 1
    umis = [n.split("_")[-1] for n, _ in _read(out2)]
    assert {len(u) for u in umis} == {7}, "variable UMI length would break umi_tools"
    assert any(u.endswith("NNNNN") for u in umis)


def test_barcode_between_two_umi_blocks_is_not_swallowed(tmp_path):
    """iCLIP2-style [UMI][barcode][UMI][footprint]: flattening the layout to totals
    would take the barcode's bases into the UMI and corrupt the dedup key."""
    call = dict(CALL, umi5_len=4, p5_layout=[
        {"role": "umi5", "offset": 0, "len": 2},
        {"role": "barcode5", "offset": 2, "len": 3, "seq": "GAT"},
        {"role": "umi5", "offset": 5, "len": 2},
    ], functional={"trim_5p": 7, "dedup_umi_len": 9})
    read = "AC" + "GAT" + "TG" + FOOT + "TTTTT" + "AGCTA" + ADAP
    out = str(tmp_path / "o.fastq")
    trim.trim_fastq(_fastq(tmp_path, [read]), call, out)
    (name, seq), = _read(out)
    assert seq == FOOT
    assert name.endswith("_ACTGTTTTT")     # the two UMI blocks, NOT the barcode GAT


def test_infer_refuses_rather_than_guesses():
    """`undetermined` is a feature: a fabricated architecture would silently
    mis-trim every read in the dataset."""
    call = infer.infer({"label": "x", "n_used": 10})
    assert call["status"] == "undetermined"
    assert "10" in call["reason"]


# --- shipped reference data ------------------------------------------------
def test_contaminant_fasta_is_bundled():
    """The QC stage cannot silently run without a contaminant filter: an unfiltered
    rRNA read maps to hundreds of genomic copies, so the library then looks like
    ~78% multimapping junk. Shipping the reference means the default just works."""
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
    # which are excised during maturation and so appear in NO mature rRNA sequence --
    # without them those fragments reach the aligner and, rDNA being a high-copy repeat,
    # come back as multimappers. Worth +1.7pp of contaminant catch on an rRNA-heavy run.
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
    """Every BAM RiboMine leaves on disk is sorted+indexed. There is deliberately no
    config key that turns that off -- an unindexed BAM is one nobody can open."""
    assert "sort_index_bam" not in cfgmod.DEFAULTS["process"]


def test_null_means_is_documented_for_every_nullable_key():
    """In JSON a `null` reads as 'nothing / off'. For these keys it means the
    opposite -- 'work it out for me'. A reader cannot tell those apart from the file,
    so every nullable default must carry an explanation in the generated config."""
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
    """`ribomine query` searches the archive; it has no use for pipeline.start's
    inputs. Validating them would make the query unusable in exactly the workflow it
    exists for: query the archive, THEN feed the accessions back in."""
    cfg = cfgmod.load(None, {"pipeline": {"start": "accessions",
                                          "accession_list": "does/not/exist.txt"}})
    cfg.validate(need_reference=False, need_inputs=False)      # must not raise
    with pytest.raises(ConfigError, match="accession list not found"):
        cfg.validate(need_reference=False)                     # but `run` still checks


# --- the 3' anchor: a hidden adapter must not become a fabricated UMI ---------
def _t3_profile(t3_match, *, n_anchored=0, frac_anchored=0.0, anchor="none"):
    """A minimal profile exercising the no-adapter 3' fallback."""
    import numpy as np
    plateau = [0.99] * 24
    return {
        "label": "x", "n_used": 50_000, "n_anchored": n_anchored,
        "frac_anchored": frac_anchored, "anchor_kind": anchor, "anchor_offset": 30,
        "p5_match": plateau, "p5_comp": [[0.25] * 4] * 24,
        "t3_match": t3_match, "t3_comp": [[0.25] * 4] * len(t3_match),
        "adap_match": [float("nan")] * 60, "adap_comp": [[float("nan")] * 4] * 60,
        "footprint_len_hist": {str(L): 100 for L in range(26, 35)},
        "read_len_hist": {"46": 900}, "top5p_locus_frac": 0.01,
    }


def test_a_flat_3p_tail_is_a_umi():
    """A genuinely trimmed deposit with a retained UMI: the walked positions sit AT
    chance for all of them, then jump. That is a step, and it is a real 5-nt UMI."""
    flat = [0.26, 0.25, 0.26, 0.25, 0.26] + [0.99] * 12      # 5 non-genomic, then plateau
    call = infer.call_3prime(_t3_profile(flat), 0.99, [], infer.Thresholds())
    res, deposit = call
    assert res["umi3_len"] == 5
    assert deposit == "adapter_trimmed_umi_retained"


def test_a_ramping_3p_tail_is_refused_not_read_as_a_long_umi():
    """SRR11945406: a McGlincy-Ingolia library whose 46-nt reads barely reach the
    adapter, so the adapter anchor was rejected and the caller fell back to the read's
    own 3' end. There the construct+adapter SMEAR across positions and the match rate
    RAMPS instead of stepping. A level-only walk read that ramp as one 11-nt UMI --
    which was really 5 nt of UMI + a 5-nt barcode + an adapter base. Refuse the ramp."""
    ramp = [0.27, 0.27, 0.27, 0.26, 0.26, 0.26, 0.28, 0.31, 0.33, 0.39,
            0.52, 0.72, 0.86, 0.93, 0.96, 0.99]
    flags: list[str] = []
    res, deposit = infer.call_3prime(_t3_profile(ramp), 0.99, flags, infer.Thresholds())
    assert deposit == "unknown", "a ramp must not be read as a fixed-length construct"
    assert any(f.startswith("3p_tail_ramps_not_steps") for f in flags)
    assert res["umi3_len"] == "unknown"


def test_a_minority_adapter_is_still_an_adapter():
    """When the insert is as long as the read, only the short-footprint minority
    reaches the adapter. Enough of those reads (and a construct that sits a FIXED
    distance from the footprint end) make it a real anchor -- rejecting it is what
    forced the fallback that fabricated the 11-nt UMI."""
    from ribomine.qc import profile as prof
    assert prof.MIN_ANCHOR_READS <= 500
    # the caller must accept an anchor carried by few reads but many of them
    thr = infer.Thresholds(min_anchor_frac=0.15)
    p = _t3_profile([0.99] * 16, n_anchored=742, frac_anchored=0.023, anchor="panel:x")
    weak = (p["frac_anchored"] < thr.min_anchor_frac
            and p["n_anchored"] < infer.MIN_ANCHORED_READS)
    assert not weak, "742 anchored reads is ample evidence, whatever the fraction"

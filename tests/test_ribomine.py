"""Tests for the parts where a silent wrong answer is possible.

The scientific core (profile / verdict / infer) is verified against its reference
implementation on real BAMs, which no unit test can substitute for. What is tested
here is the machinery that would corrupt a run *quietly*: the config typo guard,
the stage graph, the query's tiered text filter, and the trimming/UMI invariants.
"""
from __future__ import annotations

import json
import os

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


def test_the_default_dedup_tool_can_run_the_default_dedup_method():
    """Not every method exists in every backend -- `percentile` is umi_tools-only.
    Shipping a default pair that cannot run would fail only for the users who turn
    dedup on, i.e. late, on someone else's cohort."""
    from ribomine.process.dedup import TOOLS, UMICOLLAPSE_ALGO

    cfg = cfgmod.load(None)
    tool, method = cfg["process.umi_dedup_tool"], cfg["process.umi_dedup_method"]
    assert tool in TOOLS
    if tool == "umicollapse":
        assert method in UMICOLLAPSE_ALGO


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


# --- the 3' scaffold: the barcode is part of the anchor, not just something to cut ---
BC_CALL = dict(CALL, umi3_len=5, barcode3_seq="ATCGT", footprint_len_mode=32,
               functional={"trim_5p": 2, "dedup_umi_len": 7})


def test_barcode_is_trimmed_and_kept_out_of_the_umi():
    """A sample barcode is CONSTANT across the library -- measured on SRR11945406 it is
    98.6% one sequence, 0.14 bits of entropy against a 5-nt UMI's ~10. Putting it in the
    dedup key adds no information at all, while adding error-prone positions that eat
    umi_tools --directional's edit-distance-1 budget. So: trim it, never dedup on it."""
    fn = infer.functional_view(
        {"umi5_len": 2, "barcode5_seq": "none", "p5_layout": [], "rt5_len": 0,
         "rt5_penetrance": 0.0, "rt_nt": False, "ts5_len": 0, "ts5_seq": ""},
        {"umi3_len": 5, "nt3_len": 0, "barcode3_seq": "ATCGT", "polyA_tail": "none",
         "adapter3_name": "illumina_truseq", "adapter3_seq": ADAP})
    bc = [s for s in fn["segments_3p"] if s["role"] == "barcode3"][0]
    assert bc["cat"] == "fixed-templated" and bc["fate"] == "trim"   # trimmed...
    assert fn["dedup_umi_len"] == 7                                  # ...but 2+5, not 2+5+5
    assert "ATCGT" not in fn["dedup_umi_mask"]
    assert fn["trim_3p_construct"] == 10                             # umi + barcode both cut


def test_the_scaffold_anchor_recovers_reads_the_adapter_alone_cannot(tmp_path):
    """SRR11945406: a 44-nt molecule in a 46-nt read. The adapter runs off the end of
    97% of reads, so anchoring on it alone found nothing and discard_untrimmed threw the
    library away (2.3% survived). The BARCODE is just as fixed and just as known, and it
    sits 5 nt closer to the insert -- anchoring on [barcode + adapter] recovers them."""
    # a read whose TruSeq is truncated to 2 nt: adapter-only (min_overlap 7) cannot see it
    read = "GG" + FOOT[:32] + "ACGTA" + "ATCGT" + ADAP[:2]
    out = str(tmp_path / "o.fastq")
    st = trim.trim_fastq(_fastq(tmp_path, [read]), BC_CALL, out, min_len=20)
    (name, seq), = _read(out)
    assert seq == FOOT[:32]              # the footprint, exactly
    assert name.endswith("_GGACGTA")     # 2 nt 5' UMI + 5 nt 3' UMI; barcode NOT in it
    assert st["n_reads_out"] == 1


def test_scaffold_anchoring_is_a_strict_superset(tmp_path):
    """A sequencing error in the barcode blocks the scaffold match (at a 7-nt overlap the
    error budget is zero). Those reads were trimmable before, so the adapter alone must
    still be tried -- the new anchor may only ever find MORE reads, never fewer."""
    bad_bc = "GG" + FOOT[:32] + "ACGTA" + "ATCGA" + ADAP     # barcode ATCGT -> ATCGA
    out = str(tmp_path / "o.fastq")
    st = trim.trim_fastq(_fastq(tmp_path, [bad_bc]), BC_CALL, out, min_len=20)
    assert st["n_reads_out"] == 1, "a barcode typo must not lose a read the adapter can anchor"
    (_, seq), = _read(out)
    assert seq == FOOT[:32]


def test_losing_most_of_a_library_is_reported(tmp_path):
    """Silently keeping 2% of a dataset is the worst possible outcome: it looks like it
    worked. The fraction must reach the caller and the TSV."""
    reads = [ "GG" + FOOT[:32] + "ACGTA" + "ATCGT" + ADAP ] + \
            [ "GG" + FOOT + "ACGTACGTAC" ] * 9              # 9 reads with no scaffold at all
    st = trim.trim_fastq(_fastq(tmp_path, reads), BC_CALL, str(tmp_path / "o.fq"))
    assert st["frac_no_adapter"] == 0.9
    assert st["frac_kept"] == 0.1


def test_polyA_tail_is_trimmed_even_when_an_adapter_is_also_present(tmp_path):
    """[footprint][poly-A][adapter]: the tail is enzymatically added, variable-length and
    NOT genomic, and the functional view already says fate=trim. But the adapter branch
    used to cut only the FIXED construct (here 0 nt), leaving the whole tail on the
    footprint -- on SRR19641906 that left 99.6% of trimmed reads ending in a run of A's,
    which end-to-end alignment then has to explain. Both must come off."""
    # NB the footprint must not itself end in A: poly(A) trimming cannot tell a genomic
    # terminal A from the tail, and does not try to (nor does cutadapt).
    fp = FOOT[:28]
    assert not fp.endswith("A")
    call = dict(CALL, umi5_len=0, umi3_len=0, barcode3_seq="none", polyA_tail="polyA",
                footprint_len_mode=28, functional={"trim_5p": 0, "dedup_umi_len": 0})
    out = str(tmp_path / "o.fastq")
    trim.trim_fastq(_fastq(tmp_path, [fp + "AAAAAAAAAAAAAA" + ADAP]), call, out, min_len=20)
    (_, seq), = _read(out)
    assert seq == fp, "the poly(A) tail must not survive as footprint"


def test_polyA_still_trimmed_when_the_adapter_is_beyond_it(tmp_path):
    """The other poly(A) case -- adapter not visible past the tail -- must keep working."""
    fp = FOOT[:28]
    call = dict(CALL, umi5_len=0, umi3_len=0, barcode3_seq="none", polyA_tail="polyA",
                adapter3_name="none_visible", adapter3_seq="none_visible",
                footprint_len_mode=28, functional={"trim_5p": 0, "dedup_umi_len": 0})
    out = str(tmp_path / "o.fastq")
    trim.trim_fastq(_fastq(tmp_path, [fp + "AAAAAAAAAAAA"]), call, out, min_len=20)
    (_, seq), = _read(out)
    assert seq == fp


def test_a_polyA_tail_anchors_the_cut_when_the_adapter_is_gone(tmp_path):
    """A poly(A) tail is SELF-ANCHORING: it sits between the footprint and the adapter,
    so it marks the footprint's 3' end whether or not the adapter made it into the read.
    Discarding those reads instead (SRR19641906: 28% retention) threw away exactly the
    long-footprint reads, which is also what biased the surviving length distribution."""
    fp = FOOT[:28]
    call = dict(CALL, umi5_len=0, umi3_len=0, barcode3_seq="none", polyA_tail="polyA",
                footprint_len_mode=28, functional={"trim_5p": 0, "dedup_umi_len": 0})
    # the adapter is present but only 2 nt of it -- far too short to anchor on (min 7).
    # TruSeq begins with an A, so the read does not even END in an A-run.
    read = fp + "AAAAAAAAAAAA" + ADAP[:2]
    assert not read.endswith("A")
    out = str(tmp_path / "o.fastq")
    st = trim.trim_fastq(_fastq(tmp_path, [read]), call, out, min_len=20)
    assert st["n_reads_out"] == 1, "the tail locates the boundary; the read must survive"
    (_, seq), = _read(out)
    assert seq == fp

    # and with NO adapter fragment at all
    st2 = trim.trim_fastq(_fastq(tmp_path, [fp + "AAAAAAAAAAAA"]), call,
                          str(tmp_path / "p.fastq"), min_len=20)
    assert st2["n_reads_out"] == 1


def test_an_adapter_followed_by_more_sequence_is_still_found(tmp_path):
    """The adapter is not always the last thing in the read. Sequence past it -- an
    index, a second adapter, a sample barcode -- is normal, and the read must still be
    cut AT the adapter.

    Without EndSkip.QUERY_STOP the aligner demands the alignment reach the END of the
    read, so an adapter with anything after it matches in ZERO reads. Measured on
    SRR25706716 ([footprint][Ingolia linker][12nt][TruSeq]): the linker is an exact
    substring of 98% of its reads and the matcher found it in none of them; with
    discard_untrimmed on, a 53.7M-read library became a 23 KB BAM. Five of the eight
    libraries in a random 20-run cohort were hit."""
    linker = "CTGTAGGCACCATCAAT"
    truseq = "AGATCGGAAGAGCACACGTCTGAACT"
    fp = "CGGGACATGTGGCGTACGAA"                      # 20 nt "footprint"
    read = fp + linker + "GGCCGGTTTCTG" + truseq     # adapter, then 38 nt more

    # the matcher must locate it where it actually is
    assert read.find(linker) == len(fp)
    assert trim.find_adapter(read, linker, min_start=14, min_overlap=7) == len(fp)

    call = {"status": "ok", "umi5_len": 0, "umi3_len": 0, "nt3_len": 0,
            "barcode3_seq": "none", "adapter3_name": "ingolia_linker",
            "adapter3_seq": linker, "polyA_tail": "none", "footprint_len_mode": 20,
            "p5_layout": [], "functional": {"trim_5p": 0, "dedup_umi_len": 0}}
    out = str(tmp_path / "o.fastq")
    st = trim.trim_fastq(_fastq(tmp_path, [read]), call, out, min_len=15)

    assert st["n_reads_out"] == 1, "the read carries its adapter; it must not be discarded"
    assert st["frac_no_adapter"] == 0.0
    (_, seq), = _read(out)
    assert seq == fp, "everything from the adapter onwards comes off, not just the adapter"


def test_an_adapter_at_the_very_end_still_works(tmp_path):
    """The ordinary case must not regress: adapter runs to the read's end, or past it."""
    fp = "ACGTACGTACGTACGTACGTACGTACGTAC"
    for tail in (ADAP, ADAP[:9]):        # complete, and truncated by the read end
        out = str(tmp_path / f"o{len(tail)}.fastq")
        call = dict(CALL, umi5_len=0, umi3_len=0, barcode3_seq="none",
                    footprint_len_mode=30, functional={"trim_5p": 0, "dedup_umi_len": 0})
        st = trim.trim_fastq(_fastq(tmp_path, [fp + tail]), call, out, min_len=20)
        assert st["n_reads_out"] == 1
        (_, seq), = _read(out)
        assert seq == fp


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
    """A truncated stream permanently failed the run before this: measured at 2 of 20
    runs when four samples stream at once. A tenth of a mining cohort is not an
    acceptable price for a transient, and curl cannot retry it -- it is writing to a
    pipe, so its retry would splice the head of the file into the middle of the gzip
    stream rather than recover anything."""
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
    """Re-reading a corrupt file on disk fails identically every time. Retrying it is
    pure latency, and it hides the fact that the file -- not the network -- is broken."""
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


# --- gene counts / the read-count matrix ------------------------------------
READS_PER_GENE = "\n".join([
    # STAR's four bookkeeping rows, then the genes.
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
    """A ribosome footprint is a piece of the mRNA, so it maps to the transcript's own
    strand. Reading the unstranded column would fold in antisense background; reading
    the antisense column would report the background INSTEAD of the library."""
    from ribomine.process import counts

    c, stats = counts.read_counts(_counts_tab(tmp_path, "r.tab"))
    assert c == {"ENSG01": 200, "ENSG02": 0, "ENSG03": 90}   # sense column, not 210/95
    assert stats["n_in_genes"] == 290
    assert stats["n_genes_detected"] == 2                    # the zero gene is not "detected"
    # what the counts do NOT contain -- the honest denominator
    assert stats["n_no_feature"] == 35 and stats["n_ambiguous"] == 8
    assert stats["frac_in_genes"] == round(290 / (290 + 35 + 8), 4)
    # STAR's N_multimapping row is 0 whenever multimap_nmax=1 (multimappers are dropped
    # before counting), so it is deliberately not reported -- Log.final.out has the truth
    assert "n_multimapping" not in stats
    assert stats["sense_over_antisense"] == round(290 / 15, 1)


def test_a_library_that_is_not_sense_stranded_is_called_out(tmp_path, caplog):
    """If the reads are on the other strand, the sense column is a fraction of the
    library rather than a measurement of it -- and every count in the matrix is wrong
    by that factor. It must not pass silently."""
    from ribomine.process import counts

    flipped = "\n".join([
        "N_noFeature\t0\t0\t0",
        "N_ambiguous\t0\t0\t0",
        "ENSG01\t2000\t1000\t1000",     # 1:1 -- unstranded or reversed, not ribo-seq
    ]) + "\n"
    with caplog.at_level("WARNING"):
        _, stats = counts.read_counts(_counts_tab(tmp_path, "f.tab", flipped), label="SRRX")
    assert stats["sense_over_antisense"] == 1.0
    assert "does not look sense-stranded" in caplog.text


def test_the_matrix_has_a_row_for_every_gene_including_the_zero_ones(tmp_path):
    """A matrix whose row set depends on which runs are in it cannot be compared with
    the next one. Rows come from the annotation, not from the data."""
    from ribomine.process import counts
    from ribomine.utils import read_tsv

    out = str(tmp_path / "m.tsv")
    counts.matrix(out, _star_index(tmp_path),
                  [("SRR1", _counts_tab(tmp_path, "a.tab")),
                   ("SRR2", _counts_tab(tmp_path, "b.tab"))])
    rows = read_tsv(out)
    assert [r["gene_id"] for r in rows] == ["ENSG01", "ENSG02", "ENSG03"]
    assert rows[0]["gene_name"] == "AAA"          # names come free from the STAR index
    assert rows[0]["SRR1"] == "200" and rows[0]["SRR2"] == "200"
    assert rows[1]["SRR1"] == "0"                 # a gene with no reads is a 0, not a gap


def test_a_run_with_no_counts_is_left_out_rather_than_left_blank(tmp_path):
    """A blank is not a zero. Every downstream tool reads this file as a numeric table,
    so a run that produced no counts must not become a column of empty cells."""
    from ribomine.process import counts
    from ribomine.utils import read_tsv

    out = str(tmp_path / "m.tsv")
    counts.matrix(out, _star_index(tmp_path),
                  [("SRR1", _counts_tab(tmp_path, "a.tab")),
                   ("SRR_MISSING", str(tmp_path / "nope.tab"))])
    rows = read_tsv(out)
    assert "SRR_MISSING" not in rows[0]
    assert rows[0]["SRR1"] == "200"


def test_every_mapping_summary_column_is_actually_produced(tmp_path):
    """A column in the header that no row ever fills is a blank column, and a key a row
    fills that the header does not list is silently DROPPED by write_tsv. Either way the
    table lies about what was measured, so the two lists must match exactly."""
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
                        "read_len_mode": 30, "periodicity_inframe_frac": 0.6,
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

    # the headline columns the table leads with, in order
    # columns 4 and 5 (0-based): the FOOTPRINTS, then the MAPPINGS
    assert reports.PROCESS_COLUMNS[:7] == [
        "run_accession", "verdict", "architecture", "n_mapped", "mean_footprint_len",
        "mean_mapped_len", "periodicity_tvd"]
    assert row["n_mapped"] == 30
    # the footprint length is STAR's INPUT (trimmed + contaminant-free), never
    # mean_len_after_trim, which still has the contaminants in it
    assert row["mean_footprint_len"] == 31.0 != row["mean_len_after_trim"]
    assert row["mean_mapped_len"] == 29.4

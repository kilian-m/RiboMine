"""Call the read architecture from a positional profile.

Reading the profile
-------------------
`p5_match[p]` is the fraction of reads whose base p matches the genome base the
alignment implies for it. Because a soft-clipped base still has an implied
coordinate, this is close to

    P(p is inside the aligned block) + P(p is clipped) * 1/4

i.e. essentially the CDF of the 5'-clip length, lifted by a 1/4 chance-match
floor. Three consequences drive the rules below:

* A 5' UMI of length u forces a clip of >= u in nearly every read, so
  `p5_match[p] ~ 0.25` for p < u and jumps to ~0.95 at p = u. A **step**.
* Low-quality or mis-mapped reads produce a long tail of large clips, which
  lifts `p5_match` **smoothly** across many positions. Not a step. Calling a
  UMI on the level alone would mistake that tail for architecture, so the RT
  test below looks at the increment, not the level.
* A non-templated base added by the reverse transcriptase is present in only a
  fraction f of molecules, so it depresses exactly one position by
  f * (plateau - 0.25) and creates a single sharp increment.

On the 3' side the footprint length varies read to read, so read position is
meaningless. The alignment end is also useless as an anchor: local alignment
stops at a mismatch, so the base just past it disagrees with the genome *by
construction* (match rate 0.00). The two honest anchors are the adapter start
and the read's own 3' end.

Segment vocabulary, insert-first:
    [5' UMI][RT nt][=== footprint ===][3' UMI][barcode][adapter]
Any of them may be absent, and the 3' UMI and barcode can swap places.

`infer()` never guesses: when the profile carries no readable boundary it
returns `status='undetermined'` with a human-readable reason (see §9 of
PIPELINE.md). That refusal is an answer, not an error.
"""
from __future__ import annotations

import logging
import math
from dataclasses import asdict, dataclass, fields
from typing import Any

import numpy as np

from ..config import DEFAULTS, Config

LOG = logging.getLogger("ribomine.arch.infer")


# --- thresholds ------------------------------------------------------------
# Match rates are interpreted relative to this sample's own genomic plateau.
# Mis-mapped and low-quality reads match the genome at the 1/4 chance rate at
# every position, so they pull the whole profile down by a constant factor. A
# noisy library can plateau at 0.79 where a clean one reaches 0.97; an absolute
# cut-off would call the noisy library's footprint "non-genomic".
_ARCH: dict[str, Any] = DEFAULTS["architecture"]


@dataclass(frozen=True)
class Thresholds:
    """One field per key of the config's `architecture` section.

    These are the knobs a user tunes per organism / per cohort; the defaults are
    the values tuned on a 187-sample human cohort. The deep internals below
    (CHANCE, MAX_WALK, the adapter PANEL, ...) are properties of the method
    rather than of the data, and stay module constants.
    """

    min_reads: int = _ARCH["min_reads"]
    # below this deep match rate the alignments are too noisy to read
    min_plateau: float = _ARCH["min_plateau"]
    # footprint length more concentrated than this = fixed-length artefact
    footprint_uniform_max: float = _ARCH["footprint_uniform_max"]
    # more than half the reads at one 5' coordinate = a single over-represented
    # species (adapter dimer / contaminant), not genuine footprints (the busiest
    # genuine library in the cohort: 0.38)
    single_locus_max: float = _ARCH["single_locus_max"]
    genomic_frac: float = _ARCH["genomic_frac"]      # fraction of the way from chance to plateau => genomic
    umi_match_frac: float = _ARCH["umi_match_frac"]  # a 3' UMI/barcode base matches the genome no more than this
                                                     # far above chance; higher = the footprint's own 3' edge
    ent_random: float = _ARCH["ent_random"]          # bits; above this a position is a uniform random position
    ent_const: float = _ARCH["ent_const"]            # bits; below this a position is a fixed sequence
                                                     # in between: a degenerate randomer (e.g. C/G only)
    footprint_pen: float = _ARCH["footprint_pen"]    # penetrance below this (sustained) => genomic footprint
    struct_pen: float = _ARCH["struct_pen"]          # penetrance above this => structural, in ~every read
    rt_max_len: int = _ARCH["rt_max_len"]            # a non-templated RT addition is 1-2 nt, never a UMI
    rt_pen_hi: float = _ARCH["rt_pen_hi"]            # a footprint-adjacent position below this penetrance is a
                                                     # variable-length (partial) RT base, not a fixed UMI position.
                                                     # Set below the ~0.72-0.75 to which leakage (chance genomic
                                                     # matches extending the alignment) depresses a real last UMI
                                                     # base, so those are not peeled; biased RT bases above it are
                                                     # still caught by the A/T-composition test.
    at_jump: float = _ARCH["at_jump"]                # A+T share rises this much over the UMI baseline at an RT base
    at_abs: float = _ARCH["at_abs"]                  # ... or reaches this absolute A+T share
    gc_template: float = _ARCH["gc_template"]        # G+C share marking a template-switch / TSO / linker addition
    gc_mid: float = _ARCH["gc_mid"]                  # a moderate-G/C position bridges into a strong-G/C run
    linker_max: int = _ARCH["linker_max"]            # longest template-switch / linker run peeled off the footprint
    min_anchor_frac: float = _ARCH["min_anchor_frac"]  # a short read may only expose the adapter in a minority
    # not used by the caller itself: the pipeline's policy for what to do with a
    # sample the caller refuses (carried here so the section maps 1:1 onto a field)
    process_undetermined: bool = _ARCH["process_undetermined"]

    @classmethod
    def from_config(cls, cfg: Config) -> "Thresholds":
        section = cfg.get("architecture", {}) or {}
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in section.items() if k in known})


# --- deep internals: properties of the method, not knobs -------------------
CHANCE = 0.25            # a non-genomic base matches the genome 1/4 of the time
BOUNDARY_STEP = 0.15     # the footprint boundary is a step, not a slope
HOMOPOLYMER_FRAC = 0.85  # one base at this share, repeated, is a tail not a barcode
HOMOPOLYMER_MIN = 4
RT_MIN_STEP = 0.08       # minimum jump in p5_match attributable to an RT base
RT_STEP_RATIO = 3.0      # ... and it must exceed this multiple of the local drift
RT_MIN_PENETRANCE = 0.12
RT_TERMINAL_PEN_LO = 0.15  # a 5'-terminal RT base (no UMI in front) sits at
                         # intermediate penetrance; below this position 0 is
                         # genomic (a footprint 5'-ligation bias), not an RT base
MAX_WALK = 26            # positions to walk left of the adapter
MIN_ANCHORED_READS = 300  # an adapter this well evidenced is real even in a minority
                          # of reads (a long insert hides it from the rest)
TAIL_FLAT_MAX = 0.10     # a fixed-length 3' construct is FLAT at chance across all its
                         # positions; a match rate that ramps is a smear of constructs
                         # of differing length, not one construct -- refuse to read a
                         # UMI length off it (see call_3prime)
MAX_BARCODE = 8          # constant blocks longer than this are linker, not sample barcode

PANEL = {
    "illumina_truseq": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "illumina_smallrna_ra3": "TGGAATTCTCGGGTGCCAAGG",
    "ingolia_linker": "CTGTAGGCACCATCAAT",
    "nextera": "CTGTCTCTTATACACATCT",
    "nebnext_smallrna": "AGATCGGAAGAGCACACGTCT",
}
# probes used to find a known adapter inside a de-novo constant block
PANEL_CORES = {name: seq[:10] for name, seq in PANEL.items()}

# short display names for the 3' adapter in the architecture string
ADAPTER_DISPLAY = {
    "illumina_truseq": "TruSeq",
    "illumina_smallrna_ra3": "smallRNA-RA3",
    "ingolia_linker": "Ingolia-linker",
    "nextera": "Nextera",
    "nebnext_smallrna": "NEBNext",
    "other": "custom-adapter",
    "none": "no-adapter",
    "unknown": "adapter?",
    "none_visible": "beyond poly(A)?",
}

# the profile fields the caller cannot work without
_REQUIRED = ("p5_match", "p5_comp", "t3_match", "t3_comp", "adap_match", "adap_comp",
             "anchor_kind", "anchor_offset", "n_anchored", "frac_anchored",
             "footprint_len_hist", "read_len_hist")


def adapter_disp(name) -> str:
    return ADAPTER_DISPLAY.get(str(name), str(name))


# --- small numeric helpers -------------------------------------------------
def entropy(row) -> float:
    r = np.asarray(row, dtype=float)
    if np.isnan(r).any() or r.sum() <= 0:
        return float("nan")
    r = r[r > 0]
    return float(-(r * np.log2(r)).sum())


def consensus(row, bases: str = "ACGT") -> str:
    r = np.asarray(row, dtype=float)
    if np.isnan(r).any():
        return "N"
    return bases[int(r.argmax())]


def dominant(row, bases: str = "ACGT") -> tuple[str, float]:
    r = np.asarray(row, dtype=float)
    if np.isnan(r).any():
        return "N", 0.0
    i = int(r.argmax())
    return bases[i], float(r[i])


def mode_of(hist: dict | None) -> int | None:
    if not hist:
        return None
    return int(max(hist.items(), key=lambda kv: kv[1])[0])


def gc_fraction(row) -> float:
    r = np.asarray(row, dtype=float)
    return float("nan") if np.isnan(r).any() else float(r[1] + r[2])


def at_fraction(row) -> float:
    """A+T share. The reverse transcriptase's non-templated base is mostly A;
    minus-strand reads are reverse-complemented to a common orientation, so it
    reads as T. Either way the RT position is A+T-enriched over a uniform UMI."""
    r = np.asarray(row, dtype=float)
    return float("nan") if np.isnan(r).any() else float(r[0] + r[3])


def classify(ent: float, thr: Thresholds) -> str:
    """A structural position is either a fixed base or a randomised one.

    Real randomers are often skewed -- one library's 5' randomer is drawn from
    C/G only, giving 1.1 bits rather than 2.0. That is still a randomer, not an
    ambiguity, so anything above the constant threshold counts as random and the
    skew is merely reported.
    """
    if math.isnan(ent):
        return "?"
    if ent <= thr.ent_const:
        return "const"
    if ent >= thr.ent_random:
        return "random"
    return "degenerate"


def segment(kinds, seqs) -> list[tuple[str, str]]:
    """Split a run of position kinds into maximal same-kind blocks."""
    blocks: list[list] = []
    for k, b in zip(kinds, seqs):
        if blocks and blocks[-1][0] == k:
            blocks[-1][1] += b
        else:
            blocks.append([k, b])
    return [(k, s) for k, s in blocks]


def genomic_plateau(profile: dict) -> float | None:
    """What a genuinely genomic base scores in this sample.

    Deep read positions are inside the footprint for essentially every read, so
    their match rate measures the sample's own ceiling -- error rate, junk reads
    and all. Everything else is judged relative to it.
    """
    m5 = np.asarray(profile["p5_match"], dtype=float)
    good = m5[~np.isnan(m5)]
    if len(good) < 10:
        return None
    return float(np.median(good[-6:]))


# --- the 5' side -----------------------------------------------------------
def call_5prime(profile: dict, plateau: float, flags: list[str], thr: Thresholds) -> dict:
    """Segment the 5' construct into UMI / barcode / non-templated additions.

    Model, insert-first:
        [non-templated (TSO G-run / RT-added C, G/C-rich, full penetrance)]?
        [UMI (uniform, balanced G/C, full penetrance)]?
        [barcode (fixed base, full penetrance)]?
        [template-switch / linker (G/C-rich, full penetrance)]?
        [RT base (A/T-enriched, partial penetrance)]?
        [footprint (genomic)]

    Every non-UMI element is a *non-templated addition* distinguished from the
    UMI by one composition axis: template-switch/TSO/linker positions are
    G/C-biased, the RT base is A/T-biased, a barcode is a fixed base -- where a
    UMI is balanced and uniform. Penetrance sets the frame: a UMI/barcode/linker
    is present in every read (penetrance ~1); the RT base in only a fraction
    (intermediate); the footprint never (~0).
    """
    m5 = np.asarray(profile["p5_match"], dtype=float)
    c5 = profile["p5_comp"]
    n = len(m5)
    pen = np.clip((plateau - m5) / max(plateau - CHANCE, 1e-6), 0.0, 1.0)
    ent = [entropy(c5[p]) for p in range(n)]
    con = [consensus(c5[p]) for p in range(n)]
    gc = [gc_fraction(c5[p]) for p in range(n)]
    at = [at_fraction(c5[p]) for p in range(n)]

    # --- footprint start: the first position that is, and stays, genomic. Some
    #     poor libraries only reach penetrance ~0.4 in the footprint, so "genomic"
    #     is judged relative to struct_pen, not an absolute ~0.
    f = None
    for p in range(n):
        if math.isnan(m5[p]):
            break
        if pen[p] < thr.footprint_pen and (p + 1 >= n or math.isnan(m5[p + 1])
                                           or pen[p + 1] < thr.footprint_pen):
            f = p
            break
    if f is None:
        flags.append("5p_never_reaches_genomic")
        return {"umi5_len": 0, "barcode5_seq": "unknown", "p5_layout": [],
                "rt5_len": 0, "rt5_penetrance": 0.0, "rt_nt": False,
                "ts5_len": 0, "ts5_seq": ""}

    # --- a genuine 5' construct is clipped in every read, so penetrance ~1 at
    #     the read start. A ramp starting mid-range is a heterogeneous footprint
    #     5' end, not a UMI.
    has_struct = pen[0] >= thr.struct_pen

    b = f
    rt5_len, ts5_seq, rt_nt, pen_rt = 0, "", False, 0.0

    # --- peel the RT base. Two independent signatures, either one suffices:
    #     (1) PARTIAL PENETRANCE. A UMI is a fixed-length element, present in
    #     every read, so its positions sit at full penetrance; an RT base is
    #     added in only a fraction of molecules, so it sits at *intermediate*
    #     penetrance. (Penetrance is exactly the survival function of the
    #     per-read 5' construct length: a fixed 2 nt UMI + RT gives reads with
    #     either 2 or 3 non-genomic bases -> a full-penetrance run of 2 then a
    #     partial position, which pins the UMI at 2.) This catches an RT base
    #     even when it is compositionally balanced.
    #     (2) A/T BIAS. The RT-added base is mostly A (reads as T after
    #     reverse-complementing), where a UMI base is uniform -- so a biased
    #     position is an RT base even at penetrance ~1, which (1) would miss.
    while b > 0 and (f - b) < thr.rt_max_len:
        j = b - 1
        if pen[j] < thr.footprint_pen or math.isnan(at[j]):
            break
        base = [at[k] for k in range(max(0, j - 3), j) if not math.isnan(at[k])]
        baseline = np.median(base) if base else CHANCE * 2
        partial = pen[j] < thr.rt_pen_hi                                  # (1) variable-length => RT
        biased = at[j] - baseline >= thr.at_jump or at[j] >= thr.at_abs   # (2) A/T bias => RT
        if partial or biased:
            rt5_len += 1                        # enzymatic, footprint-adjacent -> KEEP
            rt_nt = True
            pen_rt = max(pen_rt, float(pen[j]))
            b = j
        else:
            break

    def result(umi5, bc5, ts_len, layout=()):
        return {
            "umi5_len": umi5,               # random-templated, TOTAL
            "barcode5_seq": bc5,            # fixed-templated
            # the construct blocks in read order (5' -> footprint). umi5_len and
            # barcode5_seq collapse this to totals and cannot express a barcode
            # sitting BETWEEN two UMI stretches, which is a real design (iCLIP2:
            # [umi5 5nt][barcode ATTGGC][umi5 4nt][footprint]); the ordered layout
            # is what the display and the dedup mask must be built from.
            "p5_layout": list(layout),
            "rt5_len": min(rt5_len, thr.rt_max_len),   # enzymatic, KEEP
            "rt5_penetrance": round(float(pen_rt), 3),
            "rt_nt": rt_nt,
            "ts5_len": ts_len,             # enzymatic (template switch), TRIM
            "ts5_seq": ts5_seq,
        }

    # --- terminal RT base with no UMI in front. Then f lands on position 0 (its
    #     penetrance is below footprint_pen) and the RT-peel loop above never
    #     fires, so a fractional A/T-biased base at the very 5' end is otherwise
    #     absorbed into the footprint. Peel it when it is clearly A/T-biased (RT
    #     signature) and sits at intermediate penetrance while the next base is
    #     cleanly genomic -- i.e. a fractional extra base, not a footprint edge.
    if (f == 0 and rt5_len == 0 and n > 1 and not math.isnan(m5[1])
            and RT_TERMINAL_PEN_LO <= pen[0] < thr.footprint_pen
            and not math.isnan(at[0]) and at[0] >= thr.at_abs
            and pen[1] < RT_TERMINAL_PEN_LO):
        rt5_len, rt_nt, pen_rt = 1, True, float(pen[0])
        return result(0, "none", 0)

    if not has_struct:
        return result(0, "none", 0)

    # --- peel a G/C-rich template-switch / linker run adjacent to the footprint
    #     (D-Plex, iCLIP2). Full penetrance, G/C-biased -- where a UMI is balanced.
    #     The G-run has moderate-G/C positions interspersed among strong ones, so
    #     bridge a moderate position when its footprint-side neighbour is strong.
    while b > 0 and pen[b - 1] >= thr.struct_pen and not math.isnan(gc[b - 1]) \
            and len(ts5_seq) < thr.linker_max:
        j = b - 1
        strong = gc[j] >= thr.gc_template
        bridge = gc[j] >= thr.gc_mid and len(ts5_seq) > 0 and gc[j + 1] >= thr.gc_template
        if strong or bridge:
            ts5_seq = con[j] + ts5_seq
            b = j
        else:
            break

    # --- peel a G/C-rich non-templated run at the read start (SMARTer TSO G-run,
    #     RT-added C opposite the cap)
    s0 = 0
    while s0 < b and pen[s0] >= thr.struct_pen and not math.isnan(gc[s0]) \
            and gc[s0] >= thr.gc_template:
        s0 += 1
    ts_len = len(ts5_seq) + s0

    # --- remaining [s0, b): UMI (balanced, uniform/degenerate) + a barcode (const)
    kinds = [classify(ent[p], thr) for p in range(s0, b)]
    seqs = [con[p] for p in range(s0, b)]
    umi5_len, barcode5 = 0, "none"
    layout: list[dict] = []
    if b - s0:
        blocks = segment(kinds, seqs)
        rand = sum(len(s) for k, s in blocks if k == "random")
        degen = sum(len(s) for k, s in blocks if k == "degenerate")
        if rand + degen == 1 and not ts5_seq:
            rt5_len += 1                        # a lone non-genomic base is RT, not a UMI
            rt_nt = True
        elif rand:
            umi5_len = rand + degen             # degenerate positions inside a UMI count
        elif degen:
            ts_len += degen                     # degenerate-only block: non-templated, trim
        const = [s for k, s in blocks if k == "const"]
        if const:
            barcode5 = "".join(const)
        # keep the blocks in read order, so a barcode between two UMI stretches
        # survives into the call instead of being flattened into totals
        if umi5_len or barcode5 != "none":
            off = s0
            for k, s in blocks:
                role = ("barcode5" if k == "const" else
                        "umi5" if k in ("random", "degenerate") else None)
                if role:
                    # random and degenerate are distinct *kinds* but the same UMI:
                    # coalesce them, or a UMI with one skewed base becomes N{3}N{1}
                    if layout and layout[-1]["role"] == role == "umi5" \
                            and layout[-1]["offset"] + layout[-1]["len"] == off:
                        layout[-1]["len"] += len(s)
                    else:
                        blk = {"role": role, "offset": off, "len": len(s)}
                        if role == "barcode5":
                            blk["seq"] = s
                        layout.append(blk)
                off += len(s)
        if len(blocks) > 2:
            flags.append(f"5p_complex_layout:{'+'.join(f'{k}{len(s)}' for k, s in blocks)}")

    return result(umi5_len, barcode5, ts_len, layout)


# --- the 3' side -----------------------------------------------------------
def walk_left(match, comp, a0: int, thr_genomic: float, thr: Thresholds):
    """Positions before the anchor that are not genomic, nearest-first.

    Stops at the footprint, and also at the start of a homopolymer run: a
    poly(A) tail is constant in composition but variable in length, so it is
    not part of the fixed construct and must not be read as a barcode.
    """
    out = []
    run_base, run_len = None, 0
    for d in range(1, MAX_WALK + 1):
        j = a0 - d
        if j < 0 or math.isnan(match[j]):
            break
        # A footprint base is genomic AND variable. A constant linker base can
        # match the genome well above chance simply because one fixed base is
        # compared against a compositionally biased genomic neighbourhood, so a
        # high match rate alone must not end the walk.
        if match[j] >= thr_genomic and entropy(comp[j]) > thr.ent_const:
            break
        b, frac = dominant(comp[j])
        if frac >= HOMOPOLYMER_FRAC and b == run_base:
            run_len += 1
        else:
            run_base, run_len = (b, 1) if frac >= HOMOPOLYMER_FRAC else (None, 0)
        if run_len >= HOMOPOLYMER_MIN:
            # drop the run's positions; they are tail, not construct
            del out[len(out) - (HOMOPOLYMER_MIN - 1):]
            return out, f"poly{run_base}"
        out.append((d, match[j], entropy(comp[j]), consensus(comp[j])))
    return out, "none"


def const_run_right(comp, a0: int, limit: int, thr: Thresholds) -> int:
    """Length of the constant stretch starting at the anchor."""
    n = 0
    while a0 + n < limit:
        e = entropy(comp[a0 + n])
        if math.isnan(e) or e > thr.ent_const:
            break
        n += 1
    return n


def find_panel_in(block: str):
    """Leftmost panel adapter whose core occurs in `block`.

    A constant stretch sitting between the insert and a known adapter is a
    sample barcode when it is short (AGCTA, TAGAC), but a custom linker when it
    is long: `CACTCGGGCACCAAGGAC` followed by TruSeq read-through is one 18 nt
    linker, not an 18 nt barcode. Reporting the whole block as the adapter is
    what a trimmer needs.
    """
    best = None
    for name, core in PANEL_CORES.items():
        i = block.find(core)
        if i >= 0 and (best is None or i < best[1]):
            best = (name, i)
    if best and best[1] > MAX_BARCODE:
        return None
    return best


def call_3prime(profile: dict, plateau: float, flags: list[str],
                thr: Thresholds) -> tuple[dict, str]:
    thr_genomic = CHANCE + thr.genomic_frac * (plateau - CHANCE)
    # a real UMI/barcode base is clearly non-genomic (match ~ chance); a position
    # with an *intermediate* match is the footprint's own 3' edge, depressed by
    # alignment-edge leakage, not a construct base.
    thr_umi = CHANCE + thr.umi_match_frac * (plateau - CHANCE)
    res = {
        "adapter3_name": "unknown",
        "adapter3_seq": "unknown",
        "umi3_len": "unknown",
        "nt3_len": 0,
        "barcode3_seq": "unknown",
        "layout3_order": "unknown",
        "polyA_tail": "none",
    }
    kind = profile["anchor_kind"]
    if kind.startswith("panel:") and kind.split(":", 1)[1] in ("polyA", "polyG"):
        htype = kind.split(":", 1)[1]
        # A poly(A/G) tail is a genuine, reportable 3' feature and is trimmed like
        # one. The real adapter, if any, lies *beyond* it and is not visible in the
        # reads, so it cannot be named -- but that is no reason to refuse the whole
        # call. Report [footprint][poly-tail][adapter?] instead.
        res.update(adapter3_name="none_visible", adapter3_seq="none_visible",
                   umi3_len=0, barcode3_seq="none",
                   layout3_order=f"{htype}_tail,adapter?", polyA_tail=htype)
        flags.append(f"polyA_tail_adapter_not_visible:{htype}")
        return res, "raw"

    # An anchor backed by enough reads is usable even when those reads are a small
    # minority: a long insert hides the adapter from most reads, not from all of them,
    # and the construct is a property of the molecule, not of how far we sequenced.
    weak = (profile["frac_anchored"] < thr.min_anchor_frac
            and profile["n_anchored"] < MIN_ANCHORED_READS)
    if kind == "none" or weak:
        # no usable adapter: fall back to the read's own 3' end
        t3m = np.asarray(profile["t3_match"], dtype=float)
        t3c = profile["t3_comp"]
        n = 0
        while n < len(t3m) and not math.isnan(t3m[n]) and t3m[n] < thr_umi:
            n += 1
        if kind != "none":
            flags.append(f"weak_anchor:{kind}:{profile['frac_anchored']:.2f}")
        if n == 0:
            res.update(adapter3_name="none", adapter3_seq="none", umi3_len=0,
                       barcode3_seq="none", layout3_order="none")
            return res, "trimmed"

        # THE BOUNDARY IS A STEP, NOT A LEVEL -- the same rule the 5' side is built on.
        # A fixed-length 3' construct sits at the chance match rate for ALL of its
        # positions and then jumps to the genomic plateau: the walked region is FLAT.
        # What is not flat is a smear: when the deposit still carries an adapter that
        # was too rare to anchor on, the read's 3' end mixes reads whose construct +
        # visible-adapter length differs, and the match rate RAMPS up gradually. Walking
        # that ramp on a level test alone counts every position until the ramp happens
        # to cross the threshold, and reports the lot as one long UMI -- fabricating,
        # for one McGlincy-Ingolia library, an "11 nt UMI" that was really 5 nt of UMI
        # plus a 5 nt barcode plus an adapter base. Refuse the ramp instead.
        span = float(np.nanmax(t3m[:n]) - np.nanmin(t3m[:n])) if n > 1 else 0.0
        if span > TAIL_FLAT_MAX:
            flags.append(f"3p_tail_ramps_not_steps:{span:.2f}")
            return res, "unknown"

        kinds = [classify(entropy(t3c[k]), thr) for k in range(n)]
        if all(k in ("random", "degenerate") for k in kinds):
            if "degenerate" in kinds:
                flags.append("3p_degenerate_randomer")
            res.update(adapter3_name="none", adapter3_seq="none", umi3_len=n,
                       barcode3_seq="none", layout3_order="umi")
            return res, "adapter_trimmed_umi_retained"
        flags.append(f"3p_tail_no_adapter_kinds:{','.join(kinds)}")
        return res, "unknown"

    if profile["n_anchored"] < 200:
        flags.append("too_few_anchored_reads")
        return res, "unknown"

    a0 = profile["anchor_offset"]
    comp = profile["adap_comp"]
    match = np.asarray(profile["adap_match"], dtype=float)

    left, tail = walk_left(match, comp, a0, thr_genomic, thr)   # nearest-first, d = 1, 2, ...

    # If nothing between the footprint and the adapter is *clearly* non-genomic
    # (match ~ chance), the walk only caught the footprint's own 3' edge --
    # positions whose match is merely depressed by alignment-edge leakage -- not a
    # real UMI/barcode. A genuine construct always has at least one clearly
    # non-genomic base; a footprint edge does not. (A real UMI whose footprint-side
    # bases leak is kept, because its adapter-side bases are still clearly
    # non-genomic.)
    if left and all(m >= thr_umi for _, m, _, _ in left):
        left = []

    n_right = const_run_right(comp, a0, len(comp), thr)

    # constant stretch immediately left of the anchor
    n_left_const = 0
    for d, _m, e, _b in left:
        if not math.isnan(e) and e <= thr.ent_const:
            n_left_const = d
        else:
            break

    block = "".join(consensus(comp[a0 - n_left_const + i]) for i in range(n_left_const + n_right))
    hit = find_panel_in(block)

    if hit:
        name, idx = hit
        adapter_start = a0 - n_left_const + idx
        res["adapter3_name"] = name
        res["adapter3_seq"] = PANEL[name]
    else:
        adapter_start = a0 - n_left_const
        res["adapter3_name"] = "other"
        res["adapter3_seq"] = block
        flags.append("denovo_adapter_no_panel_match")

    # everything between the footprint and the adapter
    n_construct = len(left) - (a0 - adapter_start)
    if n_construct < 0:
        flags.append("adapter_start_left_of_footprint")
        n_construct = 0

    segs = []
    for i in range(n_construct):
        j = adapter_start - 1 - i          # walking left from the adapter
        segs.append((classify(entropy(comp[j]), thr), consensus(comp[j])))

    # segs[0] is adjacent to the adapter; reverse to insert-first order
    segs = segs[::-1]
    kinds = [k for k, _ in segs]
    seqs = [b for _, b in segs]

    # The gap between the footprint end and the adapter has a sharp mode when
    # the construct is a fixed-length UMI/barcode, and none when a variable
    # poly(A) tail sits in between.
    gap = profile.get("fpend_to_adapter_gap_hist", {})
    if gap and segs:
        tot = sum(gap.values())
        if max(gap.values()) / tot < 0.30 and tail == "none":
            flags.append("3p_construct_length_not_fixed")

    blocks = segment(kinds, seqs)
    rand3 = sum(len(s) for k, s in blocks if k == "random")
    degen3 = sum(len(s) for k, s in blocks if k == "degenerate")

    # A degenerate-only block is not a UMI (too few states to deduplicate on) but it
    # is still non-templated: it sits between the footprint end and the adapter, so
    # it must be TRIMMED or those bases stay on the footprint. The 5' side does the
    # same (ts_len += degen) -- e.g. SRR29632406, whose linker opens with a 2-nt
    # 'WW' ligation-bias reducer (A/T-degenerate: 61%A/38%T then 79%T).
    umi3 = rand3 + degen3 if rand3 else 0
    nt3 = degen3 if (degen3 and not rand3) else 0
    if nt3:
        flags.append(f"3p_degenerate_nontemplated:{nt3}")

    const = "".join(s for k, s in blocks if k == "const")
    order = [("barcode" if k == "const" else "umi") for k, _ in blocks]
    if tail != "none":
        order = ["polyA_tail"] + order
    order.append("adapter")

    res.update(
        umi3_len=umi3,
        nt3_len=nt3,                    # degenerate non-templated: trim, do not dedup
        barcode3_seq=const or "none",
        layout3_order=",".join(order),
        polyA_tail=tail,
        # where the adapter actually starts, relative to the profile's anchor
        # (a0). The de-novo seed can begin a few nt inside the adapter, so this is
        # negative when the walk extended the adapter left through fixed bases.
        adapter3_start_rel=adapter_start - a0,
    )
    if len(blocks) > 2:
        flags.append(f"3p_complex_layout:{'+'.join(f'{k}{len(s)}' for k, s in blocks)}")

    return res, "raw"


# --- the functional view ---------------------------------------------------
def functional_view(five: dict, three: dict) -> dict:
    """Recast the call into the four functional categories and derive the two
    things that actually drive processing: the trim boundaries and the dedup
    mask. Categories: footprint | random-templated (UMI/spacer) | fixed-templated
    (barcode/linker/adapter) | enzymatic (RT nt, kept; template-switch, trimmed).

    The '+/-1 UMI vs spacer' ambiguity lives *inside* random-templated, so it
    changes none of the derived quantities below.
    """
    # --- 5' construct, read order (position 0 -> footprint)
    segs: list[dict] = []
    if five["ts5_len"]:
        segs.append({"cat": "enzymatic", "role": "template_switch",
                     "len": five["ts5_len"], "fate": "trim"})
    # walk the construct in read order, so [umi][barcode][umi] renders as it is
    for blk in five.get("p5_layout") or []:
        if blk["role"] == "barcode5":
            segs.append({"cat": "fixed-templated", "role": "barcode5",
                         "len": blk["len"], "seq": blk["seq"], "fate": "trim"})
        else:
            segs.append({"cat": "random-templated", "role": "umi5",
                         "len": blk["len"], "fate": "trim+dedup"})
    if not (five.get("p5_layout") or []):     # no ordered layout: fall back to totals
        if five["barcode5_seq"] not in ("none", "", "unknown"):
            segs.append({"cat": "fixed-templated", "role": "barcode5",
                         "len": len(five["barcode5_seq"]), "seq": five["barcode5_seq"],
                         "fate": "trim"})
        if five["umi5_len"]:
            segs.append({"cat": "random-templated", "role": "umi5",
                         "len": five["umi5_len"], "fate": "trim+dedup"})
    if five["rt5_len"]:
        segs.append({"cat": "enzymatic", "role": "rt_untemplated",
                     "len": five["rt5_len"], "penetrance": five["rt5_penetrance"],
                     "fate": "keep"})

    umi5 = five["umi5_len"]
    # trim to reach the footprint-you-keep: everything 5' of the footprint except
    # the retained RT base
    trim_5p = five["ts5_len"] + \
        (len(five["barcode5_seq"]) if five["barcode5_seq"] not in ("none", "", "unknown") else 0) + \
        umi5

    # --- 3' construct, footprint -> read end
    segs3: list[dict] = []
    umi3 = _as_int(three.get("umi3_len"), 0)
    nt3 = _as_int(three.get("nt3_len"), 0)
    bc3 = three.get("barcode3_seq", "none")
    if three.get("polyA_tail", "none") not in ("none", ""):
        segs3.append({"cat": "enzymatic", "role": "polyA_tail",
                      "len": "variable", "fate": "trim"})
    if nt3:
        segs3.append({"cat": "enzymatic", "role": "nontemplated3",
                      "len": nt3, "fate": "trim"})
    if umi3:
        segs3.append({"cat": "random-templated", "role": "umi3",
                      "len": umi3, "fate": "trim+dedup"})
    if bc3 not in ("none", "", "unknown"):
        segs3.append({"cat": "fixed-templated", "role": "barcode3",
                      "len": len(bc3), "seq": bc3, "fate": "trim"})
    # only report an adapter segment when one is actually present and identified;
    # an already-trimmed deposit still carries a valid 5' construct (and sometimes
    # a retained 3' UMI), and a poly(A) library's adapter lies beyond the tail and
    # is not visible -- neither should invent an adapter segment.
    adap_name = three.get("adapter3_name", "unknown")
    if adap_name not in ("none", "", "unknown", "none_visible"):
        segs3.append({"cat": "fixed-templated", "role": "adapter",
                      "seq": three.get("adapter3_seq", "unknown"),
                      "name": adap_name, "fate": "trim"})

    dedup_len = umi5 + umi3
    mask = []
    if umi5:
        # spell out the 5' construct when a barcode splits the UMI: N{5}ATTGGC N{4}
        # is the extractable pattern; a bare N{9} would take the barcode's bases
        lay = five.get("p5_layout") or []
        if sum(1 for b in lay if b["role"] == "umi5") > 1:
            mask.append("5':" + "".join(
                b["seq"] if b["role"] == "barcode5" else f"N{{{b['len']}}}" for b in lay))
        else:
            mask.append(f"5':N{{{umi5}}}")
    if umi3:
        mask.append(f"3':N{{{umi3}}}")

    return {
        "segments_5p": segs,
        "footprint": {"category": "footprint",
                      "keeps_enzymatic_rt_nt": five["rt5_len"]},
        "segments_3p": segs3,
        "trim_5p": trim_5p,                          # nt to clip from the 5' end
        "footprint_retains_rt_nt": five["rt5_len"],  # kept, not trimmed
        "trim_3p_adapter": three.get("adapter3_seq", "unknown"),
        # everything between the footprint end and the adapter -- all of it comes off
        "trim_3p_construct": umi3 + nt3 + (len(bc3) if bc3 not in ("none", "", "unknown") else 0),
        "dedup_umi_len": dedup_len,                  # total random-templated content
        "dedup_umi_mask": " ".join(mask) or "none",
    }


def _as_int(v, d: int = 0) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        return d


# --- the call --------------------------------------------------------------
# Anything that leaves a segment genuinely ambiguous is reported as undetermined
# rather than guessed. A de-novo adapter or a variable-length tail is
# informative, not ambiguous, so neither disqualifies a call.
HARD_FLAGS = (
    "3p_tail_no_adapter",
    "5p_never_reaches_genomic",
    "homopolymer_anchor",
    "too_few_anchored_reads",
    "adapter_start_left_of_footprint",
)


def infer(profile: dict, thr: Thresholds | None = None) -> dict:
    """Call the read architecture from one profile dict (`qc.profile.profile_bam`).

    Returns the call: `status` is 'ok' or 'undetermined' (with a `reason`); an
    'ok' call also carries the `functional` block whose `trim_5p` /
    `dedup_umi_len` drive `arch.trim`. Raises ValueError if the profile is not
    one -- a malformed input is a bug, not a refusal.
    """
    thr = thr or Thresholds()
    flags: list[str] = []
    for key in ("label", "n_used"):
        if key not in profile:
            raise ValueError(f"not a read-architecture profile: missing {key!r}")
    out: dict[str, Any] = {
        "run_accession": profile["label"],
        "n_reads_used": profile["n_used"],
        # Record what this call was decided against. The architecture figure draws
        # these as its reference lines, so it audits the decision that was actually
        # made rather than the defaults it would otherwise assume -- and a call kept
        # on disk stays interpretable after someone re-tunes the config.
        "thresholds": asdict(thr),
    }

    if profile["n_used"] < thr.min_reads:
        out.update(status="undetermined",
                   reason=f"only {profile['n_used']} usable alignments")
        return out

    missing = [k for k in _REQUIRED if k not in profile]
    if missing:
        raise ValueError(f"profile is missing field(s): {', '.join(missing)}")

    plateau = genomic_plateau(profile)
    if plateau is None:
        out.update(status="undetermined", reason="reads too short to profile")
        return out
    if plateau < thr.min_plateau:
        out.update(
            status="undetermined",
            reason=f"no genomic plateau (deep match rate {plateau:.2f}); alignments too noisy",
        )
        return out

    # Genuine ribosome footprints vary in length (the ribosome protects ~28-32 nt
    # with natural spread). A single sharp length spike means the "footprints" are
    # a fixed-length artefact -- adapter dimers, a fixed contaminant -- not real
    # footprints, and any 5' construct read off them is fabricated.
    fph = profile.get("footprint_len_hist") or {}
    if fph:
        tot = sum(fph.values())
        mode = max(fph, key=lambda k: fph[k])
        mode_frac = fph[mode] / tot
        if mode_frac > thr.footprint_uniform_max:
            out.update(
                status="undetermined",
                reason=f"footprint length is a single {mode} nt spike ({mode_frac:.0%} of "
                       f"reads); not genuine ribosome footprints (adapter-dimer / artefact)",
            )
            return out

    # Genuine footprints spread across thousands of loci, so no single 5'
    # coordinate holds more than ~1-2% of reads (the busiest genuine library here
    # tops out at 38%, and that one has 37k distinct loci -- a single hyper-
    # expressed transcript). A *majority* of reads at one exact 5' coordinate is a
    # single over-represented species -- an adapter dimer aligning to an adapter-
    # like locus, a spike-in, a contaminant -- whose architecture is fabricated.
    conc = profile.get("top5p_locus_frac")
    if conc is not None and conc > thr.single_locus_max:
        out.update(
            status="undetermined",
            reason=f"{conc:.0%} of reads share one 5' locus "
                   f"({profile.get('n_distinct_5p_loci', 0):,} distinct); a single over-"
                   f"represented species (adapter dimer / contaminant), not footprints",
        )
        return out

    five = call_5prime(profile, plateau, flags, thr)
    three, deposit = call_3prime(profile, plateau, flags, thr)

    out.update(
        genomic_plateau=round(plateau, 3),
        deposit_state=deposit,
        umi5_len=five["umi5_len"],
        p5_layout=five["p5_layout"],
        nontemplated5_len=five["rt5_len"] + five["ts5_len"],
        template_switch5_seq=five["ts5_seq"] or "none",
        barcode5_seq=five["barcode5_seq"],
        rt_untemplated_5p=five["rt_nt"],
        rt_penetrance=five["rt5_penetrance"],
        footprint_len_mode=mode_of(profile["footprint_len_hist"]),
        read_len_mode=mode_of(profile["read_len_hist"]),
        anchor=profile["anchor_kind"],
        frac_anchored=round(profile["frac_anchored"], 3),
        **three,
    )
    out["flags"] = flags
    # emit the functional view for every readable deposit, not just raw ones: a
    # deposit whose 3' adapter was already trimmed still carries a valid 5'
    # construct (UMI / RT nt) and sometimes a retained 3' UMI, and those should
    # be shown rather than hidden behind an "already trimmed" placeholder.
    if deposit in ("raw", "trimmed", "adapter_trimmed_umi_retained"):
        out["functional"] = functional_view(five, three)

    hard = [f for f in flags if f.startswith(HARD_FLAGS)]
    if hard or deposit == "unknown":
        out["status"] = "undetermined"
        out["reason"] = "; ".join(hard) or "3' side could not be segmented"
    else:
        out["status"] = "ok"
    LOG.debug("%s: %s%s", out["run_accession"], out["status"],
              f" ({', '.join(flags)})" if flags else "")
    return out


# --- the pretty-printer, shared by the TSV report and the plot title --------
def _arch_token_5p(s: dict) -> str:
    r = s["role"]
    if r == "umi5":
        return f"[UMI,{s['len']}nt]"
    if r == "rt_untemplated":
        pen = s.get("penetrance")
        return f"[RT,{pen * 100:.0f}%]" if pen else f"[RT,{s.get('len', 1)}nt]"
    if r == "barcode5":
        return f"[barcode,{s.get('seq', '?')}]"
    if r == "template_switch":
        return f"[TS,{s['len']}nt]"
    return f"[{r},{s.get('len', '')}]"


def _arch_token_3p(s: dict) -> str:
    r = s["role"]
    if r == "umi3":
        return f"[UMI,{s['len']}nt]"
    if r == "barcode3":
        return f"[barcode,{s.get('seq', '?')}]"
    if r == "polyA_tail":
        return "[polyA]"
    if r == "adapter":
        return f"[{adapter_disp(s.get('name'))}]"
    return f"[{r}]"


def architecture_string(call: dict) -> str:
    """A compact 5'->3' rendering of the inferred architecture, e.g.
       5'-[UMI,2nt]-[RT,40%]-[footprint,~28nt]-[UMI,5nt]-[barcode,GATCA]-[TruSeq]-3'
    """
    if call.get("status") != "ok":
        return f"{call.get('status', '?').upper()}: {call.get('reason', '')[:120]}"
    fp = f"[footprint,~{call.get('footprint_len_mode')}nt]"
    fn = call.get("functional") or {}
    if not fn:
        return f"5'-{fp}-[adapter trimmed]-3'   (deposit already trimmed)"
    parts = [_arch_token_5p(s) for s in fn.get("segments_5p", [])]
    parts.append(fp)
    parts += [_arch_token_3p(s) for s in fn.get("segments_3p", [])]
    arch = "5'-" + "-".join(p for p in parts if p) + "-3'"
    if call.get("deposit_state") != "raw":   # adapter was pre-trimmed off the 3'
        arch += "   (3' adapter pre-trimmed)"
    return arch

"""Trim a FASTQ using the architecture inferred by `ribomine.arch.infer`.

Given the raw reads and the sample's architecture call, this removes exactly
what the inferred architecture says is not footprint, and moves the
random-templated (UMI) content into the read name for downstream deduplication:

    5'  [template-switch][barcode5][UMI5]   <-- trimmed (UMI5 -> header)
        [RT nt]                             <-- kept (enzymatic, part of the molecule)
        [============ footprint ==========] <-- kept: the output
        [UMI3][barcode3]                    <-- trimmed (UMI3 -> header)
        [adapter | poly(A) tail]            <-- trimmed

The read name gets `_<UMI5><UMI3>` (umi_tools style). Reads left shorter than
`min_len` after trimming are dropped and counted.

The call is inferred once, by the architecture stage, and handed in here -- a
run of 10-100M reads is streamed record by record, never held in memory.
"""
from __future__ import annotations

import logging
import os
from functools import lru_cache

from cutadapt.align import Aligner, EndSkip

from ..utils import open_fastq

LOG = logging.getLogger("ribomine.arch.trim")

# values that mean "this element is absent from the call"
_ABSENT = ("none", "unknown", "", None)

# a library that cannot show its 3' scaffold to this share of its reads will lose them
NO_ADAPTER_WARN = 0.30

# cutadapt's error model: a RATE, not a count -- so a short partial adapter at the
# read end must match exactly (0.12 * 7 = 0 errors) while a full-length one tolerates
# sequencing error. Below ~9 nt of overlap the budget is zero.
ERROR_RATE = 0.12

# A 3' adapter may begin anywhere in the read (skip the read's prefix: the footprint),
# may be cut short by the read's end (skip the adapter's suffix), and -- the one that is
# easy to forget -- may be FOLLOWED by more sequence (skip the read's suffix).
#
# That last skip is not optional. Without QUERY_STOP the aligner requires the alignment
# to reach the end of the READ, i.e. the adapter must be the last thing in it. Any
# library that sequences past the adapter into an index, a barcode or a second adapter
# then matches in ZERO reads -- and `discard_untrimmed` throws the entire library away.
# Measured on SRR25706716 ([footprint][Ingolia linker][12 nt][TruSeq]): the linker is an
# exact substring of 98% of reads, and the matcher found it in 0% of them; 53.7 M reads
# became a 23 KB BAM. The three skips together ARE cutadapt's `-a` (its BACK flag), which
# is what this code always meant to be.
_ADAPTER_FLAGS = EndSkip.QUERY_START | EndSkip.QUERY_STOP | EndSkip.REFERENCE_END

# Score indels out of reach (the error budget never exceeds ~5 over a scaffold this
# long), so the match stays Hamming-only -- the same model the hand-written matcher
# this replaced used, and so the same cut points. Drop this to 1 to let the aligner
# absorb an indel in the adapter, which is a real change to where reads get cut.
_NO_INDELS = 1000


@lru_cache(maxsize=64)
def _aligner(adapter: str, min_overlap: int) -> Aligner:
    """cutadapt's C aligner for one adapter, built once and reused for every read.

    Constructing it is what costs; matching against it is ~8x faster than the
    equivalent Python loop, and a run streams tens of millions of reads past it.
    A handful of adapters exist per run (the scaffold and the bare adapter, per
    sample), so the cache never grows.
    """
    return Aligner(adapter, max_error_rate=ERROR_RATE, flags=_ADAPTER_FLAGS,
                   min_overlap=min_overlap, indel_cost=_NO_INDELS)


def find_adapter(read: str, adapter: str, min_start: int, min_overlap: int = 7) -> int:
    """Where `adapter` starts in `read` at or after `min_start`, or -1.

    Handles the adapter running off the read's end (a partial match at the end is
    still a match, down to `min_overlap` nt). Where several placements are legal
    the aligner takes the one matching the most bases -- so a full-length adapter
    late in the read beats a chance 7-mer early in it.
    """
    lo = max(0, min_start)
    m = _aligner(adapter, max(1, min_overlap)).locate(read[lo:])
    return lo + m[2] if m else -1     # m[2] = where the adapter starts in read[lo:]


def trim_polyA(read: str, min_run: int = 6, max_gap: int = 1) -> int:
    """Index where a poly(A) (or poly(T)) tail ENDING THE READ starts, or len(read).

    Terminal-only, and used where the caller has already cut everything past the tail
    off (so the tail is the last thing in the string by construction). When the read may
    still run on into a construct, use `find_polyA`.

    A conserved base embedded in the tail is stepped over: some poly(A)-primed libraries
    carry an [A..]-C-[A..] linker before the adapter (SRR18113808: a fixed C between two
    ~10 nt A-runs), and a terminal-only cut removes only the adapter-proximal run, leaving
    the footprint-proximal run AND the C on the read -- ~11 nt that STAR then soft-clips,
    so `mean_mapped_len` falls that far below `mean_footprint_len`. After the terminal run
    the search hops up to `max_gap` non-run bases, but ONLY when another run of >= min_run
    of the same base resumes just past the gap. A genomic footprint does not carry a
    >= min_run A-run abutting the tail, so that boundary is never crossed (the internal-run
    guard the footprint-protection tests rest on is preserved).
    """
    n = len(read)
    for base in ("A", "T"):
        i = n
        while i > 0 and read[i - 1] == base:
            i -= 1
        if n - i < min_run:
            continue                         # no terminal run of this base
        while True:                          # hop a short gap only if a real run resumes
            k, g = i, 0
            while k > 0 and read[k - 1] != base and g < max_gap:
                k -= 1
                g += 1
            if k == 0 or read[k - 1] != base:
                break                        # nothing but the gap: stop at the run we have
            j = k
            while j > 0 and read[j - 1] == base:
                j -= 1
            if k - j < min_run:
                break                        # the stretch past the gap is too short to be tail
            i = j
        return i
    return n


def find_polyA(read: str, *, min_start: int = 0, min_run: int = 6) -> int:
    """Index where the 3'-most poly(A)/poly(T) tail starts, or len(read) if there is none.

    Unlike `trim_polyA`, the tail need NOT end the read. It sits between the footprint
    and the construct, so a read long enough to sequence *through* the tail carries it in
    the MIDDLE -- and a terminal-only search then finds nothing and leaves the footprint
    with the tail AND the construct still attached. That is not a yield question: those
    reads are no longer genomic, so they simply fail to align. Measured on SRR30214250
    (an adapter RiboMine's panel cannot name, so there is nothing else to anchor on): 20%
    of reads kept their whole 3' end, and they mapped at 3% against 14% for the reads that
    were cut.

    The tail is taken as the RIGHTMOST qualifying run, not the leftmost: a chance A-run
    inside the footprint lies to the LEFT of the real tail, so scanning from the 3' end
    walks past the footprint's own sequence rather than cutting into it.
    """
    n = len(read)
    lo = max(0, min_start)
    for base in ("A", "T"):
        i = n
        while i > lo:
            if read[i - 1] != base:
                i -= 1
                continue
            j = i                                  # walk to the run's start
            while j > lo and read[j - 1] == base:
                j -= 1
            if i - j >= min_run:
                return j                           # first run met scanning right = 3'-most
            i = j
    return n


def plan_from_call(call: dict, *, min_overlap: int = 7) -> dict:
    """Extract the trimming plan from an `infer()` result."""
    fn = call.get("functional") or {}
    bc3 = call.get("barcode3_seq", "none")
    return {
        "status": call.get("status"),
        "trim_5p": fn.get("trim_5p", 0),
        "umi5": call.get("umi5_len", 0) or 0,
        # ordered 5' blocks: the UMI may be split around a barcode, in which case the
        # UMI bases are NOT the last umi5 nt before the footprint
        "p5_layout": call.get("p5_layout") or [],
        "umi3": _as_int(call.get("umi3_len"), 0),
        "nt3": _as_int(call.get("nt3_len"), 0),
        "bc3": bc3,
        # the barcode's actual bases -- part of the fixed 3' scaffold we anchor on
        "bc3_seq": bc3 if bc3 not in _ABSENT else "",
        "adapter_name": call.get("adapter3_name", "unknown"),
        "adapter_seq": call.get("adapter3_seq", ""),
        "polyA": call.get("polyA_tail", "none"),
        "fp_mode": call.get("footprint_len_mode") or 28,
        # nt of the fixed scaffold that must be visible before we will cut on it.
        # 7 nt of a known string, anchored at the read end, matches by chance once in
        # 16,000 reads. Lowering it trades that specificity for depth on a library
        # whose molecules barely fit the read (see trim_read).
        "min_overlap": max(1, int(min_overlap)),
    }


def _as_int(v, d: int = 0) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        return d


def trim_read(seq: str, qual: str, p: dict) -> tuple[str, str, str, bool | None] | None:
    """Apply the plan to one read.

    -> (trimmed_seq, trimmed_qual, umi, adapter_found), or None when nothing of the
    footprint is left. `adapter_found` is None when the architecture expects no
    ligated adapter (an already-trimmed deposit, or a poly(A) library whose adapter
    lies beyond the tail); otherwise it says whether this particular read actually
    carried one.

    That flag matters because a read WITHOUT the expected adapter did not sequence
    through to the end of the molecule: its 3' end is the read's end, not the
    footprint's, and any 3' UMI sitting between them was never read. So such a read
    has neither a trustworthy footprint boundary nor a complete UMI -- see
    `trim_fastq`, which decides what to do with it.
    """
    trim5 = p["trim_5p"]
    lay = [b for b in p["p5_layout"] if b["role"] == "umi5"]
    if lay:
        # take each UMI block at its own offset -- a barcode may sit between them
        umi5 = "".join(seq[b["offset"]:b["offset"] + b["len"]] for b in lay)
    else:
        umi5 = seq[max(0, trim5 - p["umi5"]):trim5] if p["umi5"] else ""
    s5 = trim5
    bc3len = len(p["bc3"]) if p["bc3"] not in _ABSENT else 0

    name = p["adapter_name"]
    umi3 = ""
    adapter_found: bool | None = None
    if name not in _ABSENT and name != "none_visible" and p["adapter_seq"]:
        # Anchor on the whole FIXED 3' scaffold -- the barcode AND the adapter -- not
        # on the adapter alone. The sample barcode is just as constant and just as
        # known, and it sits `bc3len` nt CLOSER to the insert, so it is visible in
        # reads whose adapter has already run off the end. That matters enormously
        # when the molecule is about as long as the read: on SRR11945406 (46 nt reads,
        # a 44 nt molecule) only 2.5% of reads show >= 7 nt of TruSeq, but 35% show
        # >= 7 nt of [barcode + TruSeq]. Anchoring on the adapter alone therefore threw
        # away 97% of a perfectly good library. Specificity is unchanged: a 7-nt match
        # to a known string, anchored at the read end, happens by chance once in 16,000.
        lo = max(s5, p["fp_mode"] - 6)
        scaffold = p["bc3_seq"] + p["adapter_seq"]
        s = find_adapter(seq, scaffold, min_start=lo, min_overlap=p["min_overlap"])
        if s < 0 and bc3len:
            # The scaffold match must start AT the barcode, so a sequencing error in
            # the barcode blocks it -- and at a 7-nt overlap the error budget is zero.
            # Those reads were trimmable before (the adapter alone matched), so fall
            # back to it and derive the same cut point. This keeps the scaffold anchor
            # a strict superset: it can only ever find MORE reads, never fewer.
            # search from the same `lo` the adapter-only anchor always used, or reads
            # whose adapter sits right at that bound would be lost to the tighter start
            a = find_adapter(seq, p["adapter_seq"], min_start=lo,
                             min_overlap=p["min_overlap"])
            if a >= bc3len:
                s = a - bc3len
        adapter_found = s >= 0
        if s < 0:
            # The scaffold ran off the end of the read -- but a poly(A) TAIL is
            # SELF-ANCHORING. It sits between the footprint and the construct, so it
            # marks the footprint's 3' end whether or not the adapter is visible: the
            # read simply ends inside the tail. Using it rescues exactly the reads a
            # long footprint pushed the adapter off (SRR19641906: 38% -> ~100%
            # retention), and removes the length bias that discarding them created.
            # It only works when nothing FIXED sits between tail and adapter, though --
            # a 3' UMI or barcode beyond the tail was never sequenced in these reads,
            # so their boundary genuinely cannot be placed.
            e3 = len(seq)                      # scaffold not in this read: keep to end
            if p["polyA"] not in _ABSENT and not p["umi3"] and not bc3len:
                # Peel any adapter fragment too short to have anchored on (1..min_overlap-1
                # nt), because it sits BETWEEN the tail and the read end and would
                # otherwise hide it -- TruSeq begins with an A, so a read ending
                # "...AAAAAAAG" has no trailing A-run at all. Only then read the tail.
                end = len(seq)
                adap = p["adapter_seq"]
                for k in range(min(p["min_overlap"] - 1, len(adap)), 0, -1):
                    if seq.endswith(adap[:k]):
                        end = len(seq) - k
                        break
                e3 = trim_polyA(seq[:end])
                adapter_found = e3 < end       # boundary located iff a tail was there
        else:
            # the scaffold starts at the barcode; the UMI (and any degenerate
            # non-templated block) sit between the footprint and it, and all of it
            # comes off the footprint
            umi3 = seq[s - p["umi3"]:s] if p["umi3"] else ""
            e3 = s - p["umi3"] - p["nt3"]
            # ... and a poly(A) tail sits between the FOOTPRINT and that construct
            # (insert-first: [footprint][poly-A][nt3][UMI][barcode][adapter]). It is
            # enzymatically added, variable in length, and NOT genomic, so it has to
            # come off too. Cutting only the fixed construct leaves the whole tail on
            # the footprint -- on SRR19641906, 99.6% of trimmed reads still ended in a
            # run of >= 6 A's, which end-to-end alignment would then have to explain.
            if p["polyA"] not in _ABSENT:
                e3 = trim_polyA(seq[:e3])
    elif p["polyA"] not in ("none", "", None):
        # No adapter to anchor on -- either the library has none, or (SRR30214250) it has
        # one this panel cannot name. The poly(A) tail anchors itself: it lies between the
        # footprint and whatever follows, so it marks the 3' boundary WITHOUT us having to
        # know what follows. `find_polyA`, not `trim_polyA`: the construct beyond the tail
        # was sequenced in these reads, so the tail is not at the read's end.
        e3 = find_polyA(seq, min_start=s5)
        # A read with no tail at all did not sequence through to the end of its molecule:
        # its 3' end is the read's end, not the footprint's. That is the same thing a
        # missing adapter means, and it gets the same answer -- see trim_fastq. Keeping
        # them was the bug: they carry the tail and the construct, and simply do not align.
        adapter_found = e3 < len(seq)
    else:
        # adapter already trimmed off the deposit; a retained 3' UMI sits at the end
        if p["umi3"]:
            umi3 = seq[len(seq) - p["umi3"]:]
            e3 = len(seq) - p["umi3"]
        else:
            e3 = len(seq)

    if e3 <= s5:
        return None
    return seq[s5:e3], qual[s5:e3], umi5 + umi3, adapter_found


def trim_fastq(in_fq: str, call: dict, out_fq: str, *, min_len: int = 20,
               label: str = "", discard_untrimmed: bool = True,
               min_overlap: int = 7) -> dict:
    """Trim `in_fq` according to the architecture `call`, writing `out_fq`.

    Both paths may be plain or `.gz`. The FASTQ is streamed record by record.
    An `undetermined` call is not fatal: its best-effort plan (whatever the
    caller could read off the profile) is applied, with a warning -- the
    pipeline decides whether such a sample is carried forward at all
    (`architecture.process_undetermined`).

    Reads that do not carry the expected adapter
    --------------------------------------------
    When the architecture says the library has a ligated 3' adapter but a given
    read has none, the sequencer never reached the end of the molecule. Two things
    follow, and both are silent corruption if ignored:

    * the read's 3' end is an artefact of the read length, not the footprint's own
      end -- so it is not a complete ribosome footprint at all (a 28-32 nt
      footprint in a 50 nt read essentially always shows its adapter; the ones that
      do not are long inserts: rRNA, mRNA fragments, other junk);
    * any 3' UMI sits *beyond* the read end and was never sequenced, so the UMI we
      would write into the header is **shorter than for every other read**.

    That second point is what makes this a correctness bug rather than a
    yield question: `umi_tools dedup` requires a fixed UMI length and aborts with
    `not all umis are the same length` the moment one read is short. Discarding
    these reads (`discard_untrimmed`, the default, and standard practice in
    ribo-seq -- it is cutadapt's `--discard-untrimmed`) removes the problem at its
    source. If they are kept instead, the absent UMI bases are written as `N` so
    that every UMI in the file has the length the architecture declared.
    """
    status = call.get("status")
    if status != "ok":
        LOG.warning("%s: architecture is '%s' (%s); trimming with best-effort plan",
                    label or os.path.basename(in_fq), status,
                    str(call.get("reason", ""))[:70])
    p = plan_from_call(call, min_overlap=min_overlap)
    # every UMI written must have exactly this length, or umi_tools will refuse the file
    umi_len = p["umi5"] + p["umi3"]

    n_in = n_out = n_short = n_no_adapter = n_padded = 0
    len_in = len_out = 0
    with open_fastq(in_fq) as fh, open_fastq(out_fq, "wt") as out:
        while True:
            h = fh.readline()
            if not h or not h.strip():
                break
            seq_line, plus_line, qual_line = fh.readline(), fh.readline(), fh.readline()
            if not seq_line or not plus_line or not qual_line:
                raise ValueError(f"{in_fq}: truncated FASTQ record after {n_in} reads")
            seq = seq_line.rstrip("\n")
            qual = qual_line.rstrip("\n")
            n_in += 1
            len_in += len(seq)
            r = trim_read(seq.upper(), qual, p)
            if r is None or len(r[0]) < min_len:
                n_short += 1
                continue
            ts, tq, umi, adapter_found = r
            if adapter_found is False:
                n_no_adapter += 1
                if discard_untrimmed:
                    continue
                if len(umi) < umi_len:      # the 3' UMI was never sequenced
                    umi = umi + "N" * (umi_len - len(umi))
                    n_padded += 1
            name = h[1:].split()[0].rstrip("\n")
            tag = f"_{umi}" if umi else ""
            out.write(f"@{name}{tag}\n{ts}\n+\n{tq}\n")
            n_out += 1
            len_out += len(ts)

    stats = {
        "label": label or os.path.basename(in_fq),
        "architecture_status": status,
        "trim_5p": p["trim_5p"],
        "umi5": p["umi5"],
        "umi3": p["umi3"],
        "umi_len": umi_len,
        "adapter": p["adapter_name"],
        "polyA": p["polyA"],
        "n_reads_in": n_in,
        "n_reads_out": n_out,
        "n_dropped_short": n_short,
        "n_no_adapter": n_no_adapter,
        "n_dropped_untrimmed": n_no_adapter if discard_untrimmed else 0,
        "n_umi_padded": n_padded,
        "discard_untrimmed": discard_untrimmed,
        "min_overlap": p["min_overlap"],
        "frac_no_adapter": round(n_no_adapter / max(n_in, 1), 4),
        "frac_kept": round(n_out / max(n_in, 1), 4),
        "mean_len_in": round(len_in / max(n_in, 1), 1),
        "mean_len_out": round(len_out / max(n_out, 1), 1),
    }
    extra = ""
    if n_no_adapter:
        extra = (f"; {n_no_adapter:,} without the adapter "
                 f"{'dropped' if discard_untrimmed else f'kept, {n_padded:,} UMIs N-padded'}")
    LOG.info("%s: %s/%s reads kept (%s dropped < %d nt%s); mean len %s -> %s nt; "
             "trim5=%s umi5=%s umi3=%s adapter=%s polyA=%s",
             stats["label"], f"{n_out:,}", f"{n_in:,}", f"{n_short:,}", min_len, extra,
             stats["mean_len_in"], stats["mean_len_out"], p["trim_5p"], p["umi5"],
             p["umi3"], p["adapter_name"], p["polyA"])

    # Losing most of a library must never be quiet. When the molecule is about as long
    # as the read, the fixed 3' scaffold runs off the end and the read cannot be cut at
    # the footprint boundary -- so those reads are (rightly) discarded, and the dataset
    # silently arrives at a fraction of its advertised depth. Say so, with the lever.
    if stats["frac_no_adapter"] >= NO_ADAPTER_WARN and p["adapter_name"] not in _ABSENT:
        LOG.warning(
            "%s: the 3' scaffold is missing from %.0f%% of reads -- the molecule is about "
            "as long as the read, so it runs off the end. Only %s of %s reads survive "
            "trimming, and THEY ARE LENGTH-BIASED: it is precisely a long footprint that "
            "pushes the scaffold past the read end, so the survivors skew short. That is a "
            "property of the library, not a failure -- but do not read the footprint-length "
            "distribution of this dataset off the BAM. Lowering process.adapter_min_overlap "
            "(now %d) recovers depth at the cost of specificity.",
            stats["label"], 100 * stats["frac_no_adapter"], f"{n_out:,}", f"{n_in:,}",
            p["min_overlap"])
    return stats

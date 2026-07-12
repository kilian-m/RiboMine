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

from ..utils import open_fastq

LOG = logging.getLogger("ribomine.arch.trim")

# values that mean "this element is absent from the call"
_ABSENT = ("none", "unknown", "", None)


def max_mm(overlap: int) -> int:
    """Mismatches tolerated over an adapter overlap of `overlap` nt (cutadapt's
    error model: a rate, not a count, so a short partial adapter at the read end
    is matched exactly and a full-length one tolerates sequencing errors)."""
    return int(0.12 * overlap)


def find_adapter(read: str, adapter: str, min_start: int, min_overlap: int = 7) -> int:
    """Leftmost i >= min_start where read[i:] is a prefix of `adapter` (cutadapt
    error model). Handles the adapter running off the read end."""
    n, m = len(read), len(adapter)
    for i in range(max(0, min_start), n - min_overlap + 1):
        ov = min(m, n - i)
        if ov < min_overlap:
            break
        mm = 0
        ok = True
        for a, b in zip(read[i:i + ov], adapter[:ov]):
            if a != b and (mm := mm + 1) > max_mm(ov):
                ok = False
                break
        if ok:
            return i
    return -1


def trim_polyA(read: str, min_run: int = 6) -> int:
    """Return the index where a 3' poly(A) (or poly(T)) tail starts, or len(read)."""
    n = len(read)
    for base in ("A", "T"):
        i = n
        while i > 0 and read[i - 1] == base:
            i -= 1
        if n - i >= min_run:
            return i
    return n


def plan_from_call(call: dict) -> dict:
    """Extract the trimming plan from an `infer()` result."""
    fn = call.get("functional") or {}
    return {
        "status": call.get("status"),
        "trim_5p": fn.get("trim_5p", 0),
        "umi5": call.get("umi5_len", 0) or 0,
        # ordered 5' blocks: the UMI may be split around a barcode, in which case the
        # UMI bases are NOT the last umi5 nt before the footprint
        "p5_layout": call.get("p5_layout") or [],
        "umi3": _as_int(call.get("umi3_len"), 0),
        "nt3": _as_int(call.get("nt3_len"), 0),
        "bc3": call.get("barcode3_seq", "none"),
        "adapter_name": call.get("adapter3_name", "unknown"),
        "adapter_seq": call.get("adapter3_seq", ""),
        "polyA": call.get("polyA_tail", "none"),
        "fp_mode": call.get("footprint_len_mode") or 28,
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
        # find the ligated adapter and cut it, plus the 3' construct in front of it
        a = find_adapter(seq, p["adapter_seq"], min_start=max(s5, p["fp_mode"] - 6))
        adapter_found = a >= 0
        if a < 0:
            e3 = len(seq)                      # adapter not in this read: keep to end
        else:
            # umi3 + a degenerate non-templated block + barcode3 all sit between the
            # footprint and the adapter, and all of it comes off the footprint
            cons_len = p["umi3"] + p["nt3"] + bc3len
            umi3 = seq[a - p["umi3"] - bc3len:a - bc3len] if p["umi3"] else ""
            e3 = a - cons_len
    elif p["polyA"] not in ("none", "", None) or name == "none_visible":
        e3 = trim_polyA(seq)                   # trim the poly(A) tail (adapter beyond it)
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
               label: str = "", discard_untrimmed: bool = True) -> dict:
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
    p = plan_from_call(call)
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
    return stats

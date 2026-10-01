"""The six-panel ribo-seq QC figure, drawn from the dict `ribomine.qc.verdict.qc()`
returns (or the same read back from `Sample.qc_json`). Nothing is recomputed.

  A  read-length distribution, with the footprint window
  B  start-codon metagene (5' ends; a 3-nt comb with a P-site offset of ~12)
  C  stop-codon metagene
  D  in-frame fraction per read length (chance = 1/3)
  E  region composition, with the CDS floor
  F  verdict, reasons, contaminant and mapping stats, projected usable reads
"""
from __future__ import annotations

import logging
import os

import matplotlib

matplotlib.use("Agg")           # no display needed; safe in a process pool

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

LOG = logging.getLogger("ribomine.qc.plot")

C_IN = "#2c9e57"     # in-frame / CDS / good
C_OUT = "#b8b8b8"    # off-frame / neutral
C_HL = "#2f6fb0"     # highlight
C_BAD = "#c8443b"
REGION_COLORS = {
    "CDS": "#2c9e57", "5'UTR": "#7cc4f2", "3'UTR": "#2f6fb0", "ncRNA": "#8250b0",
    "intron": "#e0a81e", "intergenic": "#b8b8b8", "mito": "#c8443b",
}

# Fallback thresholds, equal to the `qc` config defaults. The figure draws the
# thresholds recorded in the QC dict (its `thresholds` block) when present, so it
# shows the values the verdict was made with.
FP_LEN_LO, FP_LEN_HI = 25, 36     # qc.footprint_len_lo / qc.footprint_len_hi
PERIODIC_MIN = 0.42               # qc.periodic_min      (chance = 1/3)
PERIODIC_STRONG = 0.50            # qc.periodic_strong
READ_LEN_MIN_FRAC = 0.75          # qc.read_len_min_frac (hard gate, read length)
CDS_REGION_MIN = 0.50             # qc.cds_region_min    (hard gate, region composition)

# metagene window drawn around the codon (nt, 5'-end based)
META_LO, META_HI = -30, 45


def plot_qc(qc: dict, out_path: str, *, dpi: int = 110) -> str:
    """Render the six-panel QC figure for one sample. Returns `out_path`.

    The extension of `out_path` picks the format (.png / .pdf / .svg).
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5))
    try:
        _plot_len(axes[0, 0], qc)
        _plot_metagene(axes[0, 1], qc.get("metagene_start", {}),
                       "B  start-codon metagene (5' ends)")
        _plot_metagene(axes[0, 2], qc.get("metagene_stop", {}),
                       "C  stop-codon metagene (5' ends)")
        _plot_frame_by_len(axes[1, 0], qc)
        _plot_regions(axes[1, 1], qc)
        _plot_verdict(axes[1, 2], qc)

        ok = bool(qc.get("is_riboseq"))
        ti = qc.get("verdict") == "TI-SEQ"
        fig.suptitle(f"{qc.get('label', '?')}    ribo-seq QC    "
                     f"[{qc.get('n_reads_scored', 0):,} reads scored]    "
                     f"→  {qc.get('verdict', '?')}",
                     fontsize=12.5, fontweight="bold",
                     color=(C_HL if ti else (C_IN if ok else C_BAD)), y=0.99)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    finally:
        # always close: figures otherwise accumulate across samples
        plt.close(fig)
    LOG.debug("wrote %s", out_path)
    return out_path


# --- helpers ----------------------------------------------------------------
def _thr(qc: dict, name: str, default: float) -> float:
    """The threshold recorded in the QC dict (`thresholds` block or top level),
    else `default`."""
    block = qc.get("thresholds")
    if isinstance(block, dict) and isinstance(block.get(name), (int, float)):
        return float(block[name])
    if isinstance(qc.get(name), (int, float)):
        return float(qc[name])
    return float(default)


def _int_keys(h: dict | None) -> dict[int, float]:
    """Histogram with int keys (they are strings after a JSON round-trip)."""
    return {int(k): v for k, v in (h or {}).items()}


def _densify(h: dict[int, float], lo: int, hi: int) -> tuple[list[int], list[float]]:
    xs = list(range(lo, hi + 1))
    return xs, [h.get(x, 0) for x in xs]


def _empty(ax, title: str, msg: str = "no reads") -> None:
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color=C_OUT)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")


# --- A: read-length distribution -------------------------------------------
def _plot_len(ax, qc: dict) -> None:
    title = "A  read-length distribution"
    h = _int_keys(qc.get("read_len_hist"))
    if not h:
        _empty(ax, title)
        return
    flo = int(_thr(qc, "footprint_len_lo", FP_LEN_LO))
    fhi = int(_thr(qc, "footprint_len_hi", FP_LEN_HI))
    lo, hi = max(min(h), 15), min(max(h), 45)
    xs, ys = _densify(h, lo, hi)
    tot = sum(h.values()) or 1
    fr = [y / tot for y in ys]
    ax.bar(xs, fr, width=0.85, color=[C_IN if flo <= x <= fhi else C_OUT for x in xs],
           edgecolor="white", linewidth=0.3)
    # shade the footprint window
    ax.axvspan(flo - 0.5, fhi + 0.5, color=C_IN, alpha=0.06)
    ax.set_xlabel("read length (nt)")
    ax.set_ylabel("fraction of reads")
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    # hard gate: the fraction of reads inside the window must reach this
    rl_floor = _thr(qc, "read_len_min_frac", READ_LEN_MIN_FRAC)
    peak = qc.get("read_len_peak_frac", 0) or 0
    ax.text(0.98, 0.95,
            f"mode {qc.get('read_len_mode', '?')} nt\n"
            f"{peak:.0%} in {flo}-{fhi} (need {rl_floor:.0%})",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            color=C_BAD if peak < rl_floor else "0.15",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))


# --- B/C: start- and stop-codon metagene ------------------------------------
def _plot_metagene(ax, meta: dict | None, title: str, lo: int = META_LO,
                   hi: int = META_HI) -> None:
    m = _int_keys(meta)
    if not m:
        _empty(ax, title, "no reads near codons")
        return
    xs, ys = _densify(m, lo, hi)
    peak = xs[int(np.argmax(ys))]                       # dominant offset (~-12)
    # highlight every third base counted from the peak (the P-site offset), not from 0
    cols = [C_HL if (x - peak) % 3 == 0 else C_OUT for x in xs]
    ax.bar(xs, ys, width=0.85, color=cols, edgecolor="white", linewidth=0.2)
    ax.axvline(0, color=C_BAD, lw=1.2, ls="--")
    ax.text(0.02, 0.95, f"P-site offset ~{-peak} nt" if peak < 0 else f"peak {peak}",
            transform=ax.transAxes, va="top", fontsize=8, color=C_HL)
    ax.set_xlabel("5' end offset from codon (nt)")
    ax.set_ylabel("reads")
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")


# --- D: periodicity by read length ------------------------------------------
def _plot_frame_by_len(ax, qc: dict) -> None:
    fbl = _int_keys(qc.get("frame_by_len"))
    pmin = _thr(qc, "periodic_min", PERIODIC_MIN)
    pstrong = _thr(qc, "periodic_strong", PERIODIC_STRONG)
    if fbl:
        Ls = sorted(fbl)
        fr = [fbl[L] for L in Ls]
        cols = [C_IN if f >= pstrong else (C_HL if f >= pmin else C_OUT) for f in fr]
        ax.bar(Ls, fr, width=0.85, color=cols, edgecolor="white", linewidth=0.3)
    ax.axhline(1 / 3, ls="--", color=C_BAD, lw=1)
    ax.text(0.98, 1 / 3, " chance (1/3)", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=7, color=C_BAD)
    ax.set_ylim(0, 1)
    ax.set_xlabel("read length (nt)")
    ax.set_ylabel("in-frame fraction")
    ax.set_title("D  periodicity by read length", loc="left", fontsize=10, fontweight="bold")


# --- E: region composition ---------------------------------------------------
def _plot_regions(ax, qc: dict) -> None:
    rf = qc.get("region_frac", {}) or {}
    order = ["CDS", "5'UTR", "3'UTR", "ncRNA", "intron", "intergenic", "mito"]
    items = [(r, rf.get(r, 0.0)) for r in order if rf.get(r, 0) > 0]
    labels = [r for r, _ in items]
    vals = [v for _, v in items]
    cols = [REGION_COLORS.get(r, "#888") for r in labels]
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, color=cols, edgecolor="white")
    for yi, v in zip(y, vals):
        ax.text(v + 0.01, yi, f"{v:.0%}", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    # hard gate: CDS must be at least this share of all reads; drawn on the CDS bar only
    cds_floor = _thr(qc, "cds_region_min", CDS_REGION_MIN)
    xhi = max((max(vals) if vals else 0.0), cds_floor) * 1.18
    if "CDS" in labels:
        cy = int(y[labels.index("CDS")])
        ax.plot([cds_floor, cds_floor], [cy - 0.45, cy + 0.45],
                color=C_BAD, lw=1.2, ls="--")
        ax.text(cds_floor, cy + 0.5, f"floor {cds_floor:.0%}",
                color=C_BAD, fontsize=7, ha="center", va="bottom")
    ax.set_xlim(0, xhi if xhi else 1)
    ax.set_xlabel("fraction of reads")
    ax.set_title("E  region composition", loc="left", fontsize=10, fontweight="bold")


# --- F: verdict, contaminants, mapping, usable reads -------------------------
def _wrap(s: str, width: int = 58) -> list[str]:
    out, line = [], ""
    for w in s.split():
        if len(line) + len(w) + 1 > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out


def _plot_verdict(ax, qc: dict) -> None:
    ax.axis("off")
    verdict = qc.get("verdict", "?")
    ti = verdict == "TI-SEQ"
    ok = bool(qc.get("is_riboseq"))
    color = C_HL if ti else (C_IN if ok else C_BAD)
    ax.add_patch(plt.Rectangle((0.02, 0.86), 0.96, 0.13, transform=ax.transAxes,
                               fc=color, alpha=0.15, ec=color, lw=1.5))
    mark = "✓" if ok else "✗"
    fs = 15 if len(verdict) < 12 else 12
    ax.text(0.5, 0.925, f"{verdict}  {mark}", transform=ax.transAxes, ha="center",
            va="center", fontsize=fs, fontweight="bold", color=color)

    lines: list[str] = []
    if qc.get("verdict_reason"):
        lines += _wrap("reason: " + qc["verdict_reason"])
        lines.append("")
    # Wrap the reasons too: an unwrapped long line widens the whole figure, because
    # savefig crops to the bounding box of all artists.
    for r in qc.get("reasons", []):
        wrapped = _wrap(r, width=56)
        lines.append(f"• {wrapped[0]}")
        lines += [f"  {w}" for w in wrapped[1:]]      # hanging indent

    cm = qc.get("contaminants") or {}
    if cm.get("n_sampled"):
        lines.append("")
        lines.append(f"contaminants removed (of {cm['n_sampled']:,} sampled):")
        lines.append(f"  rRNA/tRNA/etc:   {cm.get('frac_rRNA_tRNA_etc', 0):.0%}")
        lines.append(f"  low-complexity:  {cm.get('frac_low_complexity', 0):.0%}")
        if cm.get("frac_position_pileup") is not None:
            lines.append(f"  position pile-up: {cm.get('frac_position_pileup', 0):.0%}")
        lines.append(f"  kept:            {cm.get('frac_kept', 0):.0%}")

    m = qc.get("mapping", {}) or {}
    if m.get("n_input"):
        lines.append("")
        lines.append(f"reads into alignment: {m.get('n_input', 0):,}")
        if m.get("frac_unique") is not None:
            low = (qc.get("signals", {}) or {}).get("low_unique_mapping")
            lines.append(f"  uniquely mapped: {m['frac_unique']:.0%}"
                         + ("  <-- LOW" if low else ""))
        if m.get("frac_multimapping") is not None:
            lines.append(f"  multimapping:    {m['frac_multimapping']:.0%}")
        if m.get("frac_unmapped") is not None:
            lines.append(f"  unmapped:        {m['frac_unmapped']:.0%}")

    u = qc.get("usable")
    if u:
        lines.append("")
        lines.append(f"projected usable reads: {u['projected_usable_reads']:,}")
        lines.append(f"  ({u['usable_frac_of_sample']:.0%} of sample x "
                     f"{u['total_dataset_reads']:,} in dataset)")

    lines.append("")
    lines.append(f"periodicity: in-frame {qc.get('periodicity_inframe_frac', 0):.0%}, "
                 f"TVD {qc.get('periodicity_tvd_uniform', 0):.2f}")
    lines.append(f"CDS: {qc.get('cds_frac_of_genic', 0):.0%} of genic reads")
    lines.append(f"top 5' locus: {qc.get('top5p_locus_frac', 0):.0%} of reads")

    ax.text(0.04, 0.80, "\n".join(lines), transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5, family="monospace")
    ax.set_title("F  verdict & mapping", loc="left", fontsize=10, fontweight="bold")

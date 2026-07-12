"""The read-architecture figure: the positional quantities the caller reads, with
the inferred segmentation drawn on top.

Every panel is one of the signals `ribomine.arch.infer` actually decides on, so
the figure is a **visual audit of the call** rather than a separate summary:

  A  5' genome-match / penetrance   -- the CDF-of-clip-length curve; the step
                                       whose position is the 5' UMI length
  B  5' base composition + entropy  -- balanced=UMI, flat=barcode, A/T=RT nt,
                                       G/C=template-switch
  C  5' clip-length histogram       -- the survival function penetrance is read
                                       from; a fixed UMI is a spike at its length
  D  3' genome-match, adapter-anchored -- walking left from the adapter start
                                       through barcode / UMI into the footprint
  E  3' base composition + entropy  -- the 3' construct, same composition axes
  F  footprint-length histogram     -- genuine footprints spread ~26-34 nt; a
                                       single spike is a fixed-length artefact

The shaded spans and the dashed reference lines are exactly the thresholds the
caller used, so where it drew a boundary you can see why.

The call is passed in, never recomputed: the pipeline infers once and the figure
must show *that* call, thresholds included.
"""
from __future__ import annotations

import logging
import math
import os

import matplotlib

matplotlib.use("Agg")           # process-pool safe: no display, no GUI thread

import matplotlib.pyplot as plt              # noqa: E402
import numpy as np                           # noqa: E402
from matplotlib.lines import Line2D          # noqa: E402
from matplotlib.patches import Patch         # noqa: E402

from ribomine.arch import infer as IA        # noqa: E402
from ribomine.arch.infer import adapter_disp, architecture_string  # noqa: E402

LOG = logging.getLogger("ribomine.arch.plot")

# category colours, shared across the match/shading panels (A, D)
C_FOOT = "#2c9e57"   # footprint (genomic)
C_UMI = "#2f6fb0"    # random-templated (UMI / spacer)
C_BC = "#e08a1e"     # fixed-templated (barcode)
C_ADAP = "#7a7a7a"   # adapter
C_RT = "#c8443b"     # enzymatic RT base (kept)
C_TS = "#8250b0"     # template-switch / linker (trimmed)
C_MATCH = "#1a1a1a"
# nucleotide colours for the composition panels (B, E). These never share a
# panel with the category shading above, and each set has its own legend.
BASE_COLORS = {"A": "#3b8c3b", "C": "#2f6fb0", "G": "#e0a81e", "T": "#c8443b"}

# a composition column is only drawn where enough reads reach that position
MIN_COMP_N = 20
# below this many adapter-anchored reads the 3' construct cannot be read at all,
# whatever the anchored fraction says
MIN_ANCHOR_READS = 200


def plot_arch(profile: dict, call: dict, out_path: str, *, dpi: int = 110) -> str:
    """Render the six-panel architecture figure for one sample. Returns `out_path`.

    `call` is the dict `ribomine.arch.infer.infer()` returned for this profile --
    including an `undetermined` one, whose panels then show why the caller
    refused. The extension of `out_path` picks the format (.png / .pdf / .svg).
    """
    # an early-refused call carries no plateau; recompute the sample's own
    # genomic ceiling so the panels still have their reference line
    plateau = call.get("genomic_plateau") or IA.genomic_plateau(profile)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.4))
    try:
        _plot_5p_match(axes[0, 0], profile, call, plateau)
        _plot_5p_comp(axes[0, 1], profile, call)
        _plot_clip5(axes[0, 2], profile, call)
        _plot_3p_match(axes[1, 0], profile, call, plateau)
        _plot_3p_comp(axes[1, 1], profile, call)
        _plot_footprint(axes[1, 2], profile, call)

        # legend 1: nucleotide colours for the composition panels (B, E)
        base_leg = [Patch(fc=BASE_COLORS[b], label=b) for b in "ACGT"] + [
            Line2D([0], [0], color="black", marker="o", ms=4, lw=1.1, label="entropy")]
        l1 = fig.legend(handles=base_leg, loc="lower center", ncol=5, fontsize=8,
                        frameon=False, bbox_to_anchor=(0.5, 0.032),
                        title="composition panels (B, E) — base fraction")
        l1.get_title().set_fontsize(8)
        fig.add_artist(l1)

        # legend 2: segment categories for the shaded spans (panels A, D) and the call
        cat_leg = [
            Patch(fc=C_FOOT, alpha=0.4, label="footprint (genomic)"),
            Patch(fc=C_UMI, alpha=0.4, label="random-templated (UMI)"),
            Patch(fc=C_BC, alpha=0.4, label="fixed-templated (barcode)"),
            Patch(fc=C_RT, alpha=0.4, label="enzymatic RT nt (kept)"),
            Patch(fc=C_TS, alpha=0.4, label="template-switch (trimmed)"),
            Patch(fc=C_ADAP, alpha=0.4, label="adapter"),
        ]
        l2 = fig.legend(handles=cat_leg, loc="lower center", ncol=6, fontsize=8,
                        frameon=False, bbox_to_anchor=(0.5, -0.01),
                        title="segments & shading (panels A, D)")
        l2.get_title().set_fontsize(8)

        status = call.get("status", "?")
        color = C_FOOT if status == "ok" else C_RT
        conc = profile.get("top5p_locus_frac")
        conc_s = f", top-locus {100 * conc:.0f}%" if conc is not None else ""
        fig.suptitle(
            f"{profile.get('label', '?')}    read-architecture profile "
            f"[{profile.get('n_used', 0):,} reads, anchor {profile.get('anchor_kind', '?')} "
            f"{100 * profile.get('frac_anchored', 0):.0f}%{conc_s}]\n"
            f"{architecture_string(call)}",
            fontsize=11.5, fontweight="bold", color=color, y=0.995)
        ts = _trim_summary(call)
        if ts:
            fig.text(0.5, 0.925, ts, ha="center", va="top", fontsize=8.5, color="#444444")

        fig.tight_layout(rect=(0, 0.075, 1, 0.905))
        os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    finally:
        # hundreds of samples run through a process pool -- a leaked figure is a
        # leaked megabyte
        plt.close(fig)
    LOG.debug("wrote %s", out_path)
    return out_path


# --- thresholds -------------------------------------------------------------
# The config-driven thresholds are read off the CALL (the caller records the
# Thresholds it was built with), so the figure draws the lines the decision was
# actually made against; `infer`'s module defaults are only the fallback for a
# call that did not record them. CHANCE is a property of the alphabet, not a
# knob, so it is always the module constant.
_THR_DEFAULTS = IA.Thresholds()


def _thr(call: dict, name: str) -> float:
    block = call.get("thresholds")
    if isinstance(block, dict) and isinstance(block.get(name), (int, float)):
        return float(block[name])
    if isinstance(call.get(name), (int, float)):
        return float(call[name])
    return float(getattr(_THR_DEFAULTS, name))


# --- small helpers ----------------------------------------------------------
def _as_int(v, d=None):
    try:
        return int(v)
    except (TypeError, ValueError):
        return d


def _int_keys(h: dict | None) -> dict[int, float]:
    """Histogram keys are ints in memory and strings after a JSON round-trip;
    the figure is drawn from either."""
    return {int(k): v for k, v in (h or {}).items()}


def _densify(h: dict[int, float], lo: int, hi: int) -> tuple[list[int], list[float]]:
    xs = list(range(lo, hi + 1))
    return xs, [h.get(x, 0) for x in xs]


def _entropy(row) -> float:
    r = np.asarray(row, dtype=float)
    if np.isnan(r).any() or r.sum() <= 0:
        return float("nan")
    r = r[r > 0]
    return float(-(r * np.log2(r)).sum())


def _shade(ax, x0, x1, color, label=None, alpha=0.16) -> None:
    ax.axvspan(x0, x1, color=color, alpha=alpha, lw=0, label=label)


def _adapter_usable(profile: dict, call: dict) -> bool:
    """Is the 3' side readable in adapter-anchored coordinates?

    A homopolymer anchor (poly-A / poly-G) is a tail, not the adapter, so it
    anchors nothing; and a construct only becomes visible once enough reads are
    long enough to expose the adapter at all.
    """
    kind = profile.get("anchor_kind", "none")
    anchored = (kind.startswith("panel:") and kind.split(":", 1)[1] not in ("polyA", "polyG")
                or kind == "denovo")
    return bool(anchored
                and profile.get("n_anchored", 0) >= MIN_ANCHOR_READS
                and profile.get("frac_anchored", 0) >= _thr(call, "min_anchor_frac"))


def _stacked_composition(ax, comp, n, positions, xlabels, call, bases="ACGT",
                         min_n=MIN_COMP_N):
    """Draw an A/C/G/T stacked bar per position and overlay Shannon entropy."""
    ent_random = _thr(call, "ent_random")
    ent_const = _thr(call, "ent_const")
    ent = []
    ax2 = ax.twinx()
    for k, (idx, _xl) in enumerate(zip(positions, xlabels)):
        if idx < 0 or idx >= len(comp):
            ent.append(float("nan"))
            continue
        row = comp[idx]
        if any(math.isnan(v) for v in row) or (n and n[idx] < min_n):
            ent.append(float("nan"))
            continue
        bottom = 0.0
        for bi, b in enumerate(bases):
            v = row[bi]
            ax.bar(k, v, bottom=bottom, width=0.86, color=BASE_COLORS[b],
                   edgecolor="white", linewidth=0.3)
            bottom += v
        ent.append(_entropy(row))
    ax.set_ylim(0, 1)
    ax.set_ylabel("base fraction")
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels, fontsize=7)
    # entropy overlay: the two dotted lines are the random / constant calls
    xs = [k for k, e in enumerate(ent) if not math.isnan(e)]
    es = [e for e in ent if not math.isnan(e)]
    ax2.plot(xs, es, "o-", color="black", ms=3, lw=1.1, label="entropy")
    ax2.axhline(ent_random, ls=":", color="black", lw=0.8)
    ax2.axhline(ent_const, ls=":", color="black", lw=0.8)
    ax2.set_ylim(0, 2.05)
    ax2.set_ylabel("entropy (bits)")
    ax2.text(len(xlabels) - 0.5, ent_random + 0.02, "random", fontsize=6,
             ha="right", va="bottom", color="black")
    ax2.text(len(xlabels) - 0.5, ent_const - 0.02, "const", fontsize=6,
             ha="right", va="top", color="black")
    return ax2


# --- A: 5' genome-match & penetrance ----------------------------------------
def _plot_5p_match(ax, prof: dict, call: dict, plateau) -> None:
    m5 = np.asarray(prof["p5_match"], dtype=float)
    npos = min(17, len(m5))
    x = np.arange(npos)
    ax.plot(x, m5[:npos], "o-", color=C_MATCH, ms=4, lw=1.4, label="genome match", zorder=5)
    ax.axhline(IA.CHANCE, ls="--", color=C_UMI, lw=1, label="chance (0.25)")
    if plateau:
        ax.axhline(plateau, ls="--", color=C_FOOT, lw=1,
                   label=f"genomic plateau ({plateau:.2f})")
        # penetrance = survival function of the 5' construct length, measured
        # against the sample's OWN genomic ceiling rather than against 1.0
        pen = np.clip((plateau - m5[:npos]) / max(plateau - IA.CHANCE, 1e-6), 0, 1)
        ax.plot(x, pen, "-", color="#b0b0b0", lw=1.0, label="penetrance", zorder=3)

    umi5 = _as_int(call.get("umi5_len"))
    fn = call.get("functional") or {}
    trim5 = fn.get("trim_5p") if fn else None
    rt5 = fn.get("footprint_retains_rt_nt", 0) if fn else 0

    if trim5 is not None:
        # construct region (trimmed), and the footprint-start boundary
        _shade(ax, -0.5, trim5 - 0.5, "#d9d9d9", alpha=0.5)
        ax.axvline(trim5 - 0.5, color=C_FOOT, lw=1.6, zorder=6)
        ax.text(trim5 - 0.4, 0.02, "footprint start\n(trim boundary)", fontsize=6.5,
                color=C_FOOT, va="bottom", ha="left")
        # the construct blocks, each at its own offset. The UMI is not always one
        # contiguous stretch ending at the footprint: a sample barcode can sit
        # BETWEEN two UMI blocks ([umi5 5][ATTGGC][umi5 4], SRR26942052), and
        # shading it as a single 9-nt UMI would paint over the barcode.
        lay = call.get("p5_layout") or []
        if lay:
            for blk in lay:
                c = C_BC if blk["role"] == "barcode5" else C_UMI
                _shade(ax, blk["offset"] - 0.5, blk["offset"] + blk["len"] - 0.5, c, alpha=0.22)
        elif umi5:
            u1 = trim5 - 0.5
            u0 = max(-0.5, trim5 - umi5 - 0.5)
            _shade(ax, u0, u1, C_UMI, alpha=0.22)
        if rt5:
            _shade(ax, trim5 - 0.5, trim5 + rt5 - 0.5, C_RT, alpha=0.3)

    ax.set_xticks(range(0, npos, 2))
    ax.set_xlabel("read position (0 = first base)")
    ax.set_ylabel("fraction")
    ax.set_ylim(0, 1.02)
    ax.set_xlim(-0.6, npos - 0.4)
    ax.set_title("A  5' genome-match & penetrance", loc="left", fontsize=10, fontweight="bold")
    txt = f"UMI = {umi5 if umi5 else 0} nt"
    if rt5:
        txt += f",  RT +{rt5} nt (kept)"
    ax.text(0.98, 0.06, txt, transform=ax.transAxes, fontsize=8, ha="right",
            va="bottom", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
    ax.legend(fontsize=6.5, loc="center right", framealpha=0.85)


# --- B: 5' composition & entropy ---------------------------------------------
def _plot_5p_comp(ax, prof: dict, call: dict):
    positions = list(range(0, 16))
    xlabels = [str(p) for p in positions]
    ax2 = _stacked_composition(ax, prof["p5_comp"], prof.get("p5_n"), positions, xlabels, call)
    fn = call.get("functional") or {}
    trim5 = fn.get("trim_5p")
    if trim5 is not None and trim5 <= len(positions):
        ax.axvline(trim5 - 0.5, color=C_FOOT, lw=1.6)
    ax.set_xlabel("read position")
    ax.set_title("B  5' composition & entropy", loc="left", fontsize=10, fontweight="bold")
    return ax2


# --- C: 5' clip-length histogram ---------------------------------------------
def _plot_clip5(ax, prof: dict, call: dict) -> None:
    xs, ys = _densify(_int_keys(prof.get("clip5_hist")), 0, 15)
    tot = sum(ys) or 1
    fr = [y / tot for y in ys]
    ax.bar(xs, fr, width=0.86, color=C_UMI, alpha=0.85, edgecolor="white", linewidth=0.3)
    umi5 = _as_int(call.get("umi5_len"), 0)
    if umi5:
        ax.axvline(umi5, color=C_FOOT, lw=1.4, ls="-")
        ax.text(umi5 + 0.5, max(fr) * 0.92, f"UMI len = {umi5}", fontsize=7.5,
                color=C_FOOT, va="top", ha="left")
    ax.set_xlabel("5' soft-clip length (nt)")
    ax.set_ylabel("fraction of reads")
    ax.set_title("C  5' clip-length histogram (penetrance source)", loc="left",
                 fontsize=10, fontweight="bold")


# --- D: 3' genome-match -------------------------------------------------------
def _adap_window(call: dict) -> list[int]:
    """k-offsets to draw around the anchor, in adapter-anchored coordinates.

    The window must reach left of the *inferred adapter start*, not of the raw
    anchor: a de-novo adapter can begin well inside the panel (SRR27968072 starts
    at A-18), and a fixed [-14, +6] window then shows nothing but adapter -- the
    construct it is meant to evidence sits off-panel entirely.
    """
    start_rel = int(call.get("adapter3_start_rel", 0) or 0)   # <= 0
    umi3 = _as_int(call.get("umi3_len"), 0) or 0
    nt3 = _as_int(call.get("nt3_len"), 0) or 0
    bc3 = call.get("barcode3_seq", "none")
    clen = umi3 + nt3 + (len(bc3) if bc3 not in ("none", "unknown", "", None) else 0)
    lo = min(-14, start_rel - clen - 4)      # construct + a few genomic bases of context
    return list(range(lo, 7))


def _plot_3p_match(ax, prof: dict, call: dict, plateau) -> None:
    if _adapter_usable(prof, call):
        _plot_3p_adapter(ax, prof, call, plateau)
    else:
        _plot_3p_readend(ax, prof, call, plateau, prof.get("anchor_kind", "none"))


def _plot_3p_adapter(ax, prof: dict, call: dict, plateau) -> None:
    a0 = prof["anchor_offset"]
    match = np.asarray(prof["adap_match"], dtype=float)
    ks = _adap_window(call)
    x = np.arange(len(ks))
    y = [match[a0 + k] if 0 <= a0 + k < len(match) else np.nan for k in ks]
    ax.plot(x, y, "o-", color=C_MATCH, ms=4, lw=1.4, zorder=5, label="genome match")
    ax.axhline(IA.CHANCE, ls="--", color=C_UMI, lw=1, label="chance")
    if plateau:
        ax.axhline(plateau, ls="--", color=C_FOOT, lw=1, label=f"plateau ({plateau:.2f})")
        # the two decision lines, both expressed as a fraction of the way from
        # chance up to this sample's own plateau
        thr_g = IA.CHANCE + _thr(call, "genomic_frac") * (plateau - IA.CHANCE)
        thr_u = IA.CHANCE + _thr(call, "umi_match_frac") * (plateau - IA.CHANCE)
        ax.axhline(thr_g, ls=":", color=C_FOOT, lw=0.9)
        ax.axhline(thr_u, ls=":", color=C_UMI, lw=0.9)
        ax.text(len(ks) - 1, thr_g + 0.01, "genomic thr", fontsize=6, ha="right", color=C_FOOT)
        ax.text(len(ks) - 1, thr_u - 0.01, "UMI thr", fontsize=6, ha="right", va="top",
                color=C_UMI)

    a0k = ks.index(0)
    # the caller may extend the adapter left of the anchor (the de-novo seed can
    # start a few nt inside the adapter), so mark the *inferred* adapter start,
    # not the anchor. adapter3_start_rel is that offset relative to A0 (<= 0).
    start_rel = int(call.get("adapter3_start_rel", 0) or 0)
    astart_k = max(0, min(len(ks) - 1, a0k + start_rel))
    # Boundaries sit BETWEEN points (at cell edges, x-0.5), never through them: a
    # line drawn on a point hides the very base it is separating. Position k spans
    # [k-0.5, k+0.5], so "the adapter starts AT astart_k" is the edge astart_k-0.5.
    if start_rel:   # keep a faint marker at the raw anchor for reference
        ax.axvline(a0k - 0.5, color=C_ADAP, lw=0.8, ls=":", alpha=0.6, zorder=5)
    ax.axvline(astart_k - 0.5, color=C_ADAP, lw=1.6, zorder=6)
    ax.text(astart_k - 0.4, 0.02, "adapter start", fontsize=6.5, color=C_ADAP, va="bottom")
    _shade(ax, astart_k - 0.5, len(ks) - 0.5, C_ADAP, alpha=0.18)

    # construct (umi3 + non-templated + barcode3) sits between footprint and adapter
    umi3 = _as_int(call.get("umi3_len"), 0) or 0
    nt3 = _as_int(call.get("nt3_len"), 0) or 0
    bc3 = call.get("barcode3_seq", "none")
    bc3len = len(bc3) if bc3 not in ("none", "unknown", "", None) else 0
    clen = umi3 + nt3 + bc3len
    if clen:
        # Shade each block in ITS OWN category colour rather than painting the whole
        # construct one hue. The 3' construct is laid out insert-first --
        # [umi3][nontemplated3][barcode3][adapter] -- so a library with both a UMI and
        # a barcode (the McGlincy-Ingolia design: 5 nt random, then AGCTA) would
        # otherwise show its barcode shaded as UMI, which is precisely the distinction
        # the figure exists to let you check.
        edge = astart_k - clen - 0.5
        for span, colour in ((umi3, C_UMI), (nt3, C_TS), (bc3len, C_BC)):
            if span:
                _shade(ax, edge, edge + span, colour, alpha=0.22)
                edge += span
        ax.axvline(astart_k - clen - 0.5, color=C_FOOT, lw=1.6, zorder=6)
        ax.text(astart_k - clen - 0.4, 0.02, "footprint end", fontsize=6.5,
                color=C_FOOT, va="bottom", ha="left")
    ax.set_xticks(x)
    ax.set_xticklabels([f"A{k:+d}" if k else "A0" for k in ks], fontsize=6.5, rotation=90)
    ax.set_ylabel("fraction")
    ax.set_ylim(0, 1.02)
    ax.set_title("D  3' genome-match (adapter-anchored)", loc="left", fontsize=10,
                 fontweight="bold")
    ax.legend(fontsize=6.5, loc="center left", framealpha=0.85)
    lay = f"umi3={umi3}, bc3={bc3 if bc3len else 'none'}, adapter={call.get('adapter3_name')}"
    ax.text(0.98, 0.98, lay, transform=ax.transAxes, fontsize=7, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))


def _plot_3p_readend(ax, prof: dict, call: dict, plateau, kind: str) -> None:
    t3 = np.asarray(prof["t3_match"], dtype=float)
    npos = min(14, len(t3))
    x = np.arange(npos)
    ax.plot(x, t3[:npos], "o-", color=C_MATCH, ms=4, lw=1.4, label="genome match")
    ax.axhline(IA.CHANCE, ls="--", color=C_UMI, lw=1, label="chance")
    if plateau:
        ax.axhline(plateau, ls="--", color=C_FOOT, lw=1, label=f"plateau ({plateau:.2f})")
    ax.set_xlabel("distance from read 3' end (0 = last base)")
    ax.set_ylabel("fraction")
    ax.set_ylim(0, 1.02)
    ax.set_title("D  3' genome-match (read-end-anchored; no adapter)", loc="left",
                 fontsize=10, fontweight="bold")
    ax.text(0.98, 0.98, f"anchor: {kind}\n(adapter not usable)", transform=ax.transAxes,
            fontsize=7, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
    ax.legend(fontsize=6.5, loc="center left", framealpha=0.85)


# --- E: 3' composition & entropy ---------------------------------------------
def _plot_3p_comp(ax, prof: dict, call: dict) -> None:
    if _adapter_usable(prof, call):
        a0 = prof["anchor_offset"]
        ks = _adap_window(call)
        positions = [a0 + k for k in ks]
        xlabels = [f"A{k:+d}" if k else "A0" for k in ks]
        _stacked_composition(ax, prof["adap_comp"], prof.get("adap_n"), positions, xlabels, call)
        # the adapter starts at the LEFT EDGE of its first base, not through it
        start_rel = int(call.get("adapter3_start_rel", 0) or 0)
        astart_k = max(0, min(len(ks) - 1, ks.index(0) + start_rel))
        ax.axvline(astart_k - 0.5, color=C_ADAP, lw=1.6)
        for lab in ax.get_xticklabels():
            lab.set_rotation(90)
            lab.set_fontsize(6.5)
        ax.set_xlabel("adapter-anchored position")
        ax.set_title("E  3' composition & entropy (adapter-anchored)", loc="left",
                     fontsize=10, fontweight="bold")
    else:
        positions = list(range(0, 12))
        xlabels = [str(p) for p in positions]
        _stacked_composition(ax, prof["t3_comp"], prof.get("t3_n"), positions, xlabels, call)
        ax.set_xlabel("distance from read 3' end")
        ax.set_title("E  3' composition & entropy (read-end-anchored)", loc="left",
                     fontsize=10, fontweight="bold")


# --- F: footprint-length distribution ----------------------------------------
def _plot_footprint(ax, prof: dict, call: dict) -> None:
    fph = _int_keys(prof.get("footprint_len_hist"))
    if fph:
        gate = _thr(call, "footprint_uniform_max")
        lo, hi = max(min(fph), 15), min(max(fph), 45)
        xs, ys = _densify(fph, lo, hi)
        tot = sum(fph.values()) or 1
        fr = [y / tot for y in ys]
        mode = max(fph, key=lambda k: fph[k])
        mode_fr = fph[mode] / tot
        # a single length above the artefact gate is not a footprint distribution
        # at all -- adapter dimers / a fixed contaminant -- and the caller refuses
        colors = [C_RT if x == mode and mode_fr > gate else C_FOOT for x in xs]
        ax.bar(xs, fr, width=0.86, color=colors, alpha=0.85, edgecolor="white", linewidth=0.3)
        ax.axhline(gate, ls="--", color=C_RT, lw=1)
        ax.text(hi, gate, f" artefact gate ({gate:.2f})",
                fontsize=6.5, color=C_RT, va="bottom", ha="right")
        ax.text(0.02, 0.98, f"mode {mode} nt ({mode_fr:.0%})", transform=ax.transAxes,
                fontsize=8, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
    ax.set_xlabel("footprint length (nt)")
    ax.set_ylabel("fraction of reads")
    ax.set_title("F  footprint-length distribution", loc="left", fontsize=10, fontweight="bold")


# --- the subtitle -------------------------------------------------------------
def _trim_summary(call: dict) -> str:
    """The load-bearing processing quantities, one compact line."""
    fn = call.get("functional") or {}
    if call.get("status") != "ok" or not fn:
        return ""
    rt = fn.get("footprint_retains_rt_nt", 0)
    rt_s = f" (keep {rt}nt RT)" if rt else ""
    adap = call.get("adapter3_name")
    if adap in ("none", "unknown", None):
        adap_s = "3' adapter already trimmed"
    elif adap == "none_visible":
        adap_s = "trim poly(A) tail (adapter beyond it, not identified)"
    else:
        adap_s = f"remove {adapter_disp(adap)} off 3'"
    return (f"trim {fn.get('trim_5p', 0)}nt off 5'{rt_s}   ·   {adap_s}   ·   "
            f"dedup UMI {fn.get('dedup_umi_len', 0)}nt")

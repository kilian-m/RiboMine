"""Merge the per-node results and report whether the run is complete.

    merge_results.py -c CONFIG --shards N

Called by `slurm/merge.sh` after the run array, however the array ended.

Merged view: each node wrote into its own `shards/work_NN/`. This builds
`<root>/merged/`, a workdir of symlinks into all of them, and runs RiboMine's own
report writers on it. The cohort tables (`qc/qc_summary.tsv`,
`architecture/architecture.tsv`, `mapping_summary.tsv`, `counts/gene_counts.tsv`)
are therefore written by the same code as on a single machine. The shards'
`failed.tsv` files are concatenated.

Status: what remains is decided from the per-sample JSONs, not from SLURM exit
codes, and written to `logs/status.txt` as COMPLETE or RESUBMIT NEEDED. A compute
node cannot call `sbatch`, so the next round is submitted by hand.
"""
from __future__ import annotations

import argparse
import copy
import os
import shutil
import sys

from ribomine import config as cfgmod, reports
from ribomine.utils import Sample, nonempty, read_json, read_tsv, setup_logging, write_tsv

from prepare_run import shard_list, shard_workdir


# ---------------------------------------------------------------------------
def owners(root: str, n: int) -> dict[str, int]:
    """acc -> the shard that owns it, from the shard lists the nodes read."""
    out: dict[str, int] = {}
    for k in range(n):
        p = shard_list(root, k)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for ln in fh:
                acc = ln.strip()
                if acc and not acc.startswith("#"):
                    out[acc] = k
    return out


def build_view(root: str, own: dict[str, int]) -> str:
    """Build `<root>/merged`: one workdir of symlinks into the shard workdirs.

    Rebuilt from scratch every round, so that no stale symlink survives.
    """
    merged = os.path.join(root, "merged")
    for sub in ("samples", "bams", "qc/plots", "architecture/plots", "meta"):
        shutil.rmtree(os.path.join(merged, sub), ignore_errors=True)
    for sub in ("samples", "bams", "qc/plots", "architecture/plots", "meta"):
        os.makedirs(os.path.join(merged, sub), exist_ok=True)

    _link(os.path.join(root, "refs"), os.path.join(merged, "refs"))
    _link(os.path.join(root, "meta", "candidates.tsv"),
          os.path.join(merged, "meta", "candidates.tsv"))

    for acc, k in own.items():
        wd = shard_workdir(root, k)
        _link(os.path.join(wd, "samples", acc),
              os.path.join(merged, "samples", acc))
        for name, dst in (
            (os.path.join("bams", f"{acc}.bam"), os.path.join("bams", f"{acc}.bam")),
            (os.path.join("bams", f"{acc}.bam.bai"), os.path.join("bams", f"{acc}.bam.bai")),
        ):
            _link(os.path.join(wd, name), os.path.join(merged, dst))
        for sub, suffix in (("qc/plots", ".qc"), ("architecture/plots", ".arch")):
            for ext in ("png", "pdf", "svg"):
                f = f"{acc}{suffix}.{ext}"
                _link(os.path.join(wd, sub, f), os.path.join(merged, sub, f))

    _merge_failures(root, merged, own)
    return merged


def _link(target: str, link: str) -> None:
    """Symlink `link` -> `target`, if the target exists at all."""
    if not os.path.exists(target):
        return
    if os.path.lexists(link):
        os.unlink(link)
    os.symlink(os.path.abspath(target), link)


def _merge_failures(root: str, merged: str, own: dict[str, int]) -> None:
    """Concatenate the shards' `failed.tsv` into one. `reports` reads it to say
    why a run has no QC verdict."""
    rows: list[dict] = []
    for k in sorted(set(own.values())):
        p = os.path.join(shard_workdir(root, k), "failed.tsv")
        if nonempty(p):
            try:
                rows += read_tsv(p)
            except OSError:
                continue
    dst = os.path.join(merged, "failed.tsv")
    if rows:
        write_tsv(dst, rows, ["run_accession", "stage", "error"])
    elif os.path.exists(dst):
        os.remove(dst)


# ---------------------------------------------------------------------------
def state(cfg, accs: list[str]) -> dict:
    """Count what is done, what is still to do and what failed.

    Read from the per-sample JSONs. "To do" follows the pipeline's gates: every
    run gets a QC verdict; only a kept verdict gets an architecture call; only a
    run with an architecture call is downloaded and mapped. A run rejected by QC
    counts as finished.
    """
    end = cfg["pipeline.end"]
    keep = set(cfg["pipeline.keep_verdicts"])
    undet_ok = bool(cfg["architecture.process_undetermined"])

    todo: list[str] = []
    done = {"qc": 0, "architecture": 0, "bam": 0}
    dropped = {"verdict": 0, "architecture": 0}

    for acc in accs:
        s = Sample(acc, cfg.workdir)
        q = read_json(s.qc_json)
        if q is None:
            todo.append(acc)
            continue
        done["qc"] += 1
        if end == "qc":
            continue
        if q.get("verdict") not in keep:
            dropped["verdict"] += 1          # QC said no: correctly finished with it
            continue

        call = read_json(s.arch_json)
        if call is None:
            todo.append(acc)
            continue
        done["architecture"] += 1
        if end == "architecture":
            continue
        if call.get("status") != "ok" and not undet_ok:
            dropped["architecture"] += 1     # no architecture, no trim: also finished
            continue

        from ribomine.pipeline import is_processed
        if is_processed(cfg, s):
            done["bam"] += 1
        else:
            todo.append(acc)

    fails = reports._failures(cfg)
    stuck = sorted(a for a in todo if a in fails)
    return {"n_total": len(accs), "done": done, "dropped": dropped,
            "todo": sorted(todo), "stuck": stuck, "end": end,
            "undetermined_processed": undet_ok}


def write_status(root: str, st: dict, cfg_path: str) -> str:
    """Write `logs/status.txt`: the state of the cohort and, if anything is left,
    the command for the next round."""
    todo, stuck = st["todo"], st["stuck"]
    retryable = [a for a in todo if a not in set(stuck)]
    complete = not retryable

    L = []
    L.append("COMPLETE" if complete else "RESUBMIT NEEDED")
    L.append("")
    L.append(f"end point            : {st['end']}")
    L.append(f"runs in the cohort   : {st['n_total']}")
    L.append(f"  QC verdict         : {st['done']['qc']}")
    if st["end"] != "qc":
        L.append(f"  architecture call  : {st['done']['architecture']}")
        L.append(f"  dropped by QC      : {st['dropped']['verdict']} "
                 f"(not ribo-seq / too low quality -- finished, not missing)")
    if st["end"] == "bam":
        L.append(f"  BAM                : {st['done']['bam']}")
        if not st["undetermined_processed"]:
            L.append(f"  dropped, no arch   : {st['dropped']['architecture']} "
                     f"(architecture undetermined -- set architecture."
                     f"process_undetermined=true to trim them best-effort)")
    L.append(f"  still to do        : {len(retryable)}")
    L.append(f"  tried and failed   : {len(stuck)}"
             + ("  <-- these will fail again; see failed.tsv" if stuck else ""))
    L.append("")

    if complete:
        L.append("Nothing left to run. The cohort tables are under:")
        L.append(f"  {os.path.join(root, 'merged')}")
    else:
        L.append("Resubmit from a LOGIN node (a compute node may not sbatch, so this")
        L.append("job could not do it). Finished runs are skipped, so each round is")
        L.append("shorter than the last:")
        L.append("")
        L.append(f"    slurm/master.sh {cfg_path} run")
        L.append("")
        if stuck:
            L.append("Failing runs, and why (they are not retried into eternity by")
            L.append("themselves -- read failed.tsv and decide):")
            for a in stuck[:20]:
                L.append(f"    {a}")
            if len(stuck) > 20:
                L.append(f"    ... and {len(stuck) - 20} more")

    path = os.path.join(root, "logs", "status.txt")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\n-> {path}")
    return path


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-c", "--config", required=True)
    p.add_argument("--shards", type=int, required=True)
    a = p.parse_args(argv)

    cfg = cfgmod.load(a.config)
    root = cfg.workdir
    setup_logging("INFO", os.path.join(cfg.dir("logs"), "merge.log"))

    own = owners(root, a.shards)
    if not own:
        print(f"no shard lists under {root}/shards -- has prep run?", file=sys.stderr)
        return 1
    accs = sorted(own)
    print(f"{len(accs)} run(s) across {a.shards} shard(s)")

    merged = build_view(root, own)
    print(f"merged view: {merged}")

    # the same Config, pointed at the merged view, so the report writers see one
    # workdir holding every sample
    mcfg = cfgmod.Config(data=copy.deepcopy(cfg.data), path=cfg.path)
    mcfg.data["project"]["workdir"] = merged

    end = mcfg["pipeline.end"]
    print(reports.qc_tsv(mcfg, accs))
    if end in ("architecture", "bam"):
        print(reports.arch_tsv(mcfg, accs))
    if end == "bam":
        print(reports.process_tsv(mcfg, accs))
        print(reports.counts_tsv(mcfg, accs) or "(no gene-count matrix)")

    print("\n" + reports.summary_line(mcfg, accs) + "\n")
    write_status(root, state(mcfg, accs), a.config)
    return 0


if __name__ == "__main__":
    sys.exit(main())

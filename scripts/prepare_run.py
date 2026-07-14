"""Split a cohort across the nodes of the SLURM run, once, and stably.

This is the small job (`slurm/prep.sh`) that the big one depends on. It does the
three things that must happen exactly once, in one process, before eight nodes
start at the same time:

1. **The query.** `ribomine query` searches ENA + NCBI and writes
   `<root>/meta/candidates.tsv`. Eight nodes each running their own copy of that
   search would be eight identical requests and eight different answers.

2. **The split.** The candidates are packed into one accession list per node.
   Packing is by `read_count`, longest-first (LPT) -- a run is 1-10 GB, and a
   node that draws the deep ones sets the wall time for the whole job.

3. **The per-node workdirs.** Each node gets its own `<root>/shards/work_NN`,
   because `ribomine run` writes cohort-level tables (`qc_summary.tsv`,
   `failed.tsv`, ...) at fixed paths under its workdir: eight nodes sharing one
   workdir would overwrite each other's. What they *may* share is seeded as a
   symlink -- `refs/` (the GTF index and the bowtie2 contaminant index, built
   once by `ribomine setup`) and `meta/candidates.tsv`.

**The split is pinned.** Re-running this after the query has grown does not
reshuffle: an accession that already has a shard keeps it, and only the new ones
are packed into the least-loaded nodes. Moving a finished accession to another
node would strand its results in a workdir nobody looks in any more, and its BAM
would be downloaded and mapped a second time.
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys

from ribomine import config as cfgmod
from ribomine.utils import read_tsv, setup_logging


def shard_dir(root: str) -> str:
    return os.path.join(root, "shards")


def shard_list(root: str, k: int) -> str:
    return os.path.join(shard_dir(root), f"shard_{k:02d}.txt")


def shard_workdir(root: str, k: int) -> str:
    return os.path.join(shard_dir(root), f"work_{k:02d}")


def manifest_path(root: str) -> str:
    return os.path.join(shard_dir(root), "manifest.json")


# ---------------------------------------------------------------------------
def accessions(cfg, override: str | None) -> tuple[list[str], dict[str, int]]:
    """The cohort, and each run's read count (the packing weight).

    Where the cohort comes from is `pipeline.start`, the same as it would be for a
    plain `ribomine run` -- so a config that names its own `accession_list` (a
    pilot, a re-run of a curated set) is split as it stands, and is NOT sent to
    search the archive behind its own back.

    A run whose read count ENA does not give gets the cohort median rather than
    zero: an unknown-size run is an average-size run, not a free one.
    """
    start = cfg["pipeline.start"]
    if start == "fastq":
        raise SystemExit(
            "pipeline.start='fastq' is not supported on the cluster: the split is over "
            "run accessions, and a directory of local FASTQs has none. Run those with "
            "`ribomine run --fastq-dir` on a single node instead.")

    if start == "accessions" and not override and not cfg["pipeline.accession_list"]:
        raise SystemExit("pipeline.start='accessions' needs pipeline.accession_list "
                         "(or pass --accessions)")

    src = override or (cfg._abs(cfg["pipeline.accession_list"])
                       if start == "accessions" else None)
    if src:
        if not os.path.isfile(src):
            raise SystemExit(f"accession list not found: {src}")
        with open(src) as fh:
            accs = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
        print(f"accessions: {len(accs)} from {src}  (pipeline.start={start})")
        # no ENA metadata for a bare list, so every run packs at the same weight
        return accs, {}

    tsv = os.path.join(cfg.dir("meta"), "candidates.tsv")
    if not (os.path.exists(tsv) and os.path.getsize(tsv) > 0):
        from ribomine.sra import query

        cfg.validate(need_reference=False, need_inputs=False)
        tsv = query.run_query(cfg)
    else:
        print(f"reusing the existing query: {tsv}  (delete it to search again)")

    rows = read_tsv(tsv)
    accs = [r["run_accession"] for r in rows if r.get("run_accession")]
    weights: dict[str, int] = {}
    for r in rows:
        try:
            weights[r["run_accession"]] = int(r["read_count"])
        except (KeyError, TypeError, ValueError):
            pass
    print(f"accessions: {len(accs)} from {tsv} "
          f"({len(weights)} with a read count)")
    return accs, weights


def pack(accs: list[str], weights: dict[str, int], n: int,
         pinned: dict[str, int]) -> list[list[str]]:
    """LPT bin-packing into `n` shards, honouring `pinned` (acc -> shard).

    Longest-processing-time-first: place the heaviest run on the lightest node,
    repeatedly. It is the standard 4/3-approximation, and here it is close to
    optimal because no single run is a meaningful fraction of a node's load --
    which is exactly the property PRICE2's stage-4 locus partition does *not*
    have, and why that one is imbalanced however it is packed.
    """
    default = int(statistics.median(weights.values())) if weights else 1
    bins: list[list[str]] = [[] for _ in range(n)]
    load = [0] * n

    for acc, k in sorted(pinned.items()):
        if 0 <= k < n and acc in set(accs):
            bins[k].append(acc)
            load[k] += weights.get(acc, default)

    fresh = [a for a in accs if a not in pinned]
    # -weight first, then the accession, so the packing is reproducible
    for acc in sorted(fresh, key=lambda a: (-weights.get(a, default), a)):
        k = min(range(n), key=lambda b: (load[b], b))
        bins[k].append(acc)
        load[k] += weights.get(acc, default)

    if fresh:
        lo, hi = min(load), max(load)
        print(f"packed {len(fresh)} new run(s) into {n} shard(s); "
              f"load spread {hi / max(lo, 1):.2f}x "
              f"(sizes {min(len(b) for b in bins)}-{max(len(b) for b in bins)} runs)")
    return bins


def existing_pins(root: str, n: int) -> dict[str, int]:
    """acc -> shard, from the shard lists already on disk. The lists are the
    truth, not the manifest: they are what the nodes actually read."""
    pins: dict[str, int] = {}
    for k in range(n):
        p = shard_list(root, k)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for ln in fh:
                acc = ln.strip()
                if acc and not acc.startswith("#"):
                    pins[acc] = k
    return pins


def seed_workdir(root: str, k: int, cfg) -> str:
    """One node's workdir: its own everything, except what may safely be shared.

    `refs/` and `meta/candidates.tsv` are symlinked to the root's, so the GTF
    index and the bowtie2 contaminant index are built once (by `ribomine setup`,
    in prep) and merely *read* by all eight nodes. RiboMine notices they are
    current and does not rebuild them -- which is the point: eight
    `bowtie2-build`s into one prefix is corruption, not a slowdown.
    """
    wd = shard_workdir(root, k)
    os.makedirs(wd, exist_ok=True)
    _link(os.path.join(root, "refs"), os.path.join(wd, "refs"))
    os.makedirs(os.path.join(wd, "meta"), exist_ok=True)
    _link(os.path.join(root, "meta", "candidates.tsv"),
          os.path.join(wd, "meta", "candidates.tsv"))
    return wd


def _link(target: str, link: str) -> None:
    if not os.path.exists(target):
        return
    if os.path.islink(link):
        if os.path.realpath(link) == os.path.realpath(target):
            return
        os.unlink(link)
    elif os.path.exists(link):
        return              # a real file/dir already there: leave it alone
    os.symlink(os.path.abspath(target), link)


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-c", "--config", required=True)
    p.add_argument("--shards", type=int, required=True,
                   help="number of nodes the run will use (jobs x nodes_per_job)")
    p.add_argument("--accessions", default=None,
                   help="skip the query and split this list instead")
    p.add_argument("--reshard", action="store_true",
                   help="repack every accession from scratch. Only safe on a workdir "
                        "with no results in it: an accession that changes shard has "
                        "its finished output stranded in a workdir nobody reads, and "
                        "will be downloaded and mapped again.")
    a = p.parse_args(argv)

    cfg = cfgmod.load(a.config)
    setup_logging("INFO", os.path.join(cfg.dir("logs"), "prep.log"))
    root = cfg.workdir
    n = a.shards
    if n < 1:
        print("--shards must be >= 1", file=sys.stderr)
        return 2

    print(f"root workdir : {root}")
    print(f"shards       : {n}")

    # `ribomine setup` must already have run: it builds refs/ (the GTF index and
    # the bowtie2 contaminant index), which every shard workdir then symlinks and
    # only reads. Without it each node discovers them missing and builds its own --
    # eight GTF parses and eight bowtie2-builds, for nothing.
    if not os.path.isdir(os.path.join(root, "refs")):
        print(f"WARNING: {root}/refs does not exist -- run `ribomine setup -c {a.config}` "
              f"first, or every node will build its own copy of the indexes.",
              file=sys.stderr)

    accs, weights = accessions(cfg, a.accessions)
    if not accs:
        print("the query returned nothing -- there is nothing to run", file=sys.stderr)
        return 1

    pins = {} if a.reshard else existing_pins(root, n)
    if pins:
        print(f"pinned       : {len(pins)} accession(s) keep the shard they already have")

    os.makedirs(shard_dir(root), exist_ok=True)
    bins = pack(accs, weights, n, pins)

    for k, accs_k in enumerate(bins):
        with open(shard_list(root, k), "w") as fh:
            fh.write("\n".join(accs_k) + ("\n" if accs_k else ""))
        seed_workdir(root, k, cfg)

    with open(manifest_path(root), "w") as fh:
        json.dump({"n_shards": n, "n_accessions": len(accs),
                   "shards": {f"shard_{k:02d}": len(b) for k, b in enumerate(bins)}},
                  fh, indent=1)
        fh.write("\n")

    print(f"\nwrote {n} shard list(s) and workdir(s) under {shard_dir(root)}")
    for k, b in enumerate(bins):
        print(f"  shard_{k:02d}  {len(b):5d} runs -> {shard_workdir(root, k)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Split a cohort across the nodes of the SLURM run. Called by `slurm/prep.sh`.

    prepare_run.py -c CONFIG --shards N [--accessions LIST] [--reshard]

1. Cohort: `pipeline.accession_list` (or `--accessions`), otherwise
   `<root>/meta/candidates.tsv`. `slurm/master.sh` writes that file on the login
   node; the query is run from here only if the file is missing.
2. Split: one accession list per node (`<root>/shards/shard_NN.txt`), packed by
   `read_count`, largest first, so that the nodes carry similar loads.
3. Workdirs: one `<root>/shards/work_NN` per node, because `ribomine run` writes
   its cohort tables at fixed paths under its workdir. `refs/` and
   `meta/candidates.tsv` are symlinked from the root and only read.

The split is pinned: on a re-run an accession that already has a shard keeps it,
and only new accessions are packed. Moving a finished accession would leave its
results in a workdir that is no longer read, and the run would be processed again.
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
    """Return the cohort and each run's read count (the packing weight).

    The source follows `pipeline.start`, as for `ribomine run`: a config with its
    own `accession_list` is split as given, without querying the archive.
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
        # A bare list has no sizes. If an earlier query left a candidates.tsv in this
        # workdir, take the read counts from it so the packing is still by size.
        w = _weights_from_candidates(cfg, accs)
        if w:
            print(f"            {len(w)}/{len(accs)} weighted by read_count from "
                  f"meta/candidates.tsv")
        return accs, w

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


def _weights_from_candidates(cfg, accs: list[str]) -> dict[str, int]:
    """read_count for the runs in `accs`, from a candidates.tsv left by an earlier
    query in this workdir. Empty if there is none; all runs then weigh the same."""
    tsv = os.path.join(cfg.workdir, "meta", "candidates.tsv")
    if not (os.path.exists(tsv) and os.path.getsize(tsv) > 0):
        return {}
    want = set(accs)
    out: dict[str, int] = {}
    try:
        for r in read_tsv(tsv):
            a = r.get("run_accession")
            if a in want:
                try:
                    out[a] = int(r["read_count"])
                except (KeyError, TypeError, ValueError):
                    pass
    except OSError:
        return {}
    return out


def pack(accs: list[str], weights: dict[str, int], n: int,
         pinned: dict[str, int]) -> list[list[str]]:
    """Pack `accs` into `n` shards, keeping `pinned` (acc -> shard) in place.

    Longest-processing-time-first: the heaviest remaining run goes to the least
    loaded shard. A run without a read count gets the median weight.
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
    """acc -> shard, from the shard lists already on disk (the files the nodes
    read), not from the manifest."""
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
    """Create one node's workdir.

    `refs/` and `meta/candidates.tsv` are symlinks to the root's, so the GTF index
    and the bowtie2 contaminant index are built once (by `ribomine setup`, in prep)
    and only read by the nodes.
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

    # refs/ is built by `ribomine setup` and symlinked into every shard workdir.
    # Without it each node would build its own copy of the indexes.
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

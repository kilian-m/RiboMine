"""Command line interface.

    ribomine init-config config.json        # write a fully-commented default config
    ribomine run -c config.json             # the pipeline, start/end from the config
    ribomine run -c config.json --to qc     # ... stopping after the QC verdict
    ribomine run -c config.json --from accessions --accessions runs.txt
    ribomine query -c config.json           # just the SRA search -> candidates.tsv
    ribomine setup -c config.json           # build the annotation + contaminant indexes
    ribomine benchmark SRR12285169          # time every download route, pick the fastest

Everything the CLI can set also lives in the config; the flags exist so a config
does not have to be edited to move an end point. A flag always wins over the file.
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__, config as cfgmod
from .config import END_POINTS, START_POINTS, ConfigError
from .utils import LOG, setup_logging


def _overrides(a: argparse.Namespace) -> dict:
    """CLI flags -> a config overlay. Only flags the user actually passed appear,
    so an unset flag never overwrites the config file."""
    o: dict = {"project": {}, "pipeline": {}, "process": {}, "plots": {}, "query": {}}
    if a.workdir:
        o["project"]["workdir"] = a.workdir
    if a.jobs:
        o["project"]["jobs"] = a.jobs
    if a.threads:
        o["project"]["threads"] = a.threads
    if getattr(a, "start", None):
        o["pipeline"]["start"] = a.start
    if getattr(a, "end", None):
        o["pipeline"]["end"] = a.end
    if getattr(a, "accessions", None):
        o["pipeline"]["start"] = "accessions"
        o["pipeline"]["accession_list"] = a.accessions
    if getattr(a, "fastq_dir", None):
        o["pipeline"]["start"] = "fastq"
        o["pipeline"]["fastq_dir"] = a.fastq_dir
    if getattr(a, "no_resume", False):
        o["pipeline"]["resume"] = False
    if getattr(a, "umi_dedup", False):
        o["process"]["umi_dedup"] = True
    if getattr(a, "no_plots", False):
        o["plots"]["enabled"] = False
    if getattr(a, "max_runs", None):
        o["query"]["max_runs"] = a.max_runs
    return {k: v for k, v in o.items() if v}


def _load(a: argparse.Namespace):
    cfg = cfgmod.load(a.config, _overrides(a))
    setup_logging(a.log_level, os.path.join(cfg.dir("logs"), "ribomine.log"))
    return cfg


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="ribomine", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"RiboMine {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, *, with_config=True):
        if with_config:
            sp.add_argument("-c", "--config", default=None,
                            help="config JSON (defaults are used for anything it omits)")
        sp.add_argument("--workdir", default=None, help="override project.workdir")
        sp.add_argument("--jobs", type=int, default=None, help="datasets in parallel")
        sp.add_argument("--threads", type=int, default=None, help="threads per dataset")
        sp.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
        return sp

    # -- run ----------------------------------------------------------------
    r = common(sub.add_parser("run", help="run the pipeline"))
    r.add_argument("--from", dest="start", choices=START_POINTS, default=None,
                   help="start point (default: from the config)")
    r.add_argument("--to", dest="end", choices=END_POINTS, default=None,
                   help="end point (default: from the config)")
    r.add_argument("--accessions", default=None,
                   help="file of run accessions, one per line (implies --from accessions)")
    r.add_argument("--fastq-dir", default=None,
                   help="directory of local FASTQs (implies --from fastq)")
    r.add_argument("--max-runs", type=int, default=None, help="cap the number of datasets")
    r.add_argument("--umi-dedup", action="store_true",
                   help="deduplicate on the UMI after mapping (off by default)")
    r.add_argument("--no-plots", action="store_true")
    r.add_argument("--no-resume", action="store_true",
                   help="recompute everything instead of reusing existing outputs")

    # -- query --------------------------------------------------------------
    q = common(sub.add_parser("query", help="search the SRA/ENA for ribo-seq runs"))
    q.add_argument("--max-runs", type=int, default=None)

    # -- setup --------------------------------------------------------------
    s = common(sub.add_parser("setup", help="build the annotation + contaminant indexes"))
    s.add_argument("--contaminant-fasta", default=None,
                   help="build the bowtie2 contaminant index from this FASTA")

    # -- init-config --------------------------------------------------------
    i = sub.add_parser("init-config", help="write a config file with every default")
    i.add_argument("path", nargs="?", default="config.json")

    # -- benchmark ----------------------------------------------------------
    b = common(sub.add_parser("benchmark", help="time each SRA download route"))
    b.add_argument("accession", help="a run accession to benchmark with (a small one is fine)")

    a = p.parse_args(argv)

    try:
        return _dispatch(a)
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


def _dispatch(a: argparse.Namespace) -> int:
    if a.cmd == "init-config":
        cfgmod.write_example(a.path)
        print(f"wrote {a.path} — every key is the pipeline's own default; edit and pass with -c")
        return 0

    cfg = _load(a)

    if a.cmd == "run":
        from . import pipeline
        pipeline.run(cfg)
        return 0

    if a.cmd == "query":
        from .sra import query
        # the query neither reads the pipeline's start point nor touches the genome
        cfg.validate(need_reference=False, need_inputs=False)
        path = query.run_query(cfg)
        print(f"\ncandidates -> {path}")
        return 0

    if a.cmd == "setup":
        from .qc import annotation, contaminants
        print(f"annotation index:  {annotation.ensure_index(cfg)}")
        fa = a.contaminant_fasta or contaminants.resolve_fasta(cfg)
        bundled = not (a.contaminant_fasta or cfg.ref("contaminant_fasta"))
        print(f"contaminant FASTA: {fa}" + ("  (bundled with RiboMine)" if bundled else ""))
        print(f"contaminant index: {contaminants.build_index(fa, cfg)}")
        return 0

    if a.cmd == "benchmark":
        from .sra import download
        rows = download.benchmark(a.accession, cfg)
        w = max(len(r["route"]) for r in rows)
        print(f"\n{'route':<{w}}  {'MB/s':>7}  {'seconds':>8}  note")
        for r in sorted(rows, key=lambda x: -(x.get("mb_per_s") or 0)):
            mbps = f"{r['mb_per_s']:.1f}" if r.get("mb_per_s") else "—"
            secs = f"{r['seconds']:.1f}" if r.get("seconds") else "—"
            print(f"{r['route']:<{w}}  {mbps:>7}  {secs:>8}  {r.get('note', '')}")
        return 0

    raise AssertionError(f"unhandled command {a.cmd}")


if __name__ == "__main__":
    sys.exit(main())

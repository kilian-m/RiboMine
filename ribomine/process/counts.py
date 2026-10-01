"""The read-count matrix: genes x runs, from STAR's `--quantMode GeneCounts`.

STAR counts reads into genes while it aligns them. Each run's
`star_final/ReadsPerGene.out.tab` is one column; this module reads and joins them.

* Counting follows htseq-count's `union` model: a read counts for a gene if it
  overlaps that gene's exons and no other gene's. Multimappers, reads in two
  genes (`ambiguous`) and reads in none (`no_feature`) are not counted. Counts
  cover all exons of a gene, so they include UTR and ncRNA signal.
* The sense column is used, because a footprint maps to the transcript's own
  strand. A library whose antisense counts rival its sense counts (reversed or
  unstranded) gets a warning: its sense column is an undercount.
* Counts are taken during alignment, before the pile-up filter and optional UMI
  deduplication, so they can exceed the reads in the final BAM.
  `mapping_summary.tsv` reports both totals.
"""
from __future__ import annotations

import os

from ..utils import LOG, nonempty, write_tsv
from .star import GENE_COUNTS, GENE_INFO

# ReadsPerGene.out.tab columns: gene_id, unstranded, sense, antisense. `sense`
# is htseq-count -s yes (read on the RNA strand), as for single-end ribo-seq.
COL_UNSTRANDED, COL_SENSE, COL_ANTISENSE = 1, 2, 3
# the four summary rows STAR puts above the genes
SUMMARY_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
# below this sense:antisense ratio the library is not considered sense-stranded
MIN_SENSE_RATIO = 2.0


def path(star_final_dir: str) -> str:
    """Path of a run's `ReadsPerGene.out.tab`."""
    return os.path.join(star_final_dir, GENE_COUNTS)


def read_counts(tab: str, *, label: str = "") -> tuple[dict[str, int], dict]:
    """Parse one `ReadsPerGene.out.tab` -> ({gene_id: sense count}, stats).

    The stats include the reads that fell in no gene or in two (`n_no_feature`,
    `n_ambiguous`), which form the denominator of `frac_in_genes`.
    """
    counts: dict[str, int] = {}
    summary: dict[str, int] = {}
    n_sense = n_anti = 0
    with open(tab) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= COL_ANTISENSE:
                continue
            gene = f[0]
            sense = int(f[COL_SENSE])
            if gene in SUMMARY_ROWS:
                summary[gene] = sense
                continue
            counts[gene] = sense
            n_sense += sense
            n_anti += int(f[COL_ANTISENSE])

    n_in_genes = sum(counts.values())
    # STAR's N_multimapping row is not reported: with `mapping.multimap_nmax: 1`
    # STAR drops multimappers before counting, so the row reads 0. The multimapping
    # rate comes from Log.final.out (the `mapping` block).
    stats = {
        "n_in_genes": n_in_genes,
        "n_genes_detected": sum(1 for v in counts.values() if v > 0),
        "n_no_feature": summary.get("N_noFeature"),
        "n_ambiguous": summary.get("N_ambiguous"),
        "n_antisense": n_anti,
        # A sense-stranded library is at 10-100x; near 1 the reads are not on the
        # sense strand and the counts are an undercount (warned about below).
        "sense_over_antisense": round(n_sense / n_anti, 1) if n_anti else None,
    }
    counted = n_in_genes + (stats["n_no_feature"] or 0) + (stats["n_ambiguous"] or 0)
    if counted:
        stats["frac_in_genes"] = round(n_in_genes / counted, 4)
    ratio = stats["sense_over_antisense"]
    if ratio is not None and ratio < MIN_SENSE_RATIO and n_sense + n_anti > 1000:
        LOG.warning(
            "[%s] gene counts are only %.1fx sense over antisense -- this library does "
            "not look sense-stranded, and the counts are of its sense strand only",
            label or os.path.basename(tab), ratio)
    return counts, stats


def gene_table(star_index: str) -> list[tuple[str, str]]:
    """[(gene_id, gene_name)] in the index's own order: the rows of the matrix.

    STAR wrote `geneInfo.tab` from the GTF it counts against, so the rows match
    the count tables. Line 1 is the gene count, then id / name / biotype.
    """
    rows: list[tuple[str, str]] = []
    with open(os.path.join(star_index, GENE_INFO)) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f and f[0]:
                rows.append((f[0], f[1] if len(f) > 1 else ""))
    return rows


def matrix(out_tsv: str, star_index: str, columns: list[tuple[str, str]]) -> str:
    """Write the genes x runs matrix. `columns` is [(accession, ReadsPerGene path)].

    Every annotated gene gets a row, including all-zero ones, so that matrices
    from different batches share one row set. Runs without counts are left out
    as columns rather than written as blanks, which would not parse as numbers;
    `mapping_summary.tsv` still has a row for them.
    """
    genes = gene_table(star_index)
    per_run = {}
    for acc, tab in columns:
        if not nonempty(tab):
            LOG.debug("no gene counts for %s (%s)", acc, tab)
            continue
        per_run[acc] = read_counts(tab, label=acc)[0]
    if not per_run:
        raise RuntimeError("no run produced gene counts; nothing to build a matrix from")

    accs = sorted(per_run)
    rows = []
    for gid, name in genes:
        row = {"gene_id": gid, "gene_name": name}
        for acc in accs:
            row[acc] = per_run[acc].get(gid, 0)
        rows.append(row)

    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)) or ".", exist_ok=True)
    write_tsv(out_tsv, rows, columns=["gene_id", "gene_name", *accs])
    LOG.info("gene_counts.tsv: %d genes x %d run(s) -> %s", len(rows), len(accs), out_tsv)
    return out_tsv

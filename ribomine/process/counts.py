"""The read-count matrix: genes x runs, straight out of STAR.

STAR counts reads into genes *while it aligns them* (`--quantMode GeneCounts`),
so the matrix costs nothing beyond the alignment we already run -- no second pass
over the BAM, no htseq/featureCounts, no gene model of our own to keep in step
with the one STAR aligned against. Each run's `star_final/ReadsPerGene.out.tab`
is one column; this module reads them and joins them.

**What a count is.** STAR's rule (htseq-count's `union` model): a read is counted
into a gene if it overlaps that gene's exons and no other gene's. A read
overlapping two genes is `ambiguous` and counted into neither; a read in no gene
at all (intron, intergenic) is `no_feature`; a multimapper is not counted. The
count is therefore over the WHOLE GENE -- every exon of every biotype -- not over
the CDS. For a ribo-seq library the CDS is where the footprints are, so the two
are close but not the same number, and UTR / ncRNA signal is inside these counts.

**The strand.** The table carries three count columns -- unstranded, sense and
antisense -- and RiboMine takes the SENSE one: a ribosome footprint is a piece of
the mRNA, so it maps to the transcript's own strand. The antisense column is read
too, and a library whose antisense counts rival its sense counts gets a warning:
that is a reversed or unstranded library, and its sense column is then an
undercount, not a measurement.

**What these counts are NOT.** They are taken during the alignment, which is
*before* the pile-up filter and before optional UMI deduplication. So a gene
sitting under an adapter-dimer pile is counted with that pile included, and with
`process.umi_dedup` on the counts are of the duplicated reads while the BAM is of
the deduplicated ones. `mapping_summary.tsv` reports both read totals side by
side, so the gap is visible rather than implied.
"""
from __future__ import annotations

import os

from ..utils import LOG, nonempty, write_tsv
from .star import GENE_COUNTS, GENE_INFO

# ReadsPerGene.out.tab: gene_id, unstranded, sense, antisense. `sense` is
# htseq-count -s yes ("1st read strand aligned with RNA"), which is what a
# single-end ribo-seq read is.
COL_UNSTRANDED, COL_SENSE, COL_ANTISENSE = 1, 2, 3
# the four bookkeeping rows STAR puts above the genes
SUMMARY_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
# below this sense:antisense ratio the library does not read as sense-stranded
MIN_SENSE_RATIO = 2.0


def path(star_final_dir: str) -> str:
    """Where STAR left this run's counts."""
    return os.path.join(star_final_dir, GENE_COUNTS)


def read_counts(tab: str, *, label: str = "") -> tuple[dict[str, int], dict]:
    """Parse one `ReadsPerGene.out.tab` -> ({gene_id: sense count}, summary).

    The summary is what the counts do NOT contain: reads that fell in no gene, and
    reads that straddled two. Those are the honest denominator for "what fraction of
    this library landed in a gene", so they are carried into the report rather than
    dropped on the floor.
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
    # STAR's N_multimapping row is NOT read. Under the default `mapping.multimap_nmax:
    # 1` STAR drops multimappers before it ever counts, so that row says 0 while the
    # library may be 90% multimapping -- a number that is worse than no number. The
    # real one comes from Log.final.out and is already in the `mapping` block.
    stats = {
        "n_in_genes": n_in_genes,
        "n_genes_detected": sum(1 for v in counts.values() if v > 0),
        "n_no_feature": summary.get("N_noFeature"),
        "n_ambiguous": summary.get("N_ambiguous"),
        "n_antisense": n_anti,
        # A sense-stranded library sits at 10-100x here. Near 1 means the strand
        # column we are reading is not the one the reads are on, and the matrix
        # would be a fraction of the real signal -- say so rather than ship it.
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
    """[(gene_id, gene_name)] in the index's own order -- the matrix's rows.

    STAR wrote `geneInfo.tab` from the same GTF it counted against, so this row set
    matches the count tables exactly; there is no join to get wrong. Line 1 is the
    gene count, then id / name / biotype.
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

    Every gene in the annotation gets a row, including the all-zero ones: a matrix
    whose row set depends on which runs happened to be in it cannot be compared with
    the next one, and every downstream tool would rather filter zeros itself than
    guess why a gene is missing.

    Runs with no counts are LEFT OUT as columns rather than filled with blanks --
    a blank is not a zero, and the difference matters to everything that reads this
    file as a numeric table. `mapping_summary.tsv` still has a row for them.
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

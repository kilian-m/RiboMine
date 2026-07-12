"""Reference data shipped with RiboMine.

Only things small enough to live in the repo and stable enough to pin a version of
belong here. The genome and the GTF do not; the contaminant reference does.

`human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt_rDNA.fa` (672 kB, 3,858
sequences: 1,900 snRNA, 943 snoRNA, 559 rRNA, 454 tRNA, 24 Mt, plus the two rDNA
entries below) is the sequence set a Ribo-seq library needs stripped before mapping.
It is *not* an adapter list -- RiboMine hard-codes no adapter anywhere; adapters are
an output of the architecture stage, and adapter/primer dimers are removed after
mapping, by position (`qc/pileups.py`).

Why the two rDNA entries matter
-------------------------------
The 3,856 mature species cover 18S / 5.8S / 28S / 5S, but a Ribo-seq library also
carries fragments of the **unprocessed precursor** -- the internal and external
transcribed spacers (ITS1/2, 5'ETS/3'ETS) that are excised during rRNA maturation and
so appear in *no* mature sequence. Two entries cover them:

    NR_046235.3   45S pre-ribosomal RNA (RNA45SN5)          13,357 bp
    U13369.1      human ribosomal DNA complete repeating unit  42,999 bp

Measured contribution, on the samples RiboMine was developed against:

    SRR618773   (rRNA-heavy)  contaminant removal  73.1% -> 74.8%   (+1.7 pp)
    SRR12285169 (clean)                            32.8% -> 32.8%
    SRR1039508  (RNA-seq)                           1.1% ->  1.4%

The gain is small but it is real, and every read it catches is a read that would
otherwise have gone to the genome aligner and, because rDNA is a high-copy repeat,
come back as a multimapper. No verdict and no architecture call changed.

Organism
--------
This is human. For another organism, point `reference.contaminant_fasta` at that
organism's rRNA / tRNA / sn(o)RNA / Mt sequences -- and include its pre-rRNA /
rDNA repeat for the same reason. `ribomine setup` builds the bowtie2 index from
whatever it is given.
"""
from __future__ import annotations

import os

HUMAN_CONTAMINANTS = "human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt_rDNA.fa"


def path(name: str) -> str:
    """Absolute path to a packaged data file.

    `importlib.resources` is the correct way to do this: it keeps working when
    RiboMine is installed as a wheel or a zip, where __file__-relative paths would
    not. The `as_file` context manager is not needed because we ship a real
    directory, never a zipped resource.
    """
    try:
        from importlib.resources import files

        p = files(__package__) / name
        if p.is_file():
            return str(p)
    except (ImportError, ModuleNotFoundError, TypeError):
        pass
    p2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    if os.path.isfile(p2):
        return p2
    raise FileNotFoundError(f"packaged data file not found: {name}")


def human_contaminants() -> str:
    """The bundled human rRNA / tRNA / snRNA / snoRNA / Mt contaminant FASTA."""
    return path(HUMAN_CONTAMINANTS)

"""Reference data shipped with RiboMine.

`human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt_rDNA.fa` holds the 3,858
sequences removed before mapping: 1,900 snRNA, 943 snoRNA, 559 rRNA and 454
tRNA (24 of them mitochondrial), plus the 45S pre-rRNA (NR_046235.3) and the
complete rDNA repeat (U13369.1). The last two cover the transcribed spacers
(ITS1/2, 5'/3' ETS), which are absent from the mature rRNA sequences and would
otherwise multimap on the genome. The file contains no adapters; those are
detected per run by the architecture stage.

The set is human. For another organism, point `reference.contaminant_fasta` at
its rRNA / tRNA / sn(o)RNA / Mt sequences, including the pre-rRNA / rDNA repeat.
"""
from __future__ import annotations

import os

HUMAN_CONTAMINANTS = "human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt_rDNA.fa"


def path(name: str) -> str:
    """Absolute path to a packaged data file (via `importlib.resources`, with a
    `__file__`-relative fallback)."""
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

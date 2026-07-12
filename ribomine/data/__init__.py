"""Reference data shipped with RiboMine.

Only things small enough to live in the repo and stable enough to pin a version of
belong here. The genome and the GTF do not; the contaminant reference does.

`human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt.fa` (629 kB, 3,856 sequences:
1,900 snRNA, 943 snoRNA, 559 rRNA, 454 tRNA, 24 Mt) is the sequence set that a
Ribo-seq library needs stripped before mapping. It is *not* an adapter list --
RiboMine hard-codes no adapter anywhere; adapters are an output of the architecture
stage, and adapter/primer dimers are removed after mapping, by position
(`qc/pileups.py`).

It is human. For another organism, point `reference.contaminant_fasta` at that
organism's rRNA/tRNA/sn(o)RNA/Mt sequences; `ribomine setup` builds the bowtie2
index from whatever it is given.

One measured caveat, if you are chasing the last percent. This reference contains
the *mature* rRNA species but not the 45S pre-rRNA (NR_046235.3) or the complete
rDNA repeating unit (U13369.1). Adding those two sequences catches a little more:
on an rRNA-heavy library (SRR618773) contaminant removal goes 73.1% -> 74.8%, and on
a clean one (SRR12285169) 32.8% -> 32.8% -- i.e. the extra 1.7 pp is pre-rRNA spacer
(ITS/ETS) fragments that the mature species do not cover. It changed no verdict and
no architecture call on the samples tested. Append them if you want the extra catch;
nothing in RiboMine depends on their absence.
"""
from __future__ import annotations

import os

HUMAN_CONTAMINANTS = "human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt.fa"


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

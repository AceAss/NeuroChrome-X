"""Minimal reader for MEME-format motif files (for example JASPAR downloads)."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np

from neurochromx.motifs.scan import Motif, motif_from_probabilities


def read_meme(path: str | Path) -> list[Motif]:
    """Parse ``MOTIF <id> [<name>]`` blocks with a ``letter-probability matrix``.

    The motif name is ``<name>`` when present, else ``<id>``. The ID is kept as
    ``Motif.name`` suffix ``" (id)"`` only when both exist, so that selection by TF name
    works (see :func:`select_motifs`).
    """
    lines = Path(path).read_text().splitlines()
    motifs: list[Motif] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MOTIF"):
            parts = line.split()
            ident = parts[1]
            name = parts[2] if len(parts) > 2 else ident
            label = name if name == ident else f"{name} ({ident})"
            while i < len(lines) and not lines[i].startswith(
                "letter-probability matrix"
            ):
                i += 1
            width = int(lines[i].split("w=")[1].split()[0])
            rows = [[float(x) for x in lines[i + 1 + k].split()] for k in range(width)]
            motifs.append(motif_from_probabilities(label, np.array(rows)))
            i += width
        i += 1
    return motifs


def select_motifs(motifs: Iterable[Motif], tf_names: Iterable[str]) -> list[Motif]:
    """Motifs whose name (before any ``" (id)"``) equals a wanted TF name (any case)."""
    wanted = {n.upper() for n in tf_names}
    return [m for m in motifs if m.name.split(" (")[0].upper() in wanted]

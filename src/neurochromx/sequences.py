"""Fetch genomic sequences for regions (needs samtools on PATH)."""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Sequence
from pathlib import Path

from neurochromx.intervals import Region


def faidx_sequences(fasta: str | Path, regions: Sequence[Region]) -> list[str]:
    """Upper-case sequence of each region (0-based, half-open), in the given order.

    Uses ``samtools faidx -r``; the FASTA index is created next to the FASTA if it is
    missing, so the directory must be writable the first time.
    """
    if not regions:
        return []
    with tempfile.TemporaryDirectory() as tmp:
        listing = Path(tmp) / "regions.txt"
        listing.write_text(
            "".join(f"{r.chrom}:{r.start + 1}-{r.end}\n" for r in regions)
        )
        done = subprocess.run(
            ["samtools", "faidx", "-r", str(listing), str(fasta)],
            capture_output=True,
            text=True,
            check=True,
        )
    sequences: list[list[str]] = []
    for line in done.stdout.splitlines():
        if line.startswith(">"):
            sequences.append([])
        elif sequences:
            sequences[-1].append(line.strip())
    out = ["".join(parts).upper() for parts in sequences]
    if len(out) != len(regions):
        raise RuntimeError(f"expected {len(regions)} sequences, got {len(out)}")
    return out

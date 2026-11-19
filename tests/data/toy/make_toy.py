#!/usr/bin/env python3
"""Deterministic synthetic toy dataset for NeuroChrom-X CI (TC-004, built by Claude).

Everything is synthetic, so CI needs no download, and the truth is known: planted
differential accessibility, planted TF motifs, and gene expression linked to the
changed regions. The PRD section 11 gates (positive controls, label-shuffle null)
can therefore be tested exactly.

Usage::

    python tests/data/toy/make_toy.py --out tests/data/toy/generated [--seed 20261105]
                                      [--depth 70000]

Outputs (under --out): toy.fa, toy.chrom.sizes, toy.gtf, toy_blacklist.bed,
toy_ccres.bed, fastq/<sample_id>_R{1,2}.fastq.gz (6 ATAC samples),
rna_counts.tsv (6 RNA samples), truth_regions.tsv, truth_genes.tsv, MANIFEST.json.

Determinism: all randomness comes from numpy ``SeedSequence(seed)`` children, and gzip
files carry no timestamp or name. MANIFEST.json holds the SHA-256 of the uncompressed
content of every file, so two runs with the same seed and numpy version agree.
(numpy does not promise identical random streams across versions: pin it.)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

GENERATOR_VERSION = "1.0.0"
DEFAULT_SEED = 20261105
DEFAULT_DEPTH = 70_000

CHR19, CHRM = "chr19", "chrM"
CHR19_LEN, CHRM_LEN = 3_000_000, 16_299
READ_LEN = 50
IEG_GENES = ("Fos", "Arc", "Npas4", "Egr1", "Junb", "Nr4a1")
MOTIFS = {"AP1": "TGAGTCA", "CREB": "TGACGTCA", "MEF2": "CTAAAAATAG"}
N_GENES, N_DISTAL = 60, 340
N_GAINED_DISTAL, N_LOST = 42, 12
CCRE_FRACTION = 0.70
LINK_WINDOW = 100_000
# Share of read pairs by origin: regions, blacklist pile-ups, uniform background, chrM.
MIX = {"region": 0.66, "blacklist": 0.04, "background": 0.18, "chrM": 0.12}
DUP_FRACTION = 0.10
PLANT_GAINED = {"AP1": 0.8, "CREB": 0.4, "MEF2": 0.3}
PLANT_UNCHANGED = 0.08
CONDITIONS = (("control", "ctrl"), ("stimulated", "stim"))
REPLICATES = (1, 2, 3)

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_COMP = np.zeros(256, dtype=np.uint8)
for _a, _b in zip(b"ACGT", b"TGCA", strict=True):
    _COMP[_a] = _b
_IDX = np.zeros(256, dtype=np.int64)
for _i, _a in enumerate(b"ACGT"):
    _IDX[_a] = _i


def sample_ids(assay: str) -> list[str]:
    """Sample IDs in the order used everywhere: control r1-r3, then stimulated r1-r3."""
    return [f"toy_{assay}_{tag}_r{rep}" for _, tag in CONDITIONS for rep in REPLICATES]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# --------------------------------------------------------------------------- #
# Genome, genes, regions
# --------------------------------------------------------------------------- #
def make_genome(rng: np.random.Generator) -> np.ndarray:
    """Random chr19 (GC drifts by 50 kb window) then chrM, as one uint8 array."""
    window = 50_000
    gc = np.clip(rng.normal(0.42, 0.04, CHR19_LEN // window), 0.30, 0.55)
    parts = []
    for value in gc:
        probs = [(1 - value) / 2, value / 2, value / 2, (1 - value) / 2]
        parts.append(_BASES[rng.choice(4, size=window, p=probs)])
    chrm = _BASES[rng.choice(4, size=CHRM_LEN, p=[0.29, 0.21, 0.21, 0.29])]
    return np.concatenate([*parts, chrm])


def make_genes(rng: np.random.Generator) -> list[dict[str, Any]]:
    slot = CHR19_LEN // N_GENES
    ieg_slots = set(rng.choice(N_GENES, size=len(IEG_GENES), replace=False).tolist())
    ieg_names = iter(rng.permutation(list(IEG_GENES)).tolist())
    genes: list[dict[str, Any]] = []
    toy_no = 0
    for i in range(N_GENES):
        start = i * slot + 10_000 + int(rng.integers(0, 8_000))
        exons, pos = [], start
        for _ in range(int(rng.integers(3, 7))):
            length = int(rng.integers(100, 401))
            exons.append((pos, pos + length))
            pos += length + int(rng.integers(200, 801))
        end = exons[-1][1]
        strand = "+" if rng.random() < 0.5 else "-"
        if i in ieg_slots:
            name = next(ieg_names)
        else:
            toy_no += 1
            name = f"ToyGene{toy_no:02d}"
        genes.append(
            {
                "gene_id": f"TOYG{i + 1:04d}",
                "gene_name": name,
                "start": start,
                "end": end,
                "strand": strand,
                "tss": start if strand == "+" else end - 1,
                "exons": exons,
            }
        )
    return genes


def promoter_interval(gene: dict[str, Any]) -> tuple[int, int]:
    tss = gene["tss"]
    return (tss - 300, tss + 200) if gene["strand"] == "+" else (tss - 200, tss + 300)


def _free(start: int, end: int, taken: list[tuple[int, int]], gap: int) -> bool:
    return not any(start < e + gap and s < end + gap for s, e in taken)


def _far_from_tss(start: int, end: int, tss: np.ndarray, distance: int) -> bool:
    gap = np.maximum(0, np.maximum(start - tss, tss - end))
    return bool(gap.min() >= distance)


def place_regions(
    rng: np.random.Generator,
    genes: list[dict[str, Any]],
    taken: list[tuple[int, int]],
    count: int,
    widths: tuple[int, int],
    tss_distance: int,
    gap: int,
) -> list[tuple[int, int]]:
    tss = np.array([g["tss"] for g in genes])
    placed: list[tuple[int, int]] = []
    for _ in range(500_000):
        if len(placed) == count:
            return placed
        width = int(rng.integers(widths[0], widths[1] + 1))
        start = int(rng.integers(5_000, CHR19_LEN - 5_000 - width))
        end = start + width
        if _far_from_tss(start, end, tss, tss_distance) and _free(
            start, end, taken, gap
        ):
            placed.append((start, end))
            taken.append((start, end))
    raise RuntimeError("could not place all regions; lower the counts")


def plant_motif(
    genome: np.ndarray,
    rng: np.random.Generator,
    region: tuple[int, int],
    motif: str,
    used: list[tuple[int, int]],
) -> int:
    """Insert 1-2 copies of ``motif`` (random strand) into ``region``."""
    copies = int(rng.integers(1, 3))
    placed = 0
    for _ in range(50):
        if placed == copies:
            break
        pos = int(rng.integers(region[0] + 20, region[1] - 20 - len(motif)))
        if not _free(pos, pos + len(motif), used, 2):
            continue
        text = motif if rng.random() < 0.5 else revcomp(motif)
        genome[pos : pos + len(motif)] = np.frombuffer(text.encode(), dtype=np.uint8)
        used.append((pos, pos + len(motif)))
        placed += 1
    return placed


def build_truth(
    rng: np.random.Generator, genome: np.ndarray, genes: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[tuple[int, int]]]:
    """Regions with truth (status, log2FC, motifs) plus the blacklist intervals."""
    promoters = [promoter_interval(g) for g in genes]
    taken = list(promoters)
    distal = place_regions(rng, genes, taken, N_DISTAL, (300, 600), 3_000, 1_000)
    blacklist = place_regions(rng, genes, taken, 2, (5_000, 5_000), 3_000, 3_000)

    regions: list[dict[str, Any]] = []
    for gene, (start, end) in zip(genes, promoters, strict=True):
        regions.append(
            {
                "class": "promoter",
                "start": start,
                "end": end,
                "gene_id": gene["gene_id"],
            }
        )
    for start, end in distal:
        regions.append({"class": "distal", "start": start, "end": end, "gene_id": None})
    regions.sort(key=lambda r: r["start"])
    for number, region in enumerate(regions, start=1):
        region["region_id"] = f"R{number:04d}"

    ieg_ids = {g["gene_id"] for g in genes if g["gene_name"] in IEG_GENES}
    distal_idx = [i for i, r in enumerate(regions) if r["class"] == "distal"]
    chosen = rng.permutation(distal_idx)
    gained_distal = set(chosen[:N_GAINED_DISTAL].tolist())
    lost_distal = set(chosen[N_GAINED_DISTAL : N_GAINED_DISTAL + N_LOST].tolist())
    for i, region in enumerate(regions):
        is_ieg_promoter = region["class"] == "promoter" and region["gene_id"] in ieg_ids
        if is_ieg_promoter or i in gained_distal:
            region["status"] = "gained"
            region["log2fc"] = float(rng.uniform(1.2, 3.0))
        elif i in lost_distal:
            region["status"] = "lost"
            region["log2fc"] = float(rng.uniform(-2.5, -1.2))
        else:
            region["status"] = "unchanged"
            region["log2fc"] = 0.0
        region["motifs"] = []
        region["base"] = float(rng.lognormal(0.0, 0.5))

    gained = [i for i, r in enumerate(regions) if r["status"] == "gained"]
    unchanged = [i for i, r in enumerate(regions) if r["status"] == "unchanged"]
    used: dict[int, list[tuple[int, int]]] = {i: [] for i in range(len(regions))}
    for name, share in PLANT_GAINED.items():
        for i in rng.choice(gained, size=round(share * len(gained)), replace=False):
            region = regions[int(i)]
            if plant_motif(
                genome,
                rng,
                (region["start"], region["end"]),
                MOTIFS[name],
                used[int(i)],
            ):
                region["motifs"].append(name)
    for name in MOTIFS:
        for i in rng.choice(
            unchanged, size=round(PLANT_UNCHANGED * len(unchanged)), replace=False
        ):
            region = regions[int(i)]
            if plant_motif(
                genome,
                rng,
                (region["start"], region["end"]),
                MOTIFS[name],
                used[int(i)],
            ):
                region["motifs"].append(name)

    tss = np.array([g["tss"] for g in genes])
    for region in regions:
        centre = (region["start"] + region["end"]) // 2
        nearest = int(np.argmin(np.abs(tss - centre)))
        within = abs(int(tss[nearest]) - centre) <= LINK_WINDOW
        if region["gene_id"] is None:
            region["gene_id"] = genes[nearest]["gene_id"] if within else "NA"
    return regions, blacklist


# --------------------------------------------------------------------------- #
# Reads
# --------------------------------------------------------------------------- #
def fragment_lengths(rng: np.random.Generator, n: int) -> np.ndarray:
    pick = rng.random(n)
    nfr = np.clip(rng.normal(80, 20, n), 50, 150)
    mono = np.clip(rng.normal(200, 25, n), 150, 320)
    di = np.clip(rng.normal(380, 40, n), 320, 600)
    return np.where(pick < 0.55, nfr, np.where(pick < 0.90, mono, di)).astype(np.int64)


def simulate_fragments(
    rng: np.random.Generator,
    regions: list[dict[str, Any]],
    weights: np.ndarray,
    blacklist: list[tuple[int, int]],
    depth: int,
) -> dict[str, np.ndarray]:
    """Unique fragments for one sample, then exact duplicates added up to ``depth``."""
    n_dup = round(DUP_FRACTION * depth)
    n_unique = depth - n_dup
    counts = {k: round(v * n_unique) for k, v in MIX.items()}
    counts["region"] += n_unique - sum(counts.values())

    starts = np.array([r["start"] for r in regions])
    widths = np.array([r["end"] - r["start"] for r in regions])
    pick = rng.choice(len(regions), size=counts["region"], p=weights / weights.sum())
    centres = [starts[pick] + rng.integers(0, widths[pick])]
    chrom = [np.zeros(counts["region"], dtype=np.int64)]
    which = rng.integers(0, len(blacklist), size=counts["blacklist"])
    bl_start = np.array([b[0] for b in blacklist])
    bl_width = np.array([b[1] - b[0] for b in blacklist])
    centres.append(bl_start[which] + rng.integers(0, bl_width[which]))
    chrom.append(np.zeros(counts["blacklist"], dtype=np.int64))
    centres.append(rng.integers(0, CHR19_LEN, size=counts["background"]))
    chrom.append(np.zeros(counts["background"], dtype=np.int64))
    centres.append(rng.integers(0, CHRM_LEN, size=counts["chrM"]))
    chrom.append(np.ones(counts["chrM"], dtype=np.int64))
    origin = np.concatenate(
        [
            np.full(counts[k], i)
            for i, k in enumerate(("region", "blacklist", "background", "chrM"))
        ]
    )

    centre = np.concatenate(centres)
    on_m = np.concatenate(chrom)
    length = fragment_lengths(rng, n_unique)
    chrom_len = np.where(on_m == 1, CHRM_LEN, CHR19_LEN)
    length = np.minimum(length, chrom_len)
    start = np.clip(centre - length // 2, 0, chrom_len - length)
    flip = rng.random(n_unique) < 0.5

    dup = rng.integers(0, n_unique, size=n_dup)
    take = np.concatenate([np.arange(n_unique), dup])
    order = rng.permutation(len(take))
    take = take[order]
    return {
        "chrom_m": on_m[take],
        "start": start[take],
        "length": length[take],
        "flip": flip[take],
        "origin": origin[take],
        "n_unique": np.array(n_unique),
    }


def extract_reads(
    rng: np.random.Generator, genome: np.ndarray, frag: dict[str, np.ndarray]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    offset = np.where(frag["chrom_m"] == 1, CHR19_LEN, 0)
    pos = offset + frag["start"]
    cols = np.arange(READ_LEN)
    forward = genome[pos[:, None] + cols[None, :]]
    tail = (pos + frag["length"] - READ_LEN)[:, None] + cols[None, :]
    reverse = _COMP[genome[tail]][:, ::-1]
    flip = frag["flip"][:, None]
    mates = []
    for read in (np.where(flip, reverse, forward), np.where(flip, forward, reverse)):
        read = read.copy()
        error = rng.random(read.shape) < 0.01
        shift = rng.integers(1, 4, size=read.shape)
        read[error] = _BASES[(_IDX[read[error]] + shift[error]) % 4]
        qual = np.where(error, ord("5"), ord("I")).astype(np.uint8)
        mates.append((read, qual))
    return mates[0][0], mates[0][1], mates[1][0], mates[1][1]


def fastq_text(sample: str, mate: int, seq: np.ndarray, qual: np.ndarray) -> bytes:
    prefix = sample.encode()
    out = [
        b"@%b:%d/%d\n%b\n+\n%b\n"
        % (prefix, i + 1, mate, seq[i].tobytes(), qual[i].tobytes())
        for i in range(len(seq))
    ]
    return b"".join(out)


def write_gz(path: Path, data: bytes) -> str:
    with path.open("wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0
    ) as handle:
        handle.write(data)
    return hashlib.sha256(data).hexdigest()


def write_text(path: Path, text: str) -> str:
    data = text.encode()
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------- #
# RNA, annotation tables
# --------------------------------------------------------------------------- #
def gene_truth(
    rng: np.random.Generator, genes: list[dict[str, Any]], regions: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    gained_link = {r["gene_id"] for r in regions if r["status"] == "gained"}
    lost_link = {r["gene_id"] for r in regions if r["status"] == "lost"}
    out = []
    for gene in genes:
        gid = gene["gene_id"]
        if gene["gene_name"] in IEG_GENES:
            lfc = float(rng.uniform(1.5, 3.0))
        elif gid in gained_link and rng.random() < 0.5:
            lfc = float(rng.uniform(0.5, 1.5))
        elif gid in lost_link and rng.random() < 0.5:
            lfc = float(rng.uniform(-1.5, -0.5))
        else:
            lfc = 0.0
        status = "gained" if lfc > 0 else "lost" if lfc < 0 else "unchanged"
        base = float(np.exp(rng.normal(np.log(300), 1.0)))
        out.append({**gene, "log2fc": lfc, "status": status, "base_mean": base})
    return out


def rna_counts(rng: np.random.Generator, truth: list[dict[str, Any]]) -> np.ndarray:
    """Negative-binomial counts (dispersion 0.1), columns as ``sample_ids("rna")``."""
    n_samples = len(CONDITIONS) * len(REPLICATES)
    size_factor = rng.lognormal(0.0, 0.1, n_samples)
    counts = np.zeros((len(truth), n_samples), dtype=np.int64)
    for j in range(n_samples):
        stimulated = j >= len(REPLICATES)
        for i, gene in enumerate(truth):
            mean = (
                gene["base_mean"]
                * size_factor[j]
                * (2 ** gene["log2fc"] if stimulated else 1.0)
            )
            r = 1.0 / 0.1
            counts[i, j] = rng.negative_binomial(r, r / (r + mean))
    return counts


def gtf_text(genes: list[dict[str, Any]]) -> str:
    lines = []
    for g in genes:
        attrs = (
            f'gene_id "{g["gene_id"]}"; gene_name "{g["gene_name"]}"; '
            'gene_type "protein_coding";'
        )
        tx = f'{attrs} transcript_id "{g["gene_id"]}.1";'
        row = f"{CHR19}\ttoy\t%s\t%d\t%d\t.\t{g['strand']}\t.\t%s"
        lines.append(row % ("gene", g["start"] + 1, g["end"], attrs))
        lines.append(row % ("transcript", g["start"] + 1, g["end"], tx))
        for start, end in g["exons"]:
            lines.append(row % ("exon", start + 1, end, tx))
    return "\n".join(lines) + "\n"


def fasta_text(genome: np.ndarray) -> str:
    out = []
    for name, lo, hi in (
        (CHR19, 0, CHR19_LEN),
        (CHRM, CHR19_LEN, CHR19_LEN + CHRM_LEN),
    ):
        seq = genome[lo:hi].tobytes().decode()
        out.append(f">{name}")
        out.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def generate(
    out_dir: Path, seed: int = DEFAULT_SEED, depth: int = DEFAULT_DEPTH
) -> dict[str, Any]:
    started = time.time()
    out_dir = Path(out_dir)
    (out_dir / "fastq").mkdir(parents=True, exist_ok=True)
    atac_ids, rna_ids = sample_ids("atac"), sample_ids("rna")
    streams = np.random.SeedSequence(seed).spawn(6 + len(atac_ids))
    rng_genome, rng_genes, rng_truth, rng_rna, rng_scale = (
        np.random.default_rng(s) for s in streams[:5]
    )

    genome = make_genome(rng_genome)
    genes = make_genes(rng_genes)
    regions, blacklist = build_truth(rng_truth, genome, genes)
    truth = gene_truth(rng_rna, genes, regions)
    counts = rna_counts(rng_rna, truth)
    depths = [
        round(depth * float(x)) for x in rng_scale.lognormal(0.0, 0.1, len(atac_ids))
    ]
    lfc = np.array([r["log2fc"] for r in regions])
    base = np.array([r["base"] for r in regions])

    hashes: dict[str, str] = {}
    hashes["toy.fa"] = write_text(out_dir / "toy.fa", fasta_text(genome))
    hashes["toy.chrom.sizes"] = write_text(
        out_dir / "toy.chrom.sizes", f"{CHR19}\t{CHR19_LEN}\n{CHRM}\t{CHRM_LEN}\n"
    )
    hashes["toy.gtf"] = write_text(out_dir / "toy.gtf", gtf_text(genes))
    hashes["toy_blacklist.bed"] = write_text(
        out_dir / "toy_blacklist.bed",
        "".join(
            f"{CHR19}\t{s}\t{e}\tblk{i + 1}\n" for i, (s, e) in enumerate(blacklist)
        ),
    )

    ccre_rng = np.random.default_rng(streams[5])
    distal_idx = [i for i, r in enumerate(regions) if r["class"] == "distal"]
    in_ccre = set(
        ccre_rng.choice(
            distal_idx, size=round(CCRE_FRACTION * len(distal_idx)), replace=False
        ).tolist()
    )
    ccre_rows = []
    for i, r in enumerate(regions):
        if r["class"] == "promoter":
            ccre_rows.append(
                f"{CHR19}\t{r['start']}\t{r['end']}\tTOYCCRE{i + 1:04d}\tPLS\n"
            )
        elif i in in_ccre:
            ccre_rows.append(
                f"{CHR19}\t{r['start']}\t{r['end']}\tTOYCCRE{i + 1:04d}\tdELS\n"
            )
        r["in_ccre"] = r["class"] == "promoter" or i in in_ccre
    hashes["toy_ccres.bed"] = write_text(out_dir / "toy_ccres.bed", "".join(ccre_rows))

    header = [
        "region_id",
        "chrom",
        "start",
        "end",
        "class",
        "linked_gene_id",
        "true_log2fc",
        "status",
        "motifs",
        "in_ccre",
    ]
    rows = ["\t".join(header)]
    for r in regions:
        rows.append(
            "\t".join(
                [
                    r["region_id"],
                    CHR19,
                    str(r["start"]),
                    str(r["end"]),
                    r["class"],
                    r["gene_id"],
                    f"{r['log2fc']:.4f}",
                    r["status"],
                    ",".join(r["motifs"]) or "NA",
                    str(r["in_ccre"]).lower(),
                ]
            )
        )
    hashes["truth_regions.tsv"] = write_text(
        out_dir / "truth_regions.tsv", "\n".join(rows) + "\n"
    )

    rows = ["gene_id\tgene_name\ttrue_log2fc\tstatus\tbase_mean"]
    rows += [
        "\t".join(
            [
                g["gene_id"],
                g["gene_name"],
                f"{g['log2fc']:.4f}",
                g["status"],
                f"{g['base_mean']:.2f}",
            ]
        )
        for g in truth
    ]
    hashes["truth_genes.tsv"] = write_text(
        out_dir / "truth_genes.tsv", "\n".join(rows) + "\n"
    )

    rows = ["gene_id\tgene_name\t" + "\t".join(rna_ids)]
    for g, line in zip(truth, counts, strict=True):
        rows.append(
            f"{g['gene_id']}\t{g['gene_name']}\t" + "\t".join(str(int(x)) for x in line)
        )
    hashes["rna_counts.tsv"] = write_text(
        out_dir / "rna_counts.tsv", "\n".join(rows) + "\n"
    )

    samples: dict[str, Any] = {}
    for k, sample in enumerate(atac_ids):
        stimulated = k >= len(REPLICATES)
        rng = np.random.default_rng(streams[6 + k])
        noise = rng.lognormal(0.0, 0.15, len(regions))
        weights = base * noise * (2.0**lfc if stimulated else 1.0)
        frag = simulate_fragments(rng, regions, weights, blacklist, depths[k])
        r1, q1, r2, q2 = extract_reads(rng, genome, frag)
        for mate, (seq, qual) in enumerate(((r1, q1), (r2, q2)), start=1):
            name = f"fastq/{sample}_R{mate}.fastq.gz"
            hashes[name] = write_gz(out_dir / name, fastq_text(sample, mate, seq, qual))
        origin = frag["origin"]
        samples[sample] = {
            "pairs": depths[k],
            "n_unique": int(frag["n_unique"]),
            "n_duplicates": depths[k] - int(frag["n_unique"]),
            "n_chrM": int((frag["chrom_m"] == 1).sum()),
            "n_blacklist": int((origin == 1).sum()),
            "n_region": int((origin == 0).sum()),
            "n_background": int((origin == 2).sum()),
        }

    manifest = {
        "generator": "make_toy.py",
        "generator_version": GENERATOR_VERSION,
        "seed": seed,
        "depth": depth,
        "numpy": np.__version__,
        "sha256_uncompressed": dict(sorted(hashes.items())),
        "atac_samples": samples,
        "rna_samples": rna_ids,
        "truth": {
            "regions": len(regions),
            "promoter": sum(r["class"] == "promoter" for r in regions),
            "distal": sum(r["class"] == "distal" for r in regions),
            "gained": sum(r["status"] == "gained" for r in regions),
            "lost": sum(r["status"] == "lost" for r in regions),
            "unchanged": sum(r["status"] == "unchanged" for r in regions),
            "genes": len(genes),
        },
        "seconds": round(time.time() - started, 1),
    }
    (out_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate the NeuroChrom-X synthetic toy dataset."
    )
    parser.add_argument("--out", type=Path, default=Path("tests/data/toy/generated"))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--depth",
        type=int,
        default=DEFAULT_DEPTH,
        help="mean read pairs per ATAC sample",
    )
    args = parser.parse_args(argv)
    manifest = generate(args.out, args.seed, args.depth)
    print(
        f"wrote {len(manifest['sha256_uncompressed'])} files to {args.out} "
        f"in {manifest['seconds']} s"
    )
    print(json.dumps(manifest["truth"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

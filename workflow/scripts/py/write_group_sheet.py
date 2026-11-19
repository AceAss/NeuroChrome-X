#!/usr/bin/env python3
"""Write the per-assay group sheet (sample_id, condition, batch, ...) for one study.

    python write_group_sheet.py --samples config/samples.tsv --study toy \
        --assay atac --out group_sheet.tsv

The output feeds consensus_peaks.py and deseq2_dar.R.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from neurochromx.io.samplesheet import load_samplesheet

COLUMNS = [
    "sample_id",
    "condition",
    "batch",
    "sex",
    "age_weeks",
    "timepoint_min",
    "replicate",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    sheet = load_samplesheet(args.samples)
    rows = [
        r
        for r in sheet.records
        if r.study_id == args.study and r.assay.value == args.assay
    ]
    if not rows:
        raise SystemExit(f"no {args.assay} samples for study {args.study}")
    with args.out.open("w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            values = [
                r.sample_id,
                r.condition,
                r.batch,
                r.sex.value,
                r.age_weeks,
                r.timepoint_min,
                r.replicate,
            ]
            handle.write(
                "\t".join("NA" if v is None else str(v) for v in values) + "\n"
            )
    print(f"{len(rows)} samples written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

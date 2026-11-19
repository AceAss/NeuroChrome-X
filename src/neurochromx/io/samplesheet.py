"""Sample-sheet validation for NeuroChrom-X (PRD v1.0 section 9, schema v0.1).

``config/samples.tsv`` is the only input the workflow trusts. Every row is checked
against controlled vocabularies and cross-field rules, then the sheet is checked as
a whole (unique IDs, unique replicate slots, group sizes, design sanity).

Conventions
-----------
* Missing values: an empty cell, ``NA``, ``NaN``, ``None``, ``null`` or ``.``
  (case-insensitive). Write ``home_cage`` or ``vehicle`` rather than ``none`` for an
  unstimulated control, because ``none`` is read as missing.
* Blank lines and lines starting with ``#`` are ignored.
* Extra columns are allowed only if their name starts with ``x_``. They are kept in
  :attr:`SampleSheet.extras` and never interpreted.
* ``path_or_uri`` for ``entry_point == "fastq"`` is one path for SE and ``R1;R2``
  for PE. It may be empty when ``srr`` is given (the fetch step resolves the run).
* ``age_weeks`` is in weeks. Convert months before entry (1 month is about 4.35 weeks).
* ``genome_native`` is required for every entry point except ``fastq``: processed
  inputs must say which build they came from (ADR-006).

Command line::

    python -m neurochromx.io.samplesheet config/samples.tsv [--check-paths]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)

if TYPE_CHECKING:
    import pandas as pd

__all__ = [
    "ALLOWED_ENTRY_POINTS",
    "REQUIRED_COLUMNS",
    "SCHEMA_VERSION",
    "Assay",
    "CellPopulation",
    "EntryPoint",
    "GenomeBuild",
    "Issue",
    "Layout",
    "SampleRecord",
    "SampleSheet",
    "SampleSheetError",
    "Sex",
    "Subregion",
    "Target",
    "Tissue",
    "group_key",
    "load_samplesheet",
    "main",
]

SCHEMA_VERSION = "0.1"
EXTENSION_PREFIX = "x_"

REQUIRED_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "study_id",
    "gsm",
    "srr",
    "assay",
    "target",
    "tissue",
    "subregion",
    "cell_population",
    "age_weeks",
    "sex",
    "strain",
    "condition",
    "stimulus",
    "timepoint_min",
    "replicate",
    "batch",
    "layout",
    "read_len",
    "genome_native",
    "entry_point",
    "path_or_uri",
    "notes",
)

NA_TOKENS = frozenset({"", "na", "nan", "none", "null", "."})
URI_PREFIXES = ("http://", "https://", "ftp://", "s3://", "gs://", "sra:")
UNSTIMULATED = frozenset({"home_cage", "vehicle", "saline", "naive", "sham", "control"})


# --------------------------------------------------------------------------- #
# Controlled vocabularies
# --------------------------------------------------------------------------- #
class Assay(StrEnum):
    ATAC = "atac"
    CHIP = "chip"
    RNA = "rna"
    WGBS = "wgbs"
    RRBS = "rrbs"
    MEDIP = "medip"


class Target(StrEnum):
    """Histone marks (and input controls) for ``assay == chip``."""

    H3K27AC = "h3k27ac"
    H3K4ME1 = "h3k4me1"
    H3K4ME3 = "h3k4me3"
    H3K9AC = "h3k9ac"
    H3K27ME3 = "h3k27me3"
    H3K9ME3 = "h3k9me3"
    H3K36ME3 = "h3k36me3"
    INPUT = "input"


class Tissue(StrEnum):
    HIPPOCAMPUS = "hippocampus"
    CORTEX = "cortex"
    MPFC = "mpfc"
    ACC = "acc"
    AMYGDALA = "amygdala"
    WHOLE_BRAIN = "whole_brain"
    OTHER = "other"


class Subregion(StrEnum):
    WHOLE = "whole"
    DG = "dg"
    CA1 = "ca1"
    CA2 = "ca2"
    CA3 = "ca3"
    DORSAL = "dorsal"
    VENTRAL = "ventral"


class CellPopulation(StrEnum):
    BULK = "bulk"
    NEUN_POS = "neun_pos"
    NEUN_NEG = "neun_neg"
    NEURONS = "neurons"
    NON_NEURONS = "non_neurons"
    ENGRAM = "engram"
    NON_ENGRAM = "non_engram"
    OTHER = "other"


class Sex(StrEnum):
    MALE = "male"
    FEMALE = "female"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class Layout(StrEnum):
    SE = "SE"
    PE = "PE"


class EntryPoint(StrEnum):
    FASTQ = "fastq"
    BAM = "bam"
    PEAKS = "peaks"
    COUNTS = "counts"
    COV = "cov"


class GenomeBuild(StrEnum):
    MM9 = "mm9"
    MM10 = "mm10"
    MM39 = "mm39"
    TOY = "toy"


ALLOWED_ENTRY_POINTS: dict[Assay, frozenset[EntryPoint]] = {
    Assay.ATAC: frozenset(
        {EntryPoint.FASTQ, EntryPoint.BAM, EntryPoint.PEAKS, EntryPoint.COUNTS}
    ),
    Assay.CHIP: frozenset(
        {EntryPoint.FASTQ, EntryPoint.BAM, EntryPoint.PEAKS, EntryPoint.COUNTS}
    ),
    Assay.RNA: frozenset({EntryPoint.FASTQ, EntryPoint.BAM, EntryPoint.COUNTS}),
    Assay.WGBS: frozenset({EntryPoint.FASTQ, EntryPoint.BAM, EntryPoint.COV}),
    Assay.RRBS: frozenset({EntryPoint.FASTQ, EntryPoint.BAM, EntryPoint.COV}),
    Assay.MEDIP: frozenset({EntryPoint.FASTQ, EntryPoint.BAM, EntryPoint.COUNTS}),
}

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_\-]{0,63}$")
STUDY_ID_RE = re.compile(r"^[a-z][a-z0-9_]{0,31}$")
SLUG_RE = re.compile(r"^[a-z0-9][a-z0-9_\-]{0,47}$")
BATCH_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_\-]{0,31}$")
STRAIN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-/+]{0,63}$")
GSM_RE = re.compile(r"^GSM\d+$")
SRR_RE = re.compile(r"^[SED]RR\d+$")


def _matches(value: str, pattern: re.Pattern[str], hint: str) -> str:
    if not pattern.fullmatch(value):
        raise ValueError(hint)
    return value


# --------------------------------------------------------------------------- #
# One row
# --------------------------------------------------------------------------- #
class SampleRecord(BaseModel):
    """One validated row of ``samples.tsv``."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True, frozen=True)

    sample_id: str
    study_id: str
    gsm: str | None = None
    srr: str | None = None
    assay: Assay
    target: Target | None = None
    tissue: Tissue
    subregion: Subregion | None = None
    cell_population: CellPopulation
    age_weeks: float | None = Field(default=None, ge=0, le=300)
    sex: Sex = Sex.UNKNOWN
    strain: str | None = None
    condition: str
    stimulus: str | None = None
    timepoint_min: int | None = Field(default=None, ge=0)
    replicate: int = Field(ge=1)
    batch: str
    layout: Layout | None = None
    read_len: int | None = Field(default=None, gt=0, le=1000)
    genome_native: GenomeBuild | None = None
    entry_point: EntryPoint
    path_or_uri: str | None = None
    notes: str | None = None

    # -- normalisation (runs before type validation) ------------------------ #
    @field_validator(
        "assay",
        "target",
        "tissue",
        "subregion",
        "cell_population",
        "sex",
        "entry_point",
        "genome_native",
        mode="before",
    )
    @classmethod
    def _lower_enum(cls, value: Any) -> Any:
        return value.strip().lower() if isinstance(value, str) else value

    @field_validator("layout", mode="before")
    @classmethod
    def _upper_layout(cls, value: Any) -> Any:
        return value.strip().upper() if isinstance(value, str) else value

    @field_validator("timepoint_min", "replicate", "read_len", mode="before")
    @classmethod
    def _integral(cls, value: Any) -> Any:
        """Accept ``"60.0"`` for an integer column; leave anything else to pydantic."""
        if isinstance(value, str):
            try:
                number = float(value.strip())
            except ValueError:
                return value
            if number.is_integer():
                return int(number)
        return value

    # -- patterns ------------------------------------------------------------ #
    @field_validator("sample_id")
    @classmethod
    def _check_sample_id(cls, value: str) -> str:
        return _matches(
            value,
            SAMPLE_ID_RE,
            "use letters, digits, '_' or '-' (no dots: IDs become file names); "
            "start with a letter or digit; at most 64 characters",
        )

    @field_validator("study_id")
    @classmethod
    def _check_study_id(cls, value: str) -> str:
        return _matches(
            value,
            STUDY_ID_RE,
            "use lower-case letters, digits and '_'; start with a letter "
            "(for example su2017)",
        )

    @field_validator("gsm")
    @classmethod
    def _check_gsm(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _matches(value, GSM_RE, "expected GSM followed by digits")

    @field_validator("srr")
    @classmethod
    def _check_srr(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _matches(value, SRR_RE, "expected an SRR, ERR or DRR run accession")

    @field_validator("condition", "stimulus")
    @classmethod
    def _check_slug(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _matches(
            value,
            SLUG_RE,
            "use lower-case letters, digits, '_' or '-' (for example stimulated)",
        )

    @field_validator("batch")
    @classmethod
    def _check_batch(cls, value: str) -> str:
        return _matches(
            value,
            BATCH_RE,
            "use letters, digits, '_' or '-' (write b1 if there is a single batch)",
        )

    @field_validator("strain")
    @classmethod
    def _check_strain(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _matches(value, STRAIN_RE, "unexpected characters in strain")

    # -- cross-field rules --------------------------------------------------- #
    @model_validator(mode="after")
    def _cross_field_rules(self) -> SampleRecord:
        problems: list[str] = []

        if self.assay is Assay.CHIP and self.target is None:
            problems.append(
                "target is required for chip "
                "(for example h3k27ac, or input for a control)"
            )
        if self.assay is not Assay.CHIP and self.target is not None:
            problems.append(f"target must be empty for assay {self.assay.value}")

        if self.entry_point not in ALLOWED_ENTRY_POINTS[self.assay]:
            allowed = ", ".join(
                sorted(e.value for e in ALLOWED_ENTRY_POINTS[self.assay])
            )
            problems.append(
                f"entry_point {self.entry_point.value} is not valid for assay "
                f"{self.assay.value} (allowed: {allowed})"
            )

        if self.entry_point is EntryPoint.FASTQ:
            if self.layout is None:
                problems.append("layout (SE or PE) is required for fastq")
            if self.srr is None and self.path_or_uri is None:
                problems.append("fastq needs srr or path_or_uri")
            if self.path_or_uri is not None and self.layout is not None:
                n_files = len([p for p in self.path_or_uri.split(";") if p.strip()])
                expected = 2 if self.layout is Layout.PE else 1
                if n_files != expected:
                    problems.append(
                        f"{self.layout.value} fastq needs {expected} path(s) separated "
                        f"by ';' (got {n_files})"
                    )
        else:
            if self.path_or_uri is None:
                problems.append(
                    f"path_or_uri is required for entry_point {self.entry_point.value}"
                )
            if self.genome_native is None:
                problems.append(
                    f"genome_native is required for entry_point "
                    f"{self.entry_point.value}"
                )
            if self.entry_point is EntryPoint.BAM and self.layout is None:
                problems.append("layout (SE or PE) is required for bam")

        if problems:
            raise ValueError("; ".join(problems))
        return self


# --------------------------------------------------------------------------- #
# Whole sheet
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Issue:
    """One validation problem, with enough context to fix it."""

    row: int | None
    sample_id: str | None
    field: str | None
    message: str

    def __str__(self) -> str:
        where = f"row {self.row}" if self.row is not None else "sheet"
        who = f" [{self.sample_id}]" if self.sample_id else ""
        what = f" {self.field}:" if self.field else ""
        return f"{where}{who}{what} {self.message}"


class SampleSheetError(ValueError):
    """Raised with every problem found, not just the first."""

    def __init__(self, issues: Sequence[Issue]) -> None:
        self.issues: tuple[Issue, ...] = tuple(issues)
        shown = "\n".join(f"  - {issue}" for issue in self.issues[:20])
        more = (
            f"\n  ... and {len(self.issues) - 20} more" if len(self.issues) > 20 else ""
        )
        super().__init__(
            f"sample sheet has {len(self.issues)} problem(s):\n{shown}{more}"
        )


GroupKey = tuple[Any, ...]


def group_key(record: SampleRecord) -> GroupKey:
    """Replicates of one biological group share this key (replicate is excluded)."""
    return (
        record.study_id,
        record.assay.value,
        record.target.value if record.target else None,
        record.tissue.value,
        record.subregion.value if record.subregion else None,
        record.cell_population.value,
        record.sex.value,
        record.age_weeks,
        record.condition,
        record.stimulus,
        record.timepoint_min,
    )


@dataclass(frozen=True)
class SampleSheet:
    """A validated sample sheet."""

    records: tuple[SampleRecord, ...]
    warnings: tuple[str, ...] = ()
    extras: dict[str, dict[str, str | None]] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.records)

    def by_study(self) -> dict[str, list[SampleRecord]]:
        out: dict[str, list[SampleRecord]] = defaultdict(list)
        for record in self.records:
            out[record.study_id].append(record)
        return dict(out)

    def by_assay(self, assay: Assay | str) -> list[SampleRecord]:
        wanted = Assay(assay)
        return [r for r in self.records if r.assay is wanted]

    def groups(self) -> dict[GroupKey, list[SampleRecord]]:
        out: dict[GroupKey, list[SampleRecord]] = defaultdict(list)
        for record in self.records:
            out[group_key(record)].append(record)
        return dict(out)

    def to_dataframe(self) -> pd.DataFrame:
        """Return the sheet as a DataFrame (needs pandas, which is optional)."""
        import pandas as pd

        rows = [record.model_dump(mode="json") for record in self.records]
        return pd.DataFrame(rows, columns=list(REQUIRED_COLUMNS))


def _clean(cell: str) -> str | None:
    text = cell.strip()
    return None if text.lower() in NA_TOKENS else text


def _read_table(path: Path) -> tuple[list[str], list[tuple[int, list[str]]]]:
    try:
        text = path.read_text(encoding="utf-8-sig")
    except OSError as exc:
        raise SampleSheetError(
            [Issue(None, None, None, f"cannot read {path}: {exc}")]
        ) from exc

    header: list[str] | None = None
    rows: list[tuple[int, list[str]]] = []
    for line_no, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        cells = next(csv.reader([line], delimiter="\t"))
        if header is None:
            header = [cell.strip() for cell in cells]
        else:
            rows.append((line_no, cells))
    if header is None:
        raise SampleSheetError([Issue(None, None, None, "file has no header line")])
    return header, rows


def _check_header(header: Sequence[str]) -> list[Issue]:
    issues: list[Issue] = []
    seen: set[str] = set()
    for name in header:
        if name in seen:
            issues.append(Issue(None, None, name, "duplicate column"))
        seen.add(name)
    missing = [c for c in REQUIRED_COLUMNS if c not in seen]
    if missing:
        issues.append(
            Issue(None, None, None, f"missing column(s): {', '.join(missing)}")
        )
    unknown = [
        c
        for c in header
        if c not in REQUIRED_COLUMNS and not c.startswith(EXTENSION_PREFIX)
    ]
    if unknown:
        issues.append(
            Issue(
                None,
                None,
                None,
                f"unknown column(s): {', '.join(unknown)} "
                f"(extra columns must start with '{EXTENSION_PREFIX}')",
            )
        )
    return issues


def _issues_from_validation(
    exc: ValidationError, row: int, sample_id: str | None
) -> list[Issue]:
    issues: list[Issue] = []
    for err in exc.errors():
        name = ".".join(str(part) for part in err["loc"]) or None
        message = str(err["msg"])
        if err["type"] == "missing":
            message = "is required (empty or NA)"
        elif message.startswith("Value error, "):
            message = message[len("Value error, ") :]
        issues.append(Issue(row, sample_id, name, message))
    return issues


def _describe(key: GroupKey) -> str:
    study, assay, target, tissue, subregion, cell_pop, sex, age, cond, stim, tp = key
    parts = [f"{study}/{assay}" + (f"/{target}" if target else "")]
    parts.append(f"{tissue}" + (f"-{subregion}" if subregion else "") + f" {cell_pop}")
    parts.append(f"condition={cond}")
    if stim:
        parts.append(f"stimulus={stim}")
    if tp is not None:
        parts.append(f"t={tp}min")
    parts.append(f"sex={sex}")
    if age is not None:
        parts.append(f"age={age:g}w")
    return " ".join(parts)


def _collect_warnings(records: Sequence[SampleRecord]) -> list[str]:
    warnings: list[str] = []
    groups: dict[GroupKey, list[SampleRecord]] = defaultdict(list)
    for record in records:
        groups[group_key(record)].append(record)

    for key, members in groups.items():
        if len(members) < 2:
            warnings.append(
                f"group [{_describe(key)}] has {len(members)} sample; "
                "PRD section 6 asks for at least 2 per group"
            )

    units: dict[tuple[Any, ...], set[GroupKey]] = defaultdict(set)
    for key in groups:
        units[key[:3]].add(key)
    for (study, assay, target), keys in units.items():
        if len(keys) < 2:
            label = f"{study}/{assay}" + (f"/{target}" if target else "")
            warnings.append(
                f"{label} has a single group, so no contrast can be computed"
            )

    no_timepoint: dict[str, list[str]] = defaultdict(list)
    needs_liftover: dict[tuple[str, str], list[str]] = defaultdict(list)
    for record in records:
        stimulated = record.stimulus is not None and record.stimulus not in UNSTIMULATED
        if stimulated and record.timepoint_min is None:
            no_timepoint[record.study_id].append(record.sample_id)
        if record.entry_point is not EntryPoint.FASTQ and record.genome_native not in (
            GenomeBuild.MM10,
            GenomeBuild.TOY,
        ):
            build = record.genome_native.value if record.genome_native else "unknown"
            needs_liftover[(record.study_id, build)].append(record.sample_id)
    for study, ids in no_timepoint.items():
        warnings.append(
            f"study {study}: {len(ids)} stimulated sample(s) have no timepoint_min "
            f"(first: {ids[0]})"
        )
    for (study, build), ids in needs_liftover.items():
        warnings.append(
            f"study {study}: {len(ids)} processed input(s) are on {build}; "
            "lift over to mm10 or exclude (ADR-006)"
        )
    return warnings


def _check_paths(
    records: Sequence[SampleRecord], base_dir: Path, rows: dict[str, int]
) -> list[Issue]:
    issues: list[Issue] = []
    for record in records:
        if not record.path_or_uri:
            continue
        for part in (p.strip() for p in record.path_or_uri.split(";") if p.strip()):
            if part.lower().startswith(URI_PREFIXES):
                continue
            candidate = Path(part)
            if not candidate.is_absolute():
                candidate = base_dir / candidate
            if not candidate.exists():
                issues.append(
                    Issue(
                        rows.get(record.sample_id),
                        record.sample_id,
                        "path_or_uri",
                        f"file not found: {part}",
                    )
                )
    return issues


def load_samplesheet(
    path: str | Path,
    *,
    check_paths: bool = False,
    base_dir: str | Path | None = None,
) -> SampleSheet:
    """Read and validate ``samples.tsv``.

    Raises :class:`SampleSheetError` listing every problem. Soft problems (small
    groups, missing time points, processed inputs on another build) are returned in
    :attr:`SampleSheet.warnings` instead.
    """
    header, raw_rows = _read_table(Path(path))
    header_issues = _check_header(header)
    if header_issues:
        raise SampleSheetError(header_issues)

    issues: list[Issue] = []
    records: list[SampleRecord] = []
    extras: dict[str, dict[str, str | None]] = {}
    first_row_of: dict[str, int] = {}

    for line_no, cells in raw_rows:
        if len(cells) > len(header):
            issues.append(
                Issue(
                    line_no,
                    None,
                    None,
                    f"row has {len(cells)} cells but the header has {len(header)}",
                )
            )
            continue
        padded = [*cells, *([""] * (len(header) - len(cells)))]
        cleaned = {
            name: _clean(cell) for name, cell in zip(header, padded, strict=True)
        }
        sample_id = cleaned.get("sample_id")

        if sample_id is not None:
            if sample_id in first_row_of:
                issues.append(
                    Issue(
                        line_no,
                        sample_id,
                        "sample_id",
                        f"duplicate (first seen on row {first_row_of[sample_id]})",
                    )
                )
                continue
            first_row_of[sample_id] = line_no

        payload = {
            name: value
            for name, value in cleaned.items()
            if name in REQUIRED_COLUMNS and value is not None
        }
        try:
            record = SampleRecord.model_validate(payload)
        except ValidationError as exc:
            issues.extend(_issues_from_validation(exc, line_no, sample_id))
            continue
        records.append(record)
        extras[record.sample_id] = {
            name: value
            for name, value in cleaned.items()
            if name.startswith(EXTENSION_PREFIX)
        }

    slots: dict[tuple[Any, ...], str] = {}
    for record in records:
        slot = (*group_key(record), record.replicate)
        if slot in slots:
            issues.append(
                Issue(
                    first_row_of.get(record.sample_id),
                    record.sample_id,
                    "replicate",
                    f"same group and replicate number as {slots[slot]}",
                )
            )
        else:
            slots[slot] = record.sample_id

    if check_paths:
        root = Path(base_dir) if base_dir is not None else Path.cwd()
        issues.extend(_check_paths(records, root, first_row_of))

    if issues:
        raise SampleSheetError(issues)

    return SampleSheet(
        records=tuple(records),
        warnings=tuple(_collect_warnings(records)),
        extras=extras,
    )


# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #
def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m neurochromx.io.samplesheet",
        description="Validate a NeuroChrom-X sample sheet (schema "
        f"{SCHEMA_VERSION}). Exit code 0 means valid; warnings do not fail.",
    )
    parser.add_argument("path", type=Path, help="path to samples.tsv")
    parser.add_argument(
        "--check-paths",
        action="store_true",
        help="also require local files in path_or_uri to exist",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=None,
        help="directory that relative paths are resolved against (default: cwd)",
    )
    args = parser.parse_args(argv)

    try:
        sheet = load_samplesheet(
            args.path, check_paths=args.check_paths, base_dir=args.base_dir
        )
    except SampleSheetError as exc:
        print(exc, file=sys.stderr)
        return 1

    groups = sheet.groups()
    print(
        f"OK: {len(sheet)} samples, {len(sheet.by_study())} study(ies), "
        f"{len(groups)} group(s) (schema {SCHEMA_VERSION})"
    )
    for key, members in groups.items():
        print(f"  n={len(members)}  {_describe(key)}")
    for message in sheet.warnings:
        print(f"warning: {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

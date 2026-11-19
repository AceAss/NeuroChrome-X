"""Tests for neurochromx.io.samplesheet (CL-003)."""

from __future__ import annotations

from pathlib import Path

import pytest

from neurochromx.io.samplesheet import (
    REQUIRED_COLUMNS,
    Assay,
    EntryPoint,
    SampleSheetError,
    load_samplesheet,
    main,
)

TOY_SHEET = Path(__file__).resolve().parents[1] / "data" / "toy" / "samples.tsv"

BASE_ROW = {
    "sample_id": "s1",
    "study_id": "toy",
    "gsm": "",
    "srr": "",
    "assay": "atac",
    "target": "",
    "tissue": "hippocampus",
    "subregion": "dg",
    "cell_population": "bulk",
    "age_weeks": "10",
    "sex": "male",
    "strain": "C57BL/6J",
    "condition": "control",
    "stimulus": "home_cage",
    "timepoint_min": "0",
    "replicate": "1",
    "batch": "b1",
    "layout": "PE",
    "read_len": "50",
    "genome_native": "toy",
    "entry_point": "fastq",
    "path_or_uri": "a_R1.fastq.gz;a_R2.fastq.gz",
    "notes": "",
}


def row(**overrides: str) -> dict[str, str]:
    return {**BASE_ROW, **overrides}


def two_groups(**shared: str) -> list[dict[str, str]]:
    """Two groups of two replicates each, so a sheet is valid and warning-free."""
    rows = []
    for cond, stim, tp in (
        ("control", "home_cage", "0"),
        ("stimulated", "synthetic", "60"),
    ):
        for rep in (1, 2):
            sid = f"{cond[:4]}_{rep}"
            rows.append(
                row(
                    sample_id=sid,
                    condition=cond,
                    stimulus=stim,
                    timepoint_min=tp,
                    replicate=str(rep),
                    path_or_uri=f"{sid}_R1.fastq.gz;{sid}_R2.fastq.gz",
                    **shared,
                )
            )
    return rows


def write_sheet(
    path: Path,
    rows: list[dict[str, str]],
    header: list[str] | None = None,
    prefix: str = "",
) -> Path:
    columns = header if header is not None else list(REQUIRED_COLUMNS)
    lines = ["\t".join(columns)]
    lines += ["\t".join(r.get(c, "") for c in columns) for r in rows]
    path.write_text(prefix + "\n".join(lines) + "\n", encoding="utf-8")
    return path


def problems(tmp_path: Path, rows: list[dict[str, str]], **kwargs: object) -> str:
    sheet = write_sheet(tmp_path / "samples.tsv", rows)
    with pytest.raises(SampleSheetError) as caught:
        load_samplesheet(sheet, **kwargs)  # type: ignore[arg-type]
    return str(caught.value)


# --------------------------------------------------------------------------- #
# Happy paths
# --------------------------------------------------------------------------- #
def test_toy_sheet_loads_without_warnings() -> None:
    sheet = load_samplesheet(TOY_SHEET)
    assert len(sheet) == 12
    assert sheet.warnings == ()
    assert {r.assay for r in sheet.records} == {Assay.ATAC, Assay.RNA}
    assert len(sheet.groups()) == 4
    assert all(len(members) == 3 for members in sheet.groups().values())
    assert len(sheet.by_assay("rna")) == 6
    assert set(sheet.by_study()) == {"toy"}


def test_valid_two_group_sheet(tmp_path: Path) -> None:
    sheet = load_samplesheet(write_sheet(tmp_path / "s.tsv", two_groups()))
    assert len(sheet) == 4
    assert sheet.warnings == ()


def test_na_tokens_and_numeric_coercion(tmp_path: Path) -> None:
    rows = two_groups()
    rows[0].update(age_weeks="10.5", timepoint_min="0.0", strain="NA", subregion=".")
    rows[1].update(strain="None", notes="null", read_len="50.0")
    record = load_samplesheet(write_sheet(tmp_path / "s.tsv", rows)).records[0]
    assert record.age_weeks == 10.5
    assert record.timepoint_min == 0
    assert record.strain is None
    assert record.subregion is None


def test_enum_case_is_normalised(tmp_path: Path) -> None:
    rows = two_groups()
    rows[0].update(assay="ATAC", layout="pe", sex="Male", tissue="Hippocampus")
    record = load_samplesheet(write_sheet(tmp_path / "s.tsv", rows)).records[0]
    assert record.assay is Assay.ATAC
    assert record.layout is not None and record.layout.value == "PE"


def test_header_extras_comments_bom_and_blank_lines(tmp_path: Path) -> None:
    header = [*REQUIRED_COLUMNS, "x_lab_note"]
    rows = [{**r, "x_lab_note": "kept"} for r in two_groups()]
    sheet_path = write_sheet(
        tmp_path / "s.tsv", rows, header=header, prefix="\ufeff# comment line\n\n"
    )
    sheet = load_samplesheet(sheet_path)
    assert len(sheet) == 4
    assert sheet.extras["cont_1"] == {"x_lab_note": "kept"}


def test_srr_only_fastq_is_valid(tmp_path: Path) -> None:
    rows = two_groups()
    for i, r in enumerate(rows):
        r.update(srr=f"SRR{1000 + i}", path_or_uri="")
    assert len(load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))) == 4


def test_chip_with_target_and_processed_entry_point(tmp_path: Path) -> None:
    rows = two_groups(
        assay="chip", target="h3k27ac", entry_point="peaks", genome_native="mm10"
    )
    for r in rows:
        r.update(layout="", read_len="", path_or_uri=f"{r['sample_id']}.narrowPeak")
    sheet = load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))
    assert sheet.records[0].entry_point is EntryPoint.PEAKS
    assert sheet.warnings == ()


def test_to_dataframe_has_schema_columns() -> None:
    pd = pytest.importorskip("pandas")
    frame = load_samplesheet(TOY_SHEET).to_dataframe()
    assert isinstance(frame, pd.DataFrame)
    assert list(frame.columns) == list(REQUIRED_COLUMNS)
    assert len(frame) == 12


# --------------------------------------------------------------------------- #
# Row-level errors
# --------------------------------------------------------------------------- #
def test_invalid_enum_value(tmp_path: Path) -> None:
    rows = two_groups()
    rows[0]["assay"] = "atacseq"
    message = problems(tmp_path, rows)
    assert "assay" in message and "row 2" in message


def test_chip_requires_target_and_atac_forbids_it(tmp_path: Path) -> None:
    chip = two_groups(assay="chip")
    assert "target is required for chip" in problems(tmp_path, chip)
    atac = two_groups()
    atac[0]["target"] = "h3k27ac"
    assert "target must be empty" in problems(tmp_path, atac)


def test_fastq_path_rules(tmp_path: Path) -> None:
    pe_one_file = two_groups()
    pe_one_file[0]["path_or_uri"] = "only_R1.fastq.gz"
    assert "needs 2 path(s)" in problems(tmp_path, pe_one_file)

    se_two_files = two_groups(layout="SE")
    assert "needs 1 path(s)" in problems(tmp_path, se_two_files)

    nothing = two_groups()
    nothing[0].update(path_or_uri="", srr="")
    assert "fastq needs srr or path_or_uri" in problems(tmp_path, nothing)

    no_layout = two_groups()
    no_layout[0]["layout"] = ""
    assert "layout (SE or PE) is required" in problems(tmp_path, no_layout)


def test_processed_input_requires_genome_and_path(tmp_path: Path) -> None:
    rows = two_groups(assay="rna", entry_point="counts")
    for r in rows:
        r.update(layout="", read_len="", path_or_uri="counts.tsv")
    rows[0]["genome_native"] = ""
    assert "genome_native is required" in problems(tmp_path, rows)

    rows = two_groups(assay="rna", entry_point="counts")
    rows[0].update(layout="", path_or_uri="")
    assert "path_or_uri is required" in problems(tmp_path, rows)


def test_assay_entry_point_compatibility(tmp_path: Path) -> None:
    rows = two_groups(assay="rna", entry_point="peaks", genome_native="mm10")
    for r in rows:
        r["path_or_uri"] = "x.bed"
    message = problems(tmp_path, rows)
    assert "not valid for assay rna" in message


def test_bad_identifiers(tmp_path: Path) -> None:
    rows = two_groups()
    rows[0].update(sample_id="bad id", study_id="Study 1", gsm="GSE1", srr="XYZ1")
    message = problems(tmp_path, rows)
    for name in ("sample_id", "study_id", "gsm", "srr"):
        assert name in message


def test_required_value_missing(tmp_path: Path) -> None:
    rows = two_groups()
    rows[0]["condition"] = "NA"
    message = problems(tmp_path, rows)
    assert "condition" in message and "is required" in message


def test_all_bad_rows_are_reported_with_row_numbers(tmp_path: Path) -> None:
    rows = two_groups()
    rows[0]["assay"] = "bogus"
    rows[2]["replicate"] = "0"
    with pytest.raises(SampleSheetError) as caught:
        load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))
    assert {issue.row for issue in caught.value.issues} == {2, 4}


# --------------------------------------------------------------------------- #
# Sheet-level errors and warnings
# --------------------------------------------------------------------------- #
def test_header_problems(tmp_path: Path) -> None:
    missing = [c for c in REQUIRED_COLUMNS if c != "batch"]
    path = write_sheet(tmp_path / "a.tsv", two_groups(), header=missing)
    with pytest.raises(SampleSheetError, match="missing column"):
        load_samplesheet(path)

    unknown = [*REQUIRED_COLUMNS, "surprise"]
    path = write_sheet(tmp_path / "b.tsv", two_groups(), header=unknown)
    with pytest.raises(SampleSheetError, match="unknown column"):
        load_samplesheet(path)


def test_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("# nothing here\n", encoding="utf-8")
    with pytest.raises(SampleSheetError, match="no header"):
        load_samplesheet(path)


def test_duplicate_sample_id(tmp_path: Path) -> None:
    rows = two_groups()
    rows[1]["sample_id"] = rows[0]["sample_id"]
    assert "duplicate (first seen on row 2)" in problems(tmp_path, rows)


def test_duplicate_replicate_slot(tmp_path: Path) -> None:
    rows = two_groups()
    rows[1]["replicate"] = rows[0]["replicate"]
    assert "same group and replicate number" in problems(tmp_path, rows)


def test_singleton_group_and_single_group_warn(tmp_path: Path) -> None:
    rows = [row(sample_id="only")]
    sheet = load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))
    text = " ".join(sheet.warnings)
    assert "has 1 sample" in text
    assert "no contrast can be computed" in text


def test_stimulated_without_timepoint_warns(tmp_path: Path) -> None:
    rows = two_groups()
    for r in rows[2:]:
        r["timepoint_min"] = ""
    sheet = load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))
    assert any("no timepoint_min" in w for w in sheet.warnings)


def test_processed_input_on_other_build_warns(tmp_path: Path) -> None:
    rows = two_groups(assay="rna", entry_point="counts", genome_native="mm9")
    for r in rows:
        r.update(layout="", read_len="", path_or_uri="counts.tsv")
    sheet = load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))
    assert any("mm9" in w and "ADR-006" in w for w in sheet.warnings)


def test_check_paths(tmp_path: Path) -> None:
    rows = two_groups()
    sheet_path = write_sheet(tmp_path / "s.tsv", rows)
    with pytest.raises(SampleSheetError, match="file not found"):
        load_samplesheet(sheet_path, check_paths=True, base_dir=tmp_path)

    for r in rows:
        for part in r["path_or_uri"].split(";"):
            (tmp_path / part).write_text("", encoding="utf-8")
    assert len(load_samplesheet(sheet_path, check_paths=True, base_dir=tmp_path)) == 4


def test_uri_paths_are_not_checked(tmp_path: Path) -> None:
    rows = two_groups()
    for r in rows:
        r["path_or_uri"] = (
            "https://example.org/a_R1.fastq.gz;https://example.org/a_R2.fastq.gz"
        )
    sheet_path = write_sheet(tmp_path / "s.tsv", rows)
    assert len(load_samplesheet(sheet_path, check_paths=True, base_dir=tmp_path)) == 4


# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #
def test_cli_ok_and_failure(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main([str(TOY_SHEET)]) == 0
    out = capsys.readouterr().out
    assert out.startswith("OK: 12 samples")

    bad = write_sheet(tmp_path / "bad.tsv", [row(assay="nope")])
    assert main([str(bad)]) == 1
    assert "sample sheet has" in capsys.readouterr().err

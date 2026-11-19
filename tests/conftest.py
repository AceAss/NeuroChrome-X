"""Shared pytest configuration and fixtures."""

import importlib.util
from pathlib import Path

import pytest

TOY_DIR = Path(__file__).resolve().parent / "data" / "toy"


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: runs real tools (bowtie2, samtools, bedtools, macs3, Rscript, "
        "snakemake); skipped when a tool is missing. "
        "Deselect with -m 'not integration'.",
    )


@pytest.fixture(scope="session")
def toy_small(tmp_path_factory):
    """A tiny generated toy dataset (depth 500) for tests that only need annotations."""
    spec = importlib.util.spec_from_file_location("make_toy", TOY_DIR / "make_toy.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    out = tmp_path_factory.mktemp("toy_small")
    module.generate(out, seed=20261105, depth=500)
    return out

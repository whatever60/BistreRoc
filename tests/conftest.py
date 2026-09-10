"""Shared deterministic expression fixtures."""

from pathlib import Path

import pytest


@pytest.fixture
def labeled_reference(tmp_path: Path) -> Path:
    """Create a small reference with two strong markers for three cell types."""

    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "Gene\tA type\tA type\tA type\tB+ type\tB+ type\tB+ type\tC type\tC type\tC type\n"
        "A1\t120\t110\t130\t1\t1\t1\t2\t2\t2\n"
        "A2\t90\t85\t95\t1\t1\t1\t2\t2\t2\n"
        "B1\t1\t1\t1\t120\t110\t130\t2\t2\t2\n"
        "B2\t1\t1\t1\t90\t85\t95\t2\t2\t2\n"
        "C1\t1\t1\t1\t2\t2\t2\t120\t110\t130\n"
        "C2\t1\t1\t1\t2\t2\t2\t90\t85\t95\n"
        "Housekeeping\t20\t20\t20\t20\t20\t20\t20\t20\t20\n",
        encoding="utf-8",
    )
    return reference


@pytest.fixture
def deconvolution_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """Create a signature and two mixtures with known dominant cell types."""

    signature = tmp_path / "signature.tsv"
    signature.write_text(
        "NAME\tA\tB\tC\n"
        "A1\t100\t1\t1\n"
        "A2\t80\t1\t1\n"
        "B1\t1\t100\t1\n"
        "B2\t1\t80\t1\n"
        "C1\t1\t1\t100\n"
        "C2\t1\t1\t80\n",
        encoding="utf-8",
    )
    mixture = tmp_path / "mixture.tsv"
    mixture.write_text(
        "Gene\tmostly_a\tmostly_c\n"
        "A1\t80.2\t10.9\n"
        "A2\t64.2\t8.9\n"
        "B1\t10.9\t10.9\n"
        "B2\t8.9\t8.9\n"
        "C1\t10.9\t80.2\n"
        "C2\t8.9\t64.2\n",
        encoding="utf-8",
    )
    return signature, mixture

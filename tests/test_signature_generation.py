"""Integration tests for standalone signature generation."""

import gzip
from pathlib import Path

import pytest

from bistreroc.compressed_io import open_text_reader
from bistreroc.signature_generation import build_signature
from bistreroc.sparse_reference import dense_integer_tsv_to_sparse_bundle


@pytest.mark.parametrize(
    ("minimum_marker_count", "maximum_marker_count", "q_value_threshold", "message"),
    [
        (0, 1, 0.01, "minimum_markers_per_cell_type must be positive"),
        (
            2,
            1,
            0.01,
            "maximum_markers_per_cell_type cannot be smaller than",
        ),
        (1, 2, 1.01, "q_value_threshold must be between zero and one"),
    ],
)
def test_signature_generation_validates_marker_search_controls(
    tmp_path: Path,
    minimum_marker_count: int,
    maximum_marker_count: int,
    q_value_threshold: float,
    message: str,
) -> None:
    """Reject invalid marker-selection controls before reading the reference."""

    with pytest.raises(ValueError, match=message):
        build_signature(
            reference_path=tmp_path / "missing_reference.tsv",
            output_path=tmp_path / "signature.tsv",
            minimum_markers_per_cell_type=minimum_marker_count,
            maximum_markers_per_cell_type=maximum_marker_count,
            q_value_threshold=q_value_threshold,
        )


def test_dense_and_sparse_references_produce_the_same_signature(
    labeled_reference: Path,
    tmp_path: Path,
) -> None:
    """Generate identical signatures from dense and sparse reference storage."""

    sparse_reference = tmp_path / "reference.bundle"
    dense_integer_tsv_to_sparse_bundle(labeled_reference, sparse_reference)
    dense_output = tmp_path / "dense_signature.tsv"
    sparse_output = tmp_path / "sparse_signature.tsv.gz"
    build_signature(
        labeled_reference,
        dense_output,
        minimum_markers_per_cell_type=1,
        maximum_markers_per_cell_type=1,
        q_value_threshold=1,
        core_count=2,
    )
    build_signature(
        sparse_reference,
        sparse_output,
        minimum_markers_per_cell_type=1,
        maximum_markers_per_cell_type=1,
        q_value_threshold=1,
        core_count=2,
    )

    with open_text_reader(sparse_output) as input_handle:
        sparse_text = input_handle.read()
    assert sparse_text == dense_output.read_text(encoding="utf-8")
    assert sparse_text.splitlines()[0] == "NAME\tA type\tB+ type\tC type"
    assert [line.split("\t", 1)[0] for line in sparse_text.splitlines()[1:]] == [
        "A1",
        "B1",
        "C1",
    ]
    assert gzip.decompress(sparse_output.read_bytes()).decode() == sparse_text


def test_signature_generation_ignores_undefined_welch_degrees_of_freedom(
    tmp_path: Path,
) -> None:
    """Map guarded Welch degrees of freedom to non-significant p-values."""

    reference_path = tmp_path / "guarded_welch_reference.tsv"
    reference_path.write_text(
        "Gene\tAlpha +\tBeta B\tGamma/C\n"
        "alpha_marker\t100000\t1000\t1014\n"
        "beta_marker\t1014\t100000\t1000\n"
        "gamma_marker\t1000\t1014\t100000\n"
        "filler\t897986\t897986\t897986\n"
        "zero\t0\t0\t0\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "signature.tsv"

    build_signature(
        reference_path,
        output_path,
        minimum_markers_per_cell_type=1,
        maximum_markers_per_cell_type=1,
        q_value_threshold=1,
    )

    assert output_path.read_text(encoding="utf-8") == (
        "NAME\tAlpha +\tBeta B\tGamma/C\n"
    )

"""Tests for dense and sparse reference handling."""

import gzip
from pathlib import Path

import numpy as np
import pytest

from bistreroc.sparse_reference import (
    dense_integer_tsv_to_sparse_bundle,
    load_reference_means,
    sparse_bundle_bytes,
)


def test_sparse_conversion_is_deterministic_and_numerically_lossless(
    labeled_reference: Path,
    tmp_path: Path,
) -> None:
    """Produce byte-stable bundles whose cell-type means equal dense input."""

    first_bundle = tmp_path / "first.bundle"
    second_bundle = tmp_path / "second.bundle"
    dense_integer_tsv_to_sparse_bundle(labeled_reference, first_bundle)
    dense_integer_tsv_to_sparse_bundle(labeled_reference, second_bundle)

    first_files = sorted(first_bundle.iterdir())
    second_files = sorted(second_bundle.iterdir())
    assert [path.name for path in first_files] == [path.name for path in second_files]
    assert [path.read_bytes() for path in first_files] == [
        path.read_bytes() for path in second_files
    ]
    dense_genes, dense_types, dense_means = load_reference_means(labeled_reference)
    sparse_genes, sparse_types, sparse_means = load_reference_means(first_bundle)
    assert sparse_genes == dense_genes
    assert sparse_types == dense_types
    np.testing.assert_allclose(sparse_means, dense_means, rtol=1e-14, atol=1e-10)
    assert sparse_bundle_bytes(first_bundle) == sum(
        path.stat().st_size for path in first_files
    )


def test_dense_gzip_reference_has_identical_means(
    labeled_reference: Path,
    tmp_path: Path,
) -> None:
    """Read compressed dense counts without creating an uncompressed copy."""

    compressed_reference = tmp_path / "reference.tsv.gz"
    compressed_reference.write_bytes(gzip.compress(labeled_reference.read_bytes()))
    plain_result = load_reference_means(labeled_reference)
    compressed_result = load_reference_means(compressed_reference)

    assert compressed_result[:2] == plain_result[:2]
    np.testing.assert_array_equal(compressed_result[2], plain_result[2])


def test_reference_conversion_rejects_duplicate_gene_names(tmp_path: Path) -> None:
    """Reject ambiguous row identifiers instead of silently merging genes."""

    duplicate_reference_path = tmp_path / "duplicate_reference.tsv"
    duplicate_reference_path.write_text(
        "Gene\tA\tB\nG1\t1\t2\nG1\t3\t4\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate gene identifier: G1"):
        dense_integer_tsv_to_sparse_bundle(
            duplicate_reference_path,
            tmp_path / "duplicate_reference.bundle",
        )


def test_reference_conversion_rejects_zero_library_size(tmp_path: Path) -> None:
    """Require every single-cell column to contain at least one count."""

    empty_cell_reference_path = tmp_path / "empty_cell_reference.tsv"
    empty_cell_reference_path.write_text(
        "Gene\tA\tB\nG1\t1\t0\nG2\t3\t0\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="each cell library size must be positive"):
        dense_integer_tsv_to_sparse_bundle(
            empty_cell_reference_path,
            tmp_path / "empty_cell_reference.bundle",
        )

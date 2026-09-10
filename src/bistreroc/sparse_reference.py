"""Sparse reference conversion and bounded-memory cell-type aggregation."""

import bz2
import csv
import gzip
import io
import lzma
from collections.abc import Iterator
from contextlib import contextmanager
from itertools import islice
from pathlib import Path
from typing import TextIO

import numpy as np

from bistreroc.compressed_io import open_text_reader

SPARSE_REFERENCE_METADATA_FILENAMES = (
    "genes.tsv.gz",
    "labels.tsv.gz",
    "library_sizes.tsv.gz",
)
MATRIX_COMPRESSION_SUFFIX_BY_FORMAT = {
    "gzip": ".gz",
    "bzip2": ".bz2",
    "xz": ".xz",
}
MATRIX_COORDINATE_CHUNK_SIZE = 500_000


@contextmanager
def open_deterministic_gzip_text_writer(path: Path) -> Iterator[TextIO]:
    """Open a deterministic gzip text writer with no filename or timestamp."""

    with (
        path.open("wb") as raw_handle,
        gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            mtime=0,
        ) as gzip_handle,
        io.TextIOWrapper(gzip_handle, encoding="utf-8") as text_handle,
    ):
        yield text_handle


def parse_dense_count_row(row: list[str], expected_column_count: int) -> np.ndarray:
    """Parse and validate one dense nonnegative integer count row."""

    if len(row) != expected_column_count + 1:
        raise ValueError("reference row width does not match its header")
    values = np.asarray(row[1:], dtype=np.int64)
    if np.any(values < 0):
        raise ValueError("single-cell counts cannot be negative")
    return values


@contextmanager
def open_compressed_matrix_writer(
    path: Path,
    compression_format: str,
) -> Iterator[TextIO]:
    """Open a deterministic text writer for one compressed coordinate matrix."""

    if compression_format == "gzip":
        with open_deterministic_gzip_text_writer(path) as output_handle:
            yield output_handle
        return
    if compression_format == "bzip2":
        with bz2.open(
            path,
            "wt",
            encoding="utf-8",
            compresslevel=9,
        ) as output_handle:
            yield output_handle
        return
    if compression_format == "xz":
        with lzma.open(
            path,
            "wt",
            encoding="utf-8",
            preset=6,
        ) as output_handle:
            yield output_handle
        return
    raise ValueError("matrix_compression must be gzip, bzip2, or xz")


def dense_integer_tsv_to_sparse_bundle(
    source_path: Path,
    destination_path: Path,
    matrix_compression: str = "xz",
) -> None:
    """Convert a dense count table to a deterministic sparse bundle in two passes."""

    if matrix_compression not in MATRIX_COMPRESSION_SUFFIX_BY_FORMAT:
        raise ValueError("matrix_compression must be gzip, bzip2, or xz")
    source_path = source_path.resolve(strict=True)
    destination_path.mkdir(parents=True)
    gene_names: list[str] = []
    observed_gene_names: set[str] = set()
    library_sizes: np.ndarray | None = None
    nonzero_count = 0

    with (
        open_text_reader(source_path) as input_handle,
        open_deterministic_gzip_text_writer(
            destination_path / "genes.tsv.gz"
        ) as gene_names_handle,
        open_deterministic_gzip_text_writer(
            destination_path / "labels.tsv.gz"
        ) as cell_type_labels_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        header = next(reader)
        cell_type_labels = header[1:]
        if not cell_type_labels:
            raise ValueError("reference must contain at least one cell")
        cell_type_labels_handle.write("\n".join(cell_type_labels))
        cell_type_labels_handle.write("\n")
        library_sizes = np.zeros(len(cell_type_labels), dtype=np.int64)
        for row in reader:
            gene_name = row[0]
            if gene_name in observed_gene_names:
                raise ValueError(f"duplicate gene identifier: {gene_name}")
            observed_gene_names.add(gene_name)
            gene_names.append(gene_name)
            gene_names_handle.write(f"{gene_name}\n")
            count_values = parse_dense_count_row(row, len(cell_type_labels))
            library_sizes += count_values
            nonzero_count += int(np.count_nonzero(count_values))

    if not gene_names:
        raise ValueError("reference must contain at least one gene")
    if np.any(library_sizes <= 0):
        raise ValueError("each cell library size must be positive")
    with open_deterministic_gzip_text_writer(
        destination_path / "library_sizes.tsv.gz"
    ) as library_sizes_handle:
        library_sizes_handle.write("\n".join(map(str, library_sizes.tolist())))
        library_sizes_handle.write("\n")

    matrix_path = destination_path / (
        "matrix.mtx" + MATRIX_COMPRESSION_SUFFIX_BY_FORMAT[matrix_compression]
    )
    with (
        open_text_reader(source_path) as input_handle,
        open_compressed_matrix_writer(
            matrix_path,
            matrix_compression,
        ) as matrix_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        next(reader)
        matrix_handle.write("%%MatrixMarket matrix coordinate integer general\n")
        matrix_handle.write("% generated by BistreRoc\n")
        matrix_handle.write(
            f"{len(gene_names)} {len(cell_type_labels)} {nonzero_count}\n"
        )
        for row_index, row in enumerate(reader, start=1):
            if row[0] != gene_names[row_index - 1]:
                raise ValueError("reference gene order changed between read passes")
            count_values = parse_dense_count_row(row, len(cell_type_labels))
            nonzero_column_indices = np.flatnonzero(count_values)
            matrix_handle.write(
                "".join(
                    f"{row_index} {column_index + 1} {count_values[column_index]}\n"
                    for column_index in nonzero_column_indices
                )
            )


def sparse_bundle_bytes(reference_bundle_path: Path) -> int:
    """Return the stored bytes in a complete sparse reference bundle."""

    stored_paths = [
        reference_bundle_path / filename
        for filename in SPARSE_REFERENCE_METADATA_FILENAMES
    ]
    stored_paths.append(resolve_matrix_path(reference_bundle_path))
    return sum(path.stat().st_size for path in stored_paths)


def resolve_matrix_path(reference_bundle_path: Path) -> Path:
    """Resolve the single supported coordinate-matrix file in a bundle."""

    candidates = [
        reference_bundle_path / f"matrix.mtx{suffix}"
        for suffix in MATRIX_COMPRESSION_SUFFIX_BY_FORMAT.values()
    ]
    present_candidates = [path for path in candidates if path.is_file()]
    if len(present_candidates) != 1:
        raise ValueError("sparse bundle must contain exactly one compressed matrix")
    return present_candidates[0]


def read_gzip_text_vector(path: Path) -> list[str]:
    """Read a one-column gzip table as an ordered string vector."""

    with gzip.open(path, "rt", encoding="utf-8") as input_handle:
        return [line.rstrip("\r\n").split("\t")[0] for line in input_handle]


def read_library_sizes(path: Path, expected_count: int) -> np.ndarray:
    """Read one positive integer library size per reference cell."""

    with gzip.open(path, "rt", encoding="ascii") as input_handle:
        library_sizes = np.asarray(
            [line.rstrip("\r\n") for line in input_handle],
            dtype=np.int64,
        )
    if len(library_sizes) != expected_count:
        raise ValueError("library size count does not match reference columns")
    if np.any(library_sizes <= 0):
        raise ValueError("each cell library size must be positive")
    return library_sizes


def encode_cell_type_labels(
    cell_type_labels: list[str],
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Encode labels in first-seen order and return per-type cell counts."""

    cell_type_index_by_name: dict[str, int] = {}
    cell_type_indices_by_cell = np.empty(len(cell_type_labels), dtype=np.int64)
    for cell_index, cell_type_name in enumerate(cell_type_labels):
        if cell_type_name not in cell_type_index_by_name:
            cell_type_index_by_name[cell_type_name] = len(cell_type_index_by_name)
        cell_type_indices_by_cell[cell_index] = cell_type_index_by_name[cell_type_name]
    cell_counts_by_type = np.bincount(
        cell_type_indices_by_cell,
        minlength=len(cell_type_index_by_name),
    )
    return (
        list(cell_type_index_by_name),
        cell_type_indices_by_cell,
        cell_counts_by_type,
    )


def load_dense_reference_means(
    reference_path: Path,
) -> tuple[list[str], list[str], np.ndarray]:
    """Compute cell-type mean CPM from a dense count table in two passes."""

    reference_path = reference_path.resolve(strict=True)
    gene_names: list[str] = []
    observed_gene_names: set[str] = set()
    with open_text_reader(reference_path) as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = next(reader)
        cell_type_labels = header[1:]
        if not cell_type_labels:
            raise ValueError("reference must contain at least one cell")
        library_sizes = np.zeros(len(cell_type_labels), dtype=np.int64)
        for row in reader:
            gene_name = row[0]
            if gene_name in observed_gene_names:
                raise ValueError(f"duplicate gene identifier: {gene_name}")
            observed_gene_names.add(gene_name)
            gene_names.append(gene_name)
            library_sizes += parse_dense_count_row(row, len(cell_type_labels))
    if not gene_names:
        raise ValueError("reference must contain at least one gene")
    if np.any(library_sizes <= 0):
        raise ValueError("each cell library size must be positive")

    cell_type_names, cell_type_indices_by_cell, cell_counts_by_type = (
        encode_cell_type_labels(cell_type_labels)
    )
    mean_expression_by_cell_type = np.empty(
        (len(gene_names), len(cell_type_names)),
        dtype=np.float64,
    )
    with open_text_reader(reference_path) as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        next(reader)
        for row_index, row in enumerate(reader):
            if row[0] != gene_names[row_index]:
                raise ValueError("reference gene order changed between read passes")
            count_values = parse_dense_count_row(row, len(cell_type_labels))
            counts_per_million = count_values / library_sizes * 1_000_000.0
            mean_expression_by_cell_type[row_index] = (
                np.bincount(
                    cell_type_indices_by_cell,
                    weights=counts_per_million,
                    minlength=len(cell_type_names),
                )
                / cell_counts_by_type
            )
    return gene_names, cell_type_names, mean_expression_by_cell_type


def iter_matrix_market_data_lines(matrix_handle: TextIO) -> Iterator[str]:
    """Yield non-comment coordinate lines from an opened Matrix Market file."""

    for raw_line in matrix_handle:
        line = raw_line.strip()
        if line and not line.startswith("%"):
            yield line


def load_sparse_reference_means(
    reference_bundle_path: Path,
) -> tuple[list[str], list[str], np.ndarray]:
    """Compute cell-type mean CPM from a row-major sparse reference bundle."""

    reference_bundle_path = reference_bundle_path.resolve(strict=True)
    gene_names = read_gzip_text_vector(reference_bundle_path / "genes.tsv.gz")
    cell_type_labels = read_gzip_text_vector(reference_bundle_path / "labels.tsv.gz")
    if not gene_names or not cell_type_labels:
        raise ValueError("sparse reference annotations cannot be empty")
    if len(set(gene_names)) != len(gene_names):
        raise ValueError("sparse reference contains duplicate gene identifiers")
    library_sizes = read_library_sizes(
        reference_bundle_path / "library_sizes.tsv.gz",
        len(cell_type_labels),
    )
    cell_type_names, cell_type_indices_by_cell, cell_counts_by_type = (
        encode_cell_type_labels(cell_type_labels)
    )
    expression_sums_by_cell_type = np.zeros(
        (len(gene_names), len(cell_type_names)),
        dtype=np.float64,
    )
    observed_library_sizes = np.zeros(len(cell_type_labels), dtype=np.int64)

    with open_text_reader(resolve_matrix_path(reference_bundle_path)) as matrix_handle:
        banner = matrix_handle.readline().lower().split()
        if banner != [
            "%%matrixmarket",
            "matrix",
            "coordinate",
            "integer",
            "general",
        ]:
            raise ValueError(
                "signature generation requires an integer coordinate matrix"
            )
        dimensions_line = next(iter_matrix_market_data_lines(matrix_handle))
        row_count, column_count, declared_nonzero_count = map(
            int,
            dimensions_line.split(),
        )
        if row_count != len(gene_names) or column_count != len(cell_type_labels):
            raise ValueError("matrix dimensions do not match reference annotations")

        coordinate_lines = iter_matrix_market_data_lines(matrix_handle)
        processed_nonzero_count = 0
        previous_coordinate_key = -1
        while True:
            coordinate_line_chunk = list(
                islice(coordinate_lines, MATRIX_COORDINATE_CHUNK_SIZE)
            )
            if not coordinate_line_chunk:
                break
            coordinates = np.loadtxt(
                coordinate_line_chunk,
                dtype=np.int64,
                ndmin=2,
            )
            if coordinates.shape[1] != 3:
                raise ValueError(
                    "each sparse coordinate requires row, column, and value"
                )
            row_indices = coordinates[:, 0] - 1
            column_indices = coordinates[:, 1] - 1
            count_values = coordinates[:, 2]
            if (
                np.any(row_indices < 0)
                or np.any(row_indices >= row_count)
                or np.any(column_indices < 0)
                or np.any(column_indices >= column_count)
                or np.any(count_values < 0)
            ):
                raise ValueError("sparse reference coordinate is out of bounds")
            coordinate_keys = row_indices * column_count + column_indices
            if coordinate_keys[0] <= previous_coordinate_key or np.any(
                np.diff(coordinate_keys) <= 0
            ):
                raise ValueError("sparse coordinates must be unique and row-major")
            previous_coordinate_key = int(coordinate_keys[-1])
            np.add.at(observed_library_sizes, column_indices, count_values)
            counts_per_million = (
                count_values / library_sizes[column_indices] * 1_000_000.0
            )
            np.add.at(
                expression_sums_by_cell_type,
                (row_indices, cell_type_indices_by_cell[column_indices]),
                counts_per_million,
            )
            processed_nonzero_count += len(coordinates)

    if processed_nonzero_count != declared_nonzero_count:
        raise ValueError("sparse coordinate count does not match its declaration")
    if not np.array_equal(observed_library_sizes, library_sizes):
        raise ValueError("stored library sizes do not match matrix column sums")
    mean_expression_by_cell_type = expression_sums_by_cell_type / cell_counts_by_type
    return gene_names, cell_type_names, mean_expression_by_cell_type


def load_reference_means(
    reference_path: Path,
) -> tuple[list[str], list[str], np.ndarray]:
    """Load mean CPM by cell type from a dense table or sparse bundle."""

    if reference_path.is_dir():
        return load_sparse_reference_means(reference_path)
    return load_dense_reference_means(reference_path)

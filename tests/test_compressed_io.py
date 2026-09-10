"""Tests for deterministic compressed text input and streaming."""

import bz2
import gzip
import lzma
import os
import threading
from collections.abc import Callable
from pathlib import Path

import pytest

from bistreroc.compressed_io import (
    compression_format,
    open_text_reader,
    open_text_writer,
    stream_decompressed_file_to_named_pipe,
)


@pytest.mark.parametrize(
    ("suffix", "expected_format"),
    [
        (".tsv", "plain"),
        (".gz", "gzip"),
        (".bz2", "bzip2"),
        (".xz", "xz"),
    ],
)
def test_compression_format_uses_filename_suffix(
    suffix: str,
    expected_format: str,
) -> None:
    """Map each supported filename suffix to its compression format."""

    assert compression_format(Path(f"matrix{suffix}")) == expected_format


@pytest.mark.parametrize("suffix", ["", ".gz", ".bz2", ".xz"])
def test_text_readers_and_writers_round_trip(
    tmp_path: Path,
    suffix: str,
) -> None:
    """Round-trip UTF-8 text through every supported storage format."""

    output_path = tmp_path / f"matrix.tsv{suffix}"
    expected_text = "Gene\tmixture_1\nGène-1\t3.25\n"
    with open_text_writer(output_path) as output_handle:
        output_handle.write(expected_text)
    with open_text_reader(output_path) as input_handle:
        observed_text = input_handle.read()

    assert observed_text == expected_text


def test_gzip_writer_is_byte_deterministic(tmp_path: Path) -> None:
    """Exclude timestamps and filenames from gzip output metadata."""

    first_path = tmp_path / "first.tsv.gz"
    second_path = tmp_path / "second.tsv.gz"
    for output_path in (first_path, second_path):
        with open_text_writer(output_path) as output_handle:
            output_handle.write("Gene\tsample\nG1\t1\n")

    assert first_path.read_bytes() == second_path.read_bytes()


@pytest.mark.parametrize(
    ("suffix", "compress"),
    [
        (".gz", gzip.compress),
        (".bz2", bz2.compress),
        (".xz", lzma.compress),
    ],
)
def test_stream_decompressed_file_to_named_pipe(
    tmp_path: Path,
    suffix: str,
    compress: Callable[[bytes], bytes],
) -> None:
    """Recover compressed bytes through a pipe without a staged regular file."""

    expected_payload = (b"Gene\tA\tB\nGENE1\t1\t2\n" * 10_000) + b"END\n"
    compressed_path = tmp_path / f"matrix.tsv{suffix}"
    compressed_path.write_bytes(compress(expected_payload))
    named_pipe_path = tmp_path / "matrix.tsv"
    os.mkfifo(named_pipe_path)
    stop_event = threading.Event()
    error_buffer: list[Exception] = []
    decompression_thread = threading.Thread(
        target=stream_decompressed_file_to_named_pipe,
        args=(compressed_path, named_pipe_path, stop_event, error_buffer),
    )
    decompression_thread.start()
    observed_payload = named_pipe_path.read_bytes()
    decompression_thread.join(timeout=5)

    assert observed_payload == expected_payload
    assert not decompression_thread.is_alive()
    assert error_buffer == []
    assert not named_pipe_path.is_file()

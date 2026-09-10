"""Bounded-memory streaming decompression for BistreRoc inputs."""

import bz2
import errno
import gzip
import io
import lzma
import os
import shutil
import threading
import zlib
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import BinaryIO, TextIO

COPY_BUFFER_BYTES = 1024 * 1024


def compression_format(path: Path) -> str:
    """Return the supported compression format encoded by a path suffix."""

    if path.suffix == ".gz":
        return "gzip"
    if path.suffix == ".bz2":
        return "bzip2"
    if path.suffix == ".xz":
        return "xz"
    return "plain"


def compressed_binary_opener(path: Path) -> Callable[..., BinaryIO]:
    """Return the binary decompressor selected from a compressed path suffix."""

    file_compression = compression_format(path)
    if file_compression == "gzip":
        return gzip.open
    if file_compression == "bzip2":
        return bz2.open
    if file_compression == "xz":
        return lzma.open
    raise ValueError(f"No decompressor for uncompressed input: {path}")


@contextmanager
def open_text_reader(path: Path) -> Iterator[TextIO]:
    """Open a plain or suffix-compressed UTF-8 text file for reading."""

    if compression_format(path) == "plain":
        with path.open("r", encoding="utf-8", newline="") as input_handle:
            yield input_handle
        return
    with compressed_binary_opener(path)(
        path,
        "rt",
        encoding="utf-8",
        newline="",
    ) as input_handle:
        yield input_handle


@contextmanager
def open_text_writer(path: Path) -> Iterator[TextIO]:
    """Open a plain or suffix-compressed UTF-8 text file for writing."""

    file_compression = compression_format(path)
    if file_compression == "plain":
        with path.open("w", encoding="utf-8", newline="") as output_handle:
            yield output_handle
        return
    if file_compression == "gzip":
        with (
            path.open("wb") as raw_handle,
            gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                mtime=0,
            ) as gzip_handle,
            io.TextIOWrapper(
                gzip_handle, encoding="utf-8", newline=""
            ) as output_handle,
        ):
            yield output_handle
        return
    opener = bz2.open if file_compression == "bzip2" else lzma.open
    with opener(path, "wt", encoding="utf-8", newline="") as output_handle:
        yield output_handle


def transcode_text_file(source_path: Path, destination_path: Path) -> None:
    """Copy text while decoding and encoding suffix-selected compression."""

    with (
        open_text_reader(source_path) as input_handle,
        open_text_writer(destination_path) as output_handle,
    ):
        shutil.copyfileobj(input_handle, output_handle, length=COPY_BUFFER_BYTES)


def open_named_pipe_writer(
    named_pipe_path: Path,
    stop_event: threading.Event,
    error_buffer: list[Exception],
) -> int | None:
    """Open a FIFO writer after its consumer appears, or stop cleanly."""

    pipe_descriptor = -1
    while pipe_descriptor == -1 and not stop_event.is_set():
        try:
            pipe_descriptor = os.open(
                named_pipe_path,
                os.O_WRONLY | os.O_NONBLOCK,
            )
        except OSError as error:
            if error.errno != errno.ENXIO:
                error_buffer.append(error)
                return None
            stop_event.wait(0.05)
    if pipe_descriptor == -1:
        return None
    os.set_blocking(pipe_descriptor, True)
    return pipe_descriptor


def stream_decompressed_file_to_named_pipe(
    source_path: Path,
    named_pipe_path: Path,
    stop_event: threading.Event,
    error_buffer: list[Exception],
) -> None:
    """Stream one compressed input into a FIFO with a bounded copy buffer."""

    pipe_descriptor = open_named_pipe_writer(
        named_pipe_path,
        stop_event,
        error_buffer,
    )
    if pipe_descriptor is None:
        return

    try:
        with (
            compressed_binary_opener(source_path)(source_path, "rb") as input_handle,
            os.fdopen(pipe_descriptor, "wb", buffering=0) as output_handle,
        ):
            while not stop_event.is_set():
                chunk = input_handle.read(COPY_BUFFER_BYTES)
                if not chunk:
                    return
                output_handle.write(chunk)
    except (EOFError, OSError, lzma.LZMAError, zlib.error) as error:
        error_buffer.append(error)

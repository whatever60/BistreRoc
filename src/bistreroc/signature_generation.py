"""Signature-matrix generation from labeled single-cell count matrices."""

import os
import signal
import subprocess
import tempfile
import threading
from importlib.resources import files
from pathlib import Path
from typing import TextIO

import numpy as np

from bistreroc.compressed_io import open_named_pipe_writer, transcode_text_file
from bistreroc.sparse_reference import load_reference_means

CELL_TYPE_PROFILE_REPLICATE_COUNT = 5


def signature_selection_script_path() -> Path:
    """Return the installed R script that selects signature genes."""

    return Path(
        str(files("bistreroc").joinpath("r/select_signature_genes.R"))
    ).resolve()


def write_cell_type_profiles(
    output_handle: TextIO,
    gene_names: list[str],
    cell_type_names: list[str],
    mean_expression_by_cell_type: np.ndarray,
) -> None:
    """Write a cell-type mean CPM table for the statistical selection kernel."""

    output_handle.write("Gene\t")
    replicated_profile_names = [
        cell_type if replicate_index == 0 else f"{cell_type}.{replicate_index}"
        for cell_type in cell_type_names
        for replicate_index in range(CELL_TYPE_PROFILE_REPLICATE_COUNT)
    ]
    output_handle.write("\t".join(replicated_profile_names))
    output_handle.write("\n")
    for gene_name, cell_type_expression in zip(
        gene_names,
        mean_expression_by_cell_type,
        strict=True,
    ):
        formatted_values = "\t".join(
            format(float(value), ".17f")
            for value in cell_type_expression
            for _ in range(CELL_TYPE_PROFILE_REPLICATE_COUNT)
        )
        output_handle.write(f"{gene_name}\t{formatted_values}\n")


def stream_cell_type_profiles_to_named_pipe(
    named_pipe_path: Path,
    gene_names: list[str],
    cell_type_names: list[str],
    mean_expression_by_cell_type: np.ndarray,
    stop_event: threading.Event,
    error_buffer: list[Exception],
) -> None:
    """Stream cell-type profiles into a FIFO and propagate writer failures."""

    pipe_descriptor = open_named_pipe_writer(
        named_pipe_path,
        stop_event,
        error_buffer,
    )
    if pipe_descriptor is None:
        return
    try:
        with os.fdopen(
            pipe_descriptor,
            "w",
            encoding="utf-8",
            newline="",
        ) as output_handle:
            write_cell_type_profiles(
                output_handle,
                gene_names,
                cell_type_names,
                mean_expression_by_cell_type,
            )
    except Exception as error:
        error_buffer.append(error)


def build_signature(
    reference_path: Path,
    output_path: Path,
    minimum_markers_per_cell_type: int = 300,
    maximum_markers_per_cell_type: int = 500,
    q_value_threshold: float = 0.01,
    core_count: int = 1,
    timeout_seconds: int = 3600,
) -> Path:
    """Build a discriminative signature matrix from labeled cell counts."""

    if minimum_markers_per_cell_type < 1:
        raise ValueError("minimum_markers_per_cell_type must be positive")
    if maximum_markers_per_cell_type < minimum_markers_per_cell_type:
        raise ValueError(
            "maximum_markers_per_cell_type cannot be smaller than "
            "minimum_markers_per_cell_type"
        )
    if not 0 <= q_value_threshold <= 1:
        raise ValueError("q_value_threshold must be between zero and one")
    if core_count < 1:
        raise ValueError("core_count must be a positive integer")
    reference_path = reference_path.resolve(strict=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path = output_path.resolve()
    gene_names, cell_type_names, mean_expression_by_cell_type = load_reference_means(
        reference_path
    )

    with tempfile.TemporaryDirectory(
        prefix="bistreroc_signature_"
    ) as temporary_directory_name:
        temporary_directory = Path(temporary_directory_name)
        profile_path = temporary_directory / "cell_type_profiles.tsv"
        os.mkfifo(profile_path)
        profile_stop_event = threading.Event()
        profile_errors: list[Exception] = []
        profile_thread = threading.Thread(
            target=stream_cell_type_profiles_to_named_pipe,
            args=(
                profile_path,
                gene_names,
                cell_type_names,
                mean_expression_by_cell_type,
                profile_stop_event,
                profile_errors,
            ),
            daemon=True,
        )
        selection_output_path = temporary_directory / "signature.tsv"
        rscript_command = [
            "Rscript",
            str(signature_selection_script_path()),
            str(profile_path),
            str(selection_output_path),
            str(minimum_markers_per_cell_type),
            str(maximum_markers_per_cell_type),
            str(q_value_threshold),
            str(core_count),
        ]
        profile_thread.start()
        try:
            rscript_process = subprocess.Popen(
                rscript_command,
                start_new_session=True,
            )
            try:
                rscript_return_code = rscript_process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                os.killpg(rscript_process.pid, signal.SIGKILL)
                rscript_process.wait()
                raise RuntimeError(
                    f"signature generation timed out after {timeout_seconds} seconds"
                ) from None
        finally:
            profile_stop_event.set()
            profile_thread.join(timeout=5)
        if profile_errors:
            raise RuntimeError("profile streaming failed") from profile_errors[0]
        if profile_thread.is_alive():
            raise RuntimeError("profile streaming did not terminate")
        if rscript_return_code != 0:
            raise RuntimeError(
                f"signature generation exited with code {rscript_return_code}"
            )
        transcode_text_file(selection_output_path, output_path)

    return output_path

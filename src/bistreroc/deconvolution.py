"""Fraction estimation with a standalone linear nu-SVR implementation."""

import os
import signal
import subprocess
import tempfile
import threading
from importlib.resources import files
from pathlib import Path

from bistreroc.compressed_io import (
    compression_format,
    stream_decompressed_file_to_named_pipe,
    transcode_text_file,
)


def fraction_estimation_script_path() -> Path:
    """Return the installed R script that performs fraction estimation."""

    return Path(str(files("bistreroc").joinpath("r/estimate_fractions.R"))).resolve()


def estimate_fractions(
    mixture_path: Path,
    signature_path: Path,
    output_path: Path,
    core_count: int = 1,
    timeout_seconds: int = 3600,
) -> Path:
    """Estimate cell fractions and write a diagnostic-rich TSV table."""

    if core_count < 1:
        raise ValueError("core_count must be a positive integer")
    mixture_path = mixture_path.resolve(strict=True)
    signature_path = signature_path.resolve(strict=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path = output_path.resolve()

    decompression_stop_event = threading.Event()
    decompression_errors: list[Exception] = []
    decompression_threads: list[threading.Thread] = []

    with tempfile.TemporaryDirectory(
        prefix="bistreroc_deconvolution_"
    ) as temporary_directory_name:
        temporary_directory = Path(temporary_directory_name)
        staged_input_paths: dict[str, Path] = {}
        for input_name, source_path in (
            ("mixture", mixture_path),
            ("signature", signature_path),
        ):
            if compression_format(source_path) == "plain":
                staged_input_paths[input_name] = source_path
                continue
            named_pipe_path = temporary_directory / f"{input_name}.tsv"
            os.mkfifo(named_pipe_path)
            decompression_thread = threading.Thread(
                target=stream_decompressed_file_to_named_pipe,
                args=(
                    source_path,
                    named_pipe_path,
                    decompression_stop_event,
                    decompression_errors,
                ),
                daemon=True,
            )
            decompression_threads.append(decompression_thread)
            staged_input_paths[input_name] = named_pipe_path

        solver_output_path = temporary_directory / "fractions.tsv"
        rscript_command = [
            "Rscript",
            str(fraction_estimation_script_path()),
            str(staged_input_paths["signature"]),
            str(staged_input_paths["mixture"]),
            str(solver_output_path),
            str(core_count),
        ]

        for decompression_thread in decompression_threads:
            decompression_thread.start()
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
                    f"fraction estimation timed out after {timeout_seconds} seconds"
                ) from None
        finally:
            decompression_stop_event.set()
            for decompression_thread in decompression_threads:
                decompression_thread.join(timeout=5)

        if decompression_errors:
            raise RuntimeError("input decompression failed") from decompression_errors[
                0
            ]
        if any(thread.is_alive() for thread in decompression_threads):
            raise RuntimeError("input decompression did not terminate")
        if rscript_return_code != 0:
            raise RuntimeError(
                f"fraction estimation exited with code {rscript_return_code}"
            )
        transcode_text_file(solver_output_path, output_path)

    return output_path

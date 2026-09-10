"""Benchmark signature generation and fraction estimation at realistic scale."""

import argparse
import bz2
import csv
import gzip
import hashlib
import json
import lzma
import os
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import TextIO

import numpy as np

from bistreroc.sparse_reference import (
    dense_integer_tsv_to_sparse_bundle,
    resolve_matrix_path,
    sparse_bundle_bytes,
)

DEFAULT_ROOT = Path("benchmark_results/realistic_8_mixtures")
CANONICAL_MARKER_COUNT = 14
CANONICAL_Q_VALUE_THRESHOLD = 0.01
REFERENCE_IMAGE = (
    "docker.io/cibersortx/fractions@"
    "sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc"
)
IMPLEMENTATIONS = ("bistreroc", "reference")
STAGES = ("signature_generation", "fraction_estimation", "end_to_end")


def utc_timestamp() -> str:
    """Return the current UTC time in a stable machine-readable form."""

    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one file without loading it into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as input_handle:
        for chunk in iter(lambda: input_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def directory_bytes(path: Path) -> int:
    """Return the retained size of all regular files below a directory."""

    return sum(
        candidate.stat().st_size for candidate in path.rglob("*") if candidate.is_file()
    )


def host_metadata() -> dict[str, int | str]:
    """Describe the operating system and CPU allocation of this process."""

    cpu_models = [
        line.split(":", 1)[1].strip()
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines()
        if line.startswith("model name")
    ]
    return {
        "system": platform.system(),
        "kernel": platform.release(),
        "machine": platform.machine(),
        "cpu_model": cpu_models[0],
        "logical_cpus": os.cpu_count() or 1,
        "available_cpus": len(os.sched_getaffinity(0)),
    }


def r_runtime_metadata() -> dict[str, str]:
    """Return versions for the R runtime and required statistical packages."""

    version_query = (
        "cat(as.character(getRversion()), "
        'as.character(packageVersion("data.table")), '
        'as.character(packageVersion("e1071")), sep="\\t")'
    )
    completed_query = subprocess.run(
        ["Rscript", "-e", version_query],
        check=True,
        capture_output=True,
        text=True,
    )
    r_version, data_table_version, e1071_version = completed_query.stdout.strip().split(
        "\t"
    )
    return {
        "r": r_version,
        "data.table": data_table_version,
        "e1071": e1071_version,
    }


def open_numeric_table(path: Path) -> TextIO:
    """Open a plain or compressed numeric table for reading."""

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    if path.suffix == ".bz2":
        return bz2.open(path, "rt", encoding="utf-8", newline="")
    if path.suffix in {".xz", ".lzma"}:
        return lzma.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def numeric_table(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Read row labels, numeric-column labels, and values from a TSV table."""

    with open_numeric_table(path) as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    if len(header) < 2 or not rows:
        raise ValueError(f"numeric table is empty: {path}")
    row_names = [row[0] for row in rows]
    values = np.asarray([row[1:] for row in rows], dtype=np.float64)
    if values.shape != (len(rows), len(header) - 1):
        raise ValueError(f"numeric table has inconsistent row widths: {path}")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"numeric table contains non-finite values: {path}")
    return row_names, header[1:], values


def maximum_numeric_difference(
    reference_path: Path,
    bistreroc_path: Path,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> float:
    """Require compatible tables and return their maximum absolute difference."""

    reference_rows, reference_columns, reference_values = numeric_table(reference_path)
    bistreroc_rows, bistreroc_columns, bistreroc_values = numeric_table(bistreroc_path)
    if reference_rows != bistreroc_rows:
        raise AssertionError("numeric result row labels differ")
    if reference_columns != bistreroc_columns:
        raise AssertionError("numeric result column labels differ")
    np.testing.assert_allclose(
        bistreroc_values,
        reference_values,
        rtol=relative_tolerance,
        atol=absolute_tolerance,
    )
    return float(np.max(np.abs(bistreroc_values - reference_values)))


def write_mixture_subset(
    source: Path, dense_output: Path, compressed_output: Path, mixture_count: int
) -> int:
    """Write the first requested mixture columns as dense TSV and XZ TSV."""

    written_rows = 0
    with (
        source.open("r", encoding="utf-8", newline="") as input_handle,
        dense_output.open("w", encoding="utf-8", newline="") as dense_handle,
        lzma.open(
            compressed_output, "wt", encoding="utf-8", newline="", preset=6
        ) as compressed_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        dense_writer = csv.writer(dense_handle, delimiter="\t", lineterminator="\n")
        compressed_writer = csv.writer(
            compressed_handle,
            delimiter="\t",
            lineterminator="\n",
        )
        header = next(reader)
        if len(header) <= mixture_count:
            raise ValueError(
                f"mixture provides {len(header) - 1} samples, fewer than {mixture_count}"
            )
        selected_width = mixture_count + 1
        dense_writer.writerow(header[:selected_width])
        compressed_writer.writerow(header[:selected_width])
        for row in reader:
            if len(row) != len(header):
                raise ValueError("mixture row width does not match its header")
            selected_row = row[:selected_width]
            dense_writer.writerow(selected_row)
            compressed_writer.writerow(selected_row)
            written_rows += 1
    if written_rows == 0:
        raise ValueError("mixture contains no genes")
    return written_rows


def sparse_dimensions(bundle: Path) -> tuple[int, int, int]:
    """Read matrix dimensions from a compressed sparse reference bundle."""

    matrix_path = resolve_matrix_path(bundle)
    with open_numeric_table(matrix_path) as matrix_handle:
        banner = matrix_handle.readline().strip().lower()
        if banner != "%%matrixmarket matrix coordinate integer general":
            raise ValueError("unexpected sparse reference Matrix Market banner")
        dimension_line = next(
            line for line in matrix_handle if line.strip() and not line.startswith("%")
        )
    return tuple(map(int, dimension_line.split()))


def prepare_inputs(
    root: Path,
    reference_source: Path,
    mixture_source: Path,
    mixture_count: int,
) -> Path:
    """Create paired input representations for the canonical benchmark."""

    if mixture_count < 1:
        raise ValueError("mixture_count must be positive")
    if root.exists():
        raise FileExistsError(f"benchmark root already exists: {root}")
    reference_source = reference_source.resolve(strict=True)
    mixture_source = mixture_source.resolve(strict=True)
    reference_input = root / "inputs" / "reference"
    bistreroc_input = root / "inputs" / "bistreroc"
    reference_input.mkdir(parents=True)
    bistreroc_input.mkdir(parents=True)

    dense_reference = reference_input / "reference.tsv"
    dense_reference.symlink_to(reference_source)
    sparse_reference = bistreroc_input / "reference.bundle"
    dense_integer_tsv_to_sparse_bundle(
        reference_source,
        sparse_reference,
        matrix_compression="xz",
    )
    dense_mixture = reference_input / "mixture.tsv"
    compressed_mixture = bistreroc_input / "mixture.tsv.xz"
    mixture_gene_count = write_mixture_subset(
        mixture_source,
        dense_mixture,
        compressed_mixture,
        mixture_count,
    )
    reference_gene_count, reference_cell_count, nonzero_count = sparse_dimensions(
        sparse_reference
    )
    with gzip.open(
        sparse_reference / "labels.tsv.gz",
        "rt",
        encoding="utf-8",
    ) as label_handle:
        labels = [line.rstrip("\r\n") for line in label_handle]

    fixture = {
        "schema_version": 1,
        "prepared_at_utc": utc_timestamp(),
        "scenario": {
            "name": (f"realistic_sparse_reference_{mixture_count}_mixtures_g14"),
            "reference_gene_count": reference_gene_count,
            "reference_cell_count": reference_cell_count,
            "cell_type_count": len(set(labels)),
            "reference_nonzero_count": nonzero_count,
            "reference_density": nonzero_count
            / (reference_gene_count * reference_cell_count),
            "mixture_gene_count": mixture_gene_count,
            "mixture_count": mixture_count,
            "markers_per_cell_type": CANONICAL_MARKER_COUNT,
            "q_value_threshold": CANONICAL_Q_VALUE_THRESHOLD,
            "expected_signature_scale": "approximately 500 unique genes",
        },
        "sources": {
            "reference": {
                "filename": reference_source.name,
                "bytes": reference_source.stat().st_size,
                "sha256": sha256_file(reference_source),
            },
            "mixture": {
                "filename": mixture_source.name,
                "bytes": mixture_source.stat().st_size,
                "sha256": sha256_file(mixture_source),
            },
        },
        "representations": {
            "reference": {
                "reference_bytes": dense_reference.stat().st_size,
                "bistreroc_bytes": sparse_bundle_bytes(sparse_reference),
            },
            "mixture": {
                "reference_bytes": dense_mixture.stat().st_size,
                "bistreroc_bytes": compressed_mixture.stat().st_size,
            },
        },
        "timing_boundary": (
            "input conversion is preparation and is excluded from timed stages"
        ),
    }
    fixture_path = root / "fixture.json"
    fixture_path.write_text(f"{json.dumps(fixture, indent=2)}\n", encoding="utf-8")
    return fixture_path


def process_snapshot() -> dict[int, tuple[int, int]]:
    """Return parent PID and resident bytes for each readable Linux process."""

    page_size = os.sysconf("SC_PAGE_SIZE")
    snapshot: dict[int, tuple[int, int]] = {}
    for process_directory in Path("/proc").glob("[0-9]*"):
        try:
            stat_fields = (
                (process_directory / "stat")
                .read_text(encoding="utf-8")
                .rsplit(")", 1)[1]
                .split()
            )
            resident_pages = int(
                (process_directory / "statm").read_text(encoding="utf-8").split()[1]
            )
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        snapshot[int(process_directory.name)] = (
            int(stat_fields[1]),
            resident_pages * page_size,
        )
    return snapshot


def process_tree_resident_bytes(root_pids: set[int]) -> int:
    """Return aggregate RSS for roots and all descendants at one instant."""

    snapshot = process_snapshot()
    selected_pids = set(root_pids)
    previous_count = -1
    while len(selected_pids) != previous_count:
        previous_count = len(selected_pids)
        selected_pids.update(
            pid
            for pid, (parent_pid, _) in snapshot.items()
            if parent_pid in selected_pids
        )
    return sum(snapshot[pid][1] for pid in selected_pids if pid in snapshot)


def running_container_pid(cidfile: Path) -> int | None:
    """Return the running container PID once Podman has populated its cidfile."""

    if not cidfile.is_file():
        return None
    container_id = cidfile.read_text(encoding="utf-8").strip()
    if not container_id:
        return None
    inspection = subprocess.run(
        ["podman", "inspect", "--format", "{{.State.Pid}}", container_id],
        check=False,
        capture_output=True,
        text=True,
    )
    if inspection.returncode != 0 or not inspection.stdout.strip():
        return None
    container_pid = int(inspection.stdout.strip())
    return container_pid or None


def measure_command(
    command: list[str],
    stdout_path: Path,
    stderr_path: Path,
    sample_interval_seconds: float,
    cidfile: Path | None = None,
) -> dict[str, float | int]:
    """Run a command and measure wall time plus aggregate process-tree peak RSS."""

    started = time.perf_counter()
    peak_memory_bytes = 0
    container_pid: int | None = None
    with (
        stdout_path.open("wb") as stdout_handle,
        stderr_path.open("wb") as stderr_handle,
    ):
        process = subprocess.Popen(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
        )
        while process.poll() is None:
            if cidfile is not None and container_pid is None:
                container_pid = running_container_pid(cidfile)
            root_pids = {process.pid}
            if container_pid is not None:
                root_pids.add(container_pid)
            peak_memory_bytes = max(
                peak_memory_bytes,
                process_tree_resident_bytes(root_pids),
            )
            time.sleep(sample_interval_seconds)
        return_code = process.wait()
    wall_seconds = time.perf_counter() - started
    if return_code != 0:
        raise subprocess.CalledProcessError(
            return_code,
            command[:1],
            stderr=f"see {stderr_path}",
        )
    if peak_memory_bytes == 0:
        raise RuntimeError("process completed before memory could be sampled")
    return {
        "wall_seconds": wall_seconds,
        "peak_memory_bytes": peak_memory_bytes,
    }


def reference_signature_path(output_directory: Path) -> Path:
    """Return the single optimized signature emitted by the reference image."""

    candidates = list(output_directory.glob("*.bm.K*.txt"))
    if len(candidates) != 1:
        raise ValueError(
            "reference signature output must contain exactly one *.bm.K*.txt file"
        )
    return candidates[0]


def reference_fraction_path(output_directory: Path) -> Path:
    """Return the fraction table emitted by the reference image."""

    candidates = list(output_directory.glob("*Results*.txt"))
    if len(candidates) != 1:
        raise ValueError(
            "reference fraction output must contain exactly one *Results*.txt file"
        )
    return candidates[0]


def podman_command(
    cidfile: Path,
    mounts: list[tuple[Path, str, bool]],
    username: str,
    token: str,
    algorithm_arguments: list[str],
) -> list[str]:
    """Build a pinned reference-image command without persisting credentials."""

    command = [
        "podman",
        "run",
        "--rm",
        "--runtime=runc",
        "--network=slirp4netns",
        "--cidfile",
        str(cidfile),
    ]
    for host_path, container_path, read_only in mounts:
        mount_value = f"{host_path.resolve()}:{container_path}"
        if read_only:
            mount_value += ":ro"
        command.extend(["-v", mount_value])
    command.extend(
        [
            REFERENCE_IMAGE,
            "--username",
            username,
            "--token",
            token,
            *algorithm_arguments,
        ]
    )
    return command


def stage_record(
    root: Path,
    measurement: dict[str, float | int],
    input_bytes: int,
    output_directory: Path,
    primary_output: Path,
) -> dict[str, float | int | str]:
    """Complete a timed measurement with retained input and output storage."""

    row_names, numeric_columns, _ = numeric_table(primary_output)
    return {
        **measurement,
        "input_bytes": input_bytes,
        "output_bytes": directory_bytes(output_directory),
        "primary_output": str(primary_output.resolve().relative_to(root)),
        "primary_output_bytes": primary_output.stat().st_size,
        "primary_output_sha256": sha256_file(primary_output),
        "primary_output_rows": len(row_names),
        "primary_output_numeric_columns": len(numeric_columns),
    }


def run_bistreroc(
    root: Path,
    python_executable: Path,
    cores: int,
    marker_count: int,
    q_value_threshold: float,
    sample_interval_seconds: float,
) -> dict:
    """Measure the two BistreRoc stages as one sequential pipeline."""

    output_root = root / "runs" / "bistreroc"
    if output_root.exists():
        raise FileExistsError(f"BistreRoc measurement already exists: {output_root}")
    signature_directory = output_root / "signature_generation"
    fractions_directory = output_root / "fraction_estimation"
    log_directory = output_root / "logs"
    signature_directory.mkdir(parents=True)
    fractions_directory.mkdir(parents=True)
    log_directory.mkdir(parents=True)
    sparse_reference = root / "inputs" / "bistreroc" / "reference.bundle"
    mixture = root / "inputs" / "bistreroc" / "mixture.tsv.xz"
    signature = signature_directory / "signature.tsv.xz"
    fractions = fractions_directory / "fractions.tsv.xz"
    common_command = [str(python_executable), "-m", "bistreroc.cli"]
    signature_command = [
        *common_command,
        "build-signature",
        "--reference",
        str(sparse_reference),
        "--output",
        str(signature),
        "--minimum-markers-per-cell-type",
        str(marker_count),
        "--maximum-markers-per-cell-type",
        str(marker_count),
        "--q-value-threshold",
        str(q_value_threshold),
        "--cores",
        str(cores),
    ]
    signature_measurement = measure_command(
        signature_command,
        log_directory / "signature.stdout.log",
        log_directory / "signature.stderr.log",
        sample_interval_seconds,
    )
    signature_record = stage_record(
        root,
        signature_measurement,
        sparse_bundle_bytes(sparse_reference),
        signature_directory,
        signature,
    )
    fractions_command = [
        *common_command,
        "deconvolve",
        "--mixture",
        str(mixture),
        "--signature",
        str(signature),
        "--output",
        str(fractions),
        "--cores",
        str(cores),
    ]
    fractions_measurement = measure_command(
        fractions_command,
        log_directory / "fractions.stdout.log",
        log_directory / "fractions.stderr.log",
        sample_interval_seconds,
    )
    fractions_record = stage_record(
        root,
        fractions_measurement,
        mixture.stat().st_size + signature.stat().st_size,
        fractions_directory,
        fractions,
    )
    end_to_end_record = {
        "wall_seconds": (
            signature_record["wall_seconds"] + fractions_record["wall_seconds"]
        ),
        "peak_memory_bytes": max(
            signature_record["peak_memory_bytes"],
            fractions_record["peak_memory_bytes"],
        ),
        "input_bytes": sparse_bundle_bytes(sparse_reference) + mixture.stat().st_size,
        "output_bytes": (
            signature_record["output_bytes"] + fractions_record["output_bytes"]
        ),
        "derivation": "sum of sequential wall time; maximum stage peak memory",
    }
    return {
        "schema_version": 1,
        "implementation": "bistreroc",
        "recorded_at_utc": utc_timestamp(),
        "software": {
            "name": "BistreRoc",
            "version": version("bistreroc"),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "statistical_runtime": r_runtime_metadata(),
        },
        "configuration": {
            "cores": cores,
            "markers_per_cell_type": marker_count,
            "q_value_threshold": q_value_threshold,
            "memory_sampling_interval_seconds": sample_interval_seconds,
            "arguments_without_environment_paths": {
                "signature_generation": [
                    "build-signature",
                    "--reference",
                    "<sparse-reference-bundle>",
                    "--output",
                    "<compressed-signature>",
                    "--minimum-markers-per-cell-type",
                    str(marker_count),
                    "--maximum-markers-per-cell-type",
                    str(marker_count),
                    "--q-value-threshold",
                    str(q_value_threshold),
                    "--cores",
                    str(cores),
                ],
                "fraction_estimation": [
                    "deconvolve",
                    "--mixture",
                    "<compressed-mixture>",
                    "--signature",
                    "<compressed-signature>",
                    "--output",
                    "<compressed-fractions>",
                    "--cores",
                    str(cores),
                ],
            },
        },
        "host": host_metadata(),
        "stages": {
            "signature_generation": signature_record,
            "fraction_estimation": fractions_record,
            "end_to_end": end_to_end_record,
        },
    }


def run_reference(
    root: Path,
    username: str,
    token: str,
    marker_count: int,
    sample_interval_seconds: float,
) -> dict:
    """Measure the two stages of the pinned CIBERSORTx reference image."""

    output_root = root / "runs" / "reference"
    if output_root.exists():
        raise FileExistsError(f"reference measurement already exists: {output_root}")
    signature_directory = output_root / "signature_generation"
    fractions_directory = output_root / "fraction_estimation"
    log_directory = output_root / "logs"
    signature_directory.mkdir(parents=True)
    fractions_directory.mkdir(parents=True)
    log_directory.mkdir(parents=True)
    reference = root / "inputs" / "reference" / "reference.tsv"
    mixture = root / "inputs" / "reference" / "mixture.tsv"

    signature_arguments = [
        "--single_cell",
        "TRUE",
        "--refsample",
        "reference.tsv",
        "--sampling",
        "1",
        "--fraction",
        "0",
        "--G.min",
        str(marker_count),
        "--G.max",
        str(marker_count),
    ]
    signature_cidfile = log_directory / "signature.cid"
    signature_command = podman_command(
        signature_cidfile,
        [
            (reference, "/src/data/reference.tsv", True),
            (signature_directory, "/src/outdir", False),
        ],
        username,
        token,
        signature_arguments,
    )
    signature_measurement = measure_command(
        signature_command,
        log_directory / "signature.stdout.log",
        log_directory / "signature.stderr.log",
        sample_interval_seconds,
        signature_cidfile,
    )
    signature = reference_signature_path(signature_directory)
    signature_record = stage_record(
        root,
        signature_measurement,
        reference.stat().st_size,
        signature_directory,
        signature,
    )

    fractions_arguments = [
        "--mixture",
        "mixture.tsv",
        "--sigmatrix",
        "signature.tsv",
    ]
    fractions_cidfile = log_directory / "fractions.cid"
    fractions_command = podman_command(
        fractions_cidfile,
        [
            (mixture, "/src/data/mixture.tsv", True),
            (signature, "/src/data/signature.tsv", True),
            (fractions_directory, "/src/outdir", False),
        ],
        username,
        token,
        fractions_arguments,
    )
    fractions_measurement = measure_command(
        fractions_command,
        log_directory / "fractions.stdout.log",
        log_directory / "fractions.stderr.log",
        sample_interval_seconds,
        fractions_cidfile,
    )
    fractions = reference_fraction_path(fractions_directory)
    fractions_record = stage_record(
        root,
        fractions_measurement,
        mixture.stat().st_size + signature.stat().st_size,
        fractions_directory,
        fractions,
    )
    end_to_end_record = {
        "wall_seconds": (
            signature_record["wall_seconds"] + fractions_record["wall_seconds"]
        ),
        "peak_memory_bytes": max(
            signature_record["peak_memory_bytes"],
            fractions_record["peak_memory_bytes"],
        ),
        "input_bytes": reference.stat().st_size + mixture.stat().st_size,
        "output_bytes": (
            signature_record["output_bytes"] + fractions_record["output_bytes"]
        ),
        "derivation": "sum of sequential wall time; maximum stage peak memory",
    }
    return {
        "schema_version": 1,
        "implementation": "reference",
        "recorded_at_utc": utc_timestamp(),
        "software": {
            "name": "CIBERSORTx",
            "image": REFERENCE_IMAGE,
            "statistical_runtime": {
                "r": "3.6.3",
                "data.table": "1.12.8",
                "e1071": "1.7.3",
            },
        },
        "configuration": {
            "worker_limit": 3,
            "markers_per_cell_type": marker_count,
            "q_value_threshold": CANONICAL_Q_VALUE_THRESHOLD,
            "memory_sampling_interval_seconds": sample_interval_seconds,
            "arguments_without_credentials": {
                "signature_generation": signature_arguments,
                "fraction_estimation": fractions_arguments,
            },
        },
        "host": host_metadata(),
        "stages": {
            "signature_generation": signature_record,
            "fraction_estimation": fractions_record,
            "end_to_end": end_to_end_record,
        },
    }


def write_measurement(root: Path, implementation: str, report: dict) -> Path:
    """Serialize one complete implementation measurement to JSON."""

    output_path = root / f"{implementation}_measurement.json"
    output_path.write_text(f"{json.dumps(report, indent=2)}\n", encoding="utf-8")
    return output_path


def improvement_ratio(reference_value: float, bistreroc_value: float) -> float:
    """Return the reference resource use divided by BistreRoc resource use."""

    if bistreroc_value <= 0 or reference_value <= 0:
        raise ValueError("resource measurements must be positive")
    return float(reference_value / bistreroc_value)


def summarize(
    root: Path,
    output: Path,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> Path:
    """Validate results and write the stage-separated improvement matrix."""

    fixture = json.loads((root / "fixture.json").read_text(encoding="utf-8"))
    reports = {
        implementation: json.loads(
            (root / f"{implementation}_measurement.json").read_text(encoding="utf-8")
        )
        for implementation in IMPLEMENTATIONS
    }
    if reports["bistreroc"]["host"] != reports["reference"]["host"]:
        raise AssertionError("both implementations must be measured on the same host")
    for implementation, report in reports.items():
        configuration = report["configuration"]
        if (
            configuration["markers_per_cell_type"]
            != fixture["scenario"]["markers_per_cell_type"]
        ):
            raise AssertionError(
                f"{implementation} marker count differs from the prepared scenario"
            )
        if (
            configuration["q_value_threshold"]
            != fixture["scenario"]["q_value_threshold"]
        ):
            raise AssertionError(
                f"{implementation} q-value threshold differs from the prepared scenario"
            )
    comparison_columns = {
        "faster": "wall_seconds",
        "less_peak_memory": "peak_memory_bytes",
        "less_input_disk": "input_bytes",
        "less_output_disk": "output_bytes",
    }
    improvement_matrix = {
        stage: {
            label: improvement_ratio(
                reports["reference"]["stages"][stage][resource_field],
                reports["bistreroc"]["stages"][stage][resource_field],
            )
            for label, resource_field in comparison_columns.items()
        }
        for stage in STAGES
    }
    numerical_equivalence = {
        "relative_tolerance": relative_tolerance,
        "absolute_tolerance": absolute_tolerance,
        "maximum_absolute_difference": {
            "signature_generation": maximum_numeric_difference(
                root
                / Path(
                    reports["reference"]["stages"]["signature_generation"][
                        "primary_output"
                    ]
                ),
                root
                / Path(
                    reports["bistreroc"]["stages"]["signature_generation"][
                        "primary_output"
                    ]
                ),
                relative_tolerance,
                absolute_tolerance,
            ),
            "fraction_estimation": maximum_numeric_difference(
                root
                / Path(
                    reports["reference"]["stages"]["fraction_estimation"][
                        "primary_output"
                    ]
                ),
                root
                / Path(
                    reports["bistreroc"]["stages"]["fraction_estimation"][
                        "primary_output"
                    ]
                ),
                relative_tolerance,
                absolute_tolerance,
            ),
        },
    }
    summary = {
        "schema_version": 1,
        "recorded_at_utc": utc_timestamp(),
        "fixture": fixture,
        "improvement_matrix": improvement_matrix,
        "matrix_rows": list(STAGES),
        "matrix_columns": list(comparison_columns),
        "numerical_equivalence": numerical_equivalence,
        "raw": reports,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"{json.dumps(summary, indent=2)}\n", encoding="utf-8")
    return output


def parse_arguments(arguments: list[str] | None = None) -> argparse.Namespace:
    """Parse benchmark preparation, execution, and summarization commands."""

    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    prepare_parser.add_argument("--reference", type=Path, required=True)
    prepare_parser.add_argument("--mixture", type=Path, required=True)
    prepare_parser.add_argument("--mixture-count", type=int, default=8)

    run_parser = commands.add_parser("run")
    run_parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    run_parser.add_argument(
        "--implementation",
        choices=IMPLEMENTATIONS,
        required=True,
    )
    run_parser.add_argument(
        "--python-executable",
        type=Path,
        default=Path(sys.executable),
    )
    run_parser.add_argument("--cores", type=int, default=8)
    run_parser.add_argument("--sample-interval-ms", type=float, default=20.0)
    run_parser.add_argument("--reference-username")
    run_parser.add_argument("--reference-token")

    summarize_parser = commands.add_parser("summarize")
    summarize_parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    summarize_parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/results/realistic_8_mixtures.json"),
    )
    summarize_parser.add_argument("--relative-tolerance", type=float, default=5e-14)
    summarize_parser.add_argument("--absolute-tolerance", type=float, default=1e-12)
    return parser.parse_args(arguments)


def main() -> None:
    """Dispatch the requested benchmark operation."""

    arguments = parse_arguments()
    if arguments.command == "prepare":
        output_path = prepare_inputs(
            arguments.root,
            arguments.reference,
            arguments.mixture,
            arguments.mixture_count,
        )
    elif arguments.command == "run":
        root = arguments.root.resolve(strict=True)
        sample_interval_seconds = arguments.sample_interval_ms / 1000.0
        if sample_interval_seconds <= 0:
            raise ValueError("sample interval must be positive")
        if arguments.cores < 1:
            raise ValueError("cores must be positive")
        if arguments.implementation == "bistreroc":
            report = run_bistreroc(
                root,
                arguments.python_executable.absolute(),
                arguments.cores,
                CANONICAL_MARKER_COUNT,
                CANONICAL_Q_VALUE_THRESHOLD,
                sample_interval_seconds,
            )
        else:
            if (
                arguments.reference_username is None
                or arguments.reference_token is None
            ):
                raise ValueError(
                    "reference runs require --reference-username and --reference-token"
                )
            report = run_reference(
                root,
                arguments.reference_username,
                arguments.reference_token,
                CANONICAL_MARKER_COUNT,
                sample_interval_seconds,
            )
        output_path = write_measurement(root, arguments.implementation, report)
    else:
        output_path = summarize(
            arguments.root.resolve(strict=True),
            arguments.output,
            arguments.relative_tolerance,
            arguments.absolute_tolerance,
        )
    print(output_path.resolve())


if __name__ == "__main__":
    main()

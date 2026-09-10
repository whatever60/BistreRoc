"""Tests for stage-separated benchmark calculations and table validation."""

import json
import lzma
import sys
from pathlib import Path

import pytest

from benchmarks.pipeline_benchmark import (
    CANONICAL_MARKER_COUNT,
    CANONICAL_Q_VALUE_THRESHOLD,
    improvement_ratio,
    main,
    maximum_numeric_difference,
    parse_arguments,
    prepare_inputs,
    summarize,
    write_mixture_subset,
)


def test_write_mixture_subset_creates_matching_dense_and_xz_tables(
    tmp_path: Path,
) -> None:
    """Retain only the requested mixture columns in both representations."""

    source = tmp_path / "source.tsv"
    source.write_text(
        "Gene\tM1\tM2\tM3\nG1\t1\t2\t3\nG2\t4\t5\t6\n",
        encoding="utf-8",
    )
    dense = tmp_path / "mixture.tsv"
    compressed = tmp_path / "mixture.tsv.xz"

    assert write_mixture_subset(source, dense, compressed, 2) == 2
    expected = "Gene\tM1\tM2\nG1\t1\t2\nG2\t4\t5\n"
    assert dense.read_text(encoding="utf-8") == expected
    with lzma.open(compressed, "rt", encoding="utf-8") as input_handle:
        assert input_handle.read() == expected


def test_prepare_inputs_records_sparse_realistic_representations(
    tmp_path: Path,
) -> None:
    """Prepare equivalent dense and sparse inputs plus fixture metadata."""

    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "Gene\tA\tB\tB\nG1\t1\t0\t2\nG2\t0\t3\t0\n",
        encoding="utf-8",
    )
    mixture = tmp_path / "mixture.tsv"
    mixture.write_text(
        "Gene\tM1\tM2\tM3\nG1\t4\t5\t6\nG2\t7\t8\t9\n",
        encoding="utf-8",
    )
    root = tmp_path / "benchmark"

    fixture_path = prepare_inputs(root, reference, mixture, mixture_count=2)
    fixture = json.loads(fixture_path.read_text(encoding="utf-8"))

    assert fixture["scenario"]["reference_gene_count"] == 2
    assert fixture["scenario"]["reference_cell_count"] == 3
    assert fixture["scenario"]["cell_type_count"] == 2
    assert fixture["scenario"]["mixture_count"] == 2
    assert fixture["scenario"]["markers_per_cell_type"] == CANONICAL_MARKER_COUNT
    assert fixture["scenario"]["q_value_threshold"] == CANONICAL_Q_VALUE_THRESHOLD
    assert fixture["sources"]["reference"]["filename"] == "reference.tsv"
    assert fixture["sources"]["mixture"]["filename"] == "mixture.tsv"
    assert "path" not in fixture["sources"]["reference"]
    assert fixture["representations"]["reference"]["bistreroc_bytes"] > 0
    assert (root / "inputs/reference/reference.tsv").is_symlink()
    assert (root / "inputs/bistreroc/reference.bundle/matrix.mtx.xz").is_file()


def test_maximum_numeric_difference_accepts_compressed_float_drift(
    tmp_path: Path,
) -> None:
    """Compare compressed and plain tables within the declared tolerance."""

    reference = tmp_path / "reference.tsv"
    reference.write_text("ID\tA\nrow\t0.25\n", encoding="utf-8")
    candidate = tmp_path / "candidate.tsv.xz"
    with lzma.open(candidate, "wt", encoding="utf-8") as output_handle:
        output_handle.write("ID\tA\nrow\t0.2500000005\n")

    difference = maximum_numeric_difference(reference, candidate, 1e-12, 1e-9)

    assert difference == pytest.approx(5e-10)


def test_summarize_writes_three_by_four_improvement_matrix(tmp_path: Path) -> None:
    """Create all requested ratios and retain raw measurements in JSON."""

    root = tmp_path / "benchmark"
    root.mkdir()
    (root / "fixture.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "scenario": {
                    "mixture_count": 8,
                    "markers_per_cell_type": CANONICAL_MARKER_COUNT,
                    "q_value_threshold": CANONICAL_Q_VALUE_THRESHOLD,
                },
            }
        ),
        encoding="utf-8",
    )
    primary_outputs: dict[str, dict[str, Path]] = {}
    for implementation in ("reference", "bistreroc"):
        primary_outputs[implementation] = {}
        for stage in ("signature_generation", "fraction_estimation"):
            output = root / f"{implementation}_{stage}.tsv"
            output.write_text("ID\tA\nrow\t0.25\n", encoding="utf-8")
            primary_outputs[implementation][stage] = output
        multiplier = 4 if implementation == "reference" else 1
        report = {
            "host": {"system": "test-host"},
            "configuration": {
                "markers_per_cell_type": CANONICAL_MARKER_COUNT,
                "q_value_threshold": CANONICAL_Q_VALUE_THRESHOLD,
            },
            "stages": {
                stage: {
                    "wall_seconds": 8 * multiplier,
                    "peak_memory_bytes": 16 * multiplier,
                    "input_bytes": 32 * multiplier,
                    "output_bytes": 64 * multiplier,
                    **(
                        {
                            "primary_output": str(
                                primary_outputs[implementation][stage].relative_to(root)
                            )
                        }
                        if stage != "end_to_end"
                        else {}
                    ),
                }
                for stage in (
                    "signature_generation",
                    "fraction_estimation",
                    "end_to_end",
                )
            },
        }
        (root / f"{implementation}_measurement.json").write_text(
            json.dumps(report),
            encoding="utf-8",
        )

    output = tmp_path / "summary.json"
    summarize(root, output, 1e-12, 1e-9)
    result = json.loads(output.read_text(encoding="utf-8"))

    assert result["matrix_rows"] == [
        "signature_generation",
        "fraction_estimation",
        "end_to_end",
    ]
    assert result["matrix_columns"] == [
        "faster",
        "less_peak_memory",
        "less_input_disk",
        "less_output_disk",
    ]
    assert result["improvement_matrix"]["end_to_end"] == {
        "faster": 4.0,
        "less_peak_memory": 4.0,
        "less_input_disk": 4.0,
        "less_output_disk": 4.0,
    }

    reference_report_path = root / "reference_measurement.json"
    mismatched_reference_report = json.loads(
        reference_report_path.read_text(encoding="utf-8")
    )
    mismatched_reference_report["configuration"]["markers_per_cell_type"] += 1
    reference_report_path.write_text(
        json.dumps(mismatched_reference_report),
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="marker count differs"):
        summarize(root, output, 1e-12, 1e-9)


def test_reference_credentials_are_optional_for_other_commands() -> None:
    """Expose credential arguments only as optional reference-run inputs."""

    arguments = parse_arguments(["run", "--implementation", "bistreroc"])

    assert arguments.reference_username is None
    assert arguments.reference_token is None


def test_main_preserves_virtual_environment_python_symlink(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Keep the selected virtual environment when launching a benchmark."""

    benchmark_root = tmp_path / "benchmark"
    benchmark_root.mkdir()
    virtual_environment_python = tmp_path / "venv" / "bin" / "python"
    virtual_environment_python.parent.mkdir(parents=True)
    virtual_environment_python.symlink_to(Path(sys.executable).resolve())
    captured_python_executables: list[Path] = []

    def capture_run(
        root: Path,
        python_executable: Path,
        cores: int,
        marker_count: int,
        q_value_threshold: float,
        sample_interval_seconds: float,
    ) -> dict:
        """Capture the executable selected by the benchmark dispatcher."""

        captured_python_executables.append(python_executable)
        return {"implementation": "bistreroc"}

    monkeypatch.setattr("benchmarks.pipeline_benchmark.run_bistreroc", capture_run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pipeline_benchmark",
            "run",
            "--root",
            str(benchmark_root),
            "--implementation",
            "bistreroc",
            "--python-executable",
            str(virtual_environment_python),
        ],
    )

    main()

    assert captured_python_executables == [virtual_environment_python.absolute()]


def test_improvement_ratio_rejects_zero_measurements() -> None:
    """Reject undefined ratios instead of emitting misleading infinity."""

    with pytest.raises(ValueError):
        improvement_ratio(10, 0)

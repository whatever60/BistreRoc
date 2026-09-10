"""Integration tests for standalone fraction estimation."""

import gzip
import lzma
from pathlib import Path

import numpy as np
import pytest

from bistreroc.compressed_io import open_text_reader
from bistreroc.deconvolution import (
    estimate_fractions,
    fraction_estimation_script_path,
)


def test_fraction_estimation_rejects_nonpositive_core_count(
    tmp_path: Path,
) -> None:
    """Reject invalid resource controls before starting the R process."""

    with pytest.raises(ValueError, match="core_count must be a positive integer"):
        estimate_fractions(
            mixture_path=tmp_path / "missing_mixture.tsv",
            signature_path=tmp_path / "missing_signature.tsv",
            output_path=tmp_path / "fractions.tsv",
            core_count=0,
        )


def test_compressed_inputs_and_output_preserve_fraction_estimates(
    deconvolution_inputs: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    """Stream compressed matrices through the solver and compress its result."""

    signature, mixture = deconvolution_inputs
    compressed_signature = tmp_path / "signature.tsv.gz"
    compressed_mixture = tmp_path / "mixture.tsv.xz"
    compressed_signature.write_bytes(gzip.compress(signature.read_bytes()))
    compressed_mixture.write_bytes(lzma.compress(mixture.read_bytes()))
    output = tmp_path / "fractions.tsv.gz"
    estimate_fractions(
        mixture_path=compressed_mixture,
        signature_path=compressed_signature,
        output_path=output,
        core_count=2,
    )

    with open_text_reader(output) as input_handle:
        rows = [line.rstrip("\n").split("\t") for line in input_handle]
    assert rows[0] == ["Mixture", "A", "B", "C", "P-value", "Correlation", "RMSE"]
    estimates = np.asarray([row[1:4] for row in rows[1:]], dtype=float)
    np.testing.assert_allclose(estimates.sum(axis=1), 1, rtol=0, atol=1e-14)
    assert estimates[0, 0] > 0.7
    assert estimates[1, 2] > 0.7


def test_fraction_estimation_preserves_signature_gene_order(
    tmp_path: Path,
) -> None:
    """Keep solver observations in signature order during gene alignment."""

    signature = tmp_path / "signature.tsv"
    signature.write_text(
        "NAME\tA\tB\tC\n"
        "C1\t1\t1\t100\n"
        "A1\t100\t1\t1\n"
        "B1\t1\t100\t1\n"
        "C2\t1\t1\t80\n"
        "A2\t80\t1\t1\n"
        "B2\t1\t80\t1\n",
        encoding="utf-8",
    )
    mixture = tmp_path / "mixture.tsv"
    mixture.write_text(
        "Gene\tmostly_a\tmostly_c\n"
        "B2\t8.9\t8.9\n"
        "C2\t8.9\t64.2\n"
        "A2\t64.2\t8.9\n"
        "B1\t10.9\t10.9\n"
        "C1\t10.9\t80.2\n"
        "A1\t80.2\t10.9\n",
        encoding="utf-8",
    )
    output = tmp_path / "fractions.tsv"

    estimate_fractions(mixture, signature, output)

    rows = [
        line.rstrip("\n").split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
    ]
    assert rows[0] == ["Mixture", "A", "B", "C", "P-value", "Correlation", "RMSE"]
    assert [row[0] for row in rows[1:]] == ["mostly_a", "mostly_c"]
    np.testing.assert_allclose(
        np.asarray([row[1:] for row in rows[1:]], dtype=float),
        np.asarray(
            [
                [
                    0.93753608768318764,
                    0.03123195556642290,
                    0.03123195675038935,
                    9999.0,
                    0.99967694244964822,
                    0.03889261279918781,
                ],
                [
                    0.03123195675038945,
                    0.03123195556642289,
                    0.93753608768318764,
                    9999.0,
                    0.99967694244964822,
                    0.03889261279918792,
                ],
            ]
        ),
        rtol=5e-14,
        atol=1e-12,
    )


def test_fraction_kernel_uses_deterministic_accumulation_order() -> None:
    """Protect floating-point accumulation order required for reproducibility."""

    script = fraction_estimation_script_path().read_text(encoding="utf-8")
    assert "row_major_signature_values <- c(t(signature_matrix))" in script
    assert script.count("Reduce(`+`,") == 4
    assert "mean(signature_matrix)" not in script
    assert "scale(mixture_matrix)" not in script

"""Tests for the public command-line surface."""

from pathlib import Path

from bistreroc.cli import parse_arguments


def test_deconvolution_command_uses_explicit_output() -> None:
    """Parse a complete fraction-estimation request."""

    parsed_arguments = parse_arguments(
        [
            "deconvolve",
            "--mixture",
            "mixture.tsv.xz",
            "--signature",
            "signature.tsv.gz",
            "--output",
            "fractions.tsv.gz",
            "--cores",
            "8",
        ]
    )

    assert parsed_arguments.command == "deconvolve"
    assert parsed_arguments.mixture_path == Path("mixture.tsv.xz")
    assert parsed_arguments.output_path == Path("fractions.tsv.gz")
    assert parsed_arguments.core_count == 8


def test_signature_command_exposes_marker_search() -> None:
    """Parse descriptive marker-count controls."""

    parsed_arguments = parse_arguments(
        [
            "build-signature",
            "--reference",
            "reference.bundle",
            "--output",
            "signature.tsv",
            "--minimum-markers-per-cell-type",
            "8",
            "--maximum-markers-per-cell-type",
            "12",
        ]
    )

    assert parsed_arguments.minimum_markers_per_cell_type == 8
    assert parsed_arguments.maximum_markers_per_cell_type == 12


def test_reference_conversion_command() -> None:
    """Expose sparse conversion as one self-contained subcommand."""

    parsed_arguments = parse_arguments(
        ["convert-reference", "reference.tsv.gz", "reference.bundle"]
    )

    assert parsed_arguments.source_path == Path("reference.tsv.gz")
    assert parsed_arguments.destination_path == Path("reference.bundle")

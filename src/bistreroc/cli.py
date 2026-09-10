"""Command-line interface for BistreRoc."""

import argparse
from importlib.metadata import version
from pathlib import Path

from bistreroc.deconvolution import estimate_fractions
from bistreroc.signature_generation import build_signature
from bistreroc.sparse_reference import dense_integer_tsv_to_sparse_bundle


def add_resource_control_arguments(command_parser: argparse.ArgumentParser) -> None:
    """Add resource controls shared by computational commands."""

    command_parser.add_argument("--cores", dest="core_count", type=int, default=3)
    command_parser.add_argument("--timeout-seconds", type=int, default=3600)


def parse_arguments(arguments: list[str] | None = None) -> argparse.Namespace:
    """Parse BistreRoc command-line arguments."""

    root_parser = argparse.ArgumentParser(
        prog="bistreroc",
        description=(
            "Build cell-type signatures and estimate cell fractions from "
            "bulk expression mixtures."
        ),
    )
    root_parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('bistreroc')}",
    )
    command_parsers = root_parser.add_subparsers(dest="command", required=True)

    deconvolution_parser = command_parsers.add_parser(
        "deconvolve",
        help="estimate cell fractions in one or more expression mixtures",
    )
    deconvolution_parser.add_argument(
        "--mixture",
        dest="mixture_path",
        type=Path,
        required=True,
    )
    deconvolution_parser.add_argument(
        "--signature",
        dest="signature_path",
        type=Path,
        required=True,
    )
    deconvolution_parser.add_argument(
        "--output",
        dest="output_path",
        type=Path,
        required=True,
    )
    add_resource_control_arguments(deconvolution_parser)

    signature_parser = command_parsers.add_parser(
        "build-signature",
        help="learn a signature matrix from labeled single-cell counts",
    )
    signature_parser.add_argument(
        "--reference",
        dest="reference_path",
        type=Path,
        required=True,
    )
    signature_parser.add_argument(
        "--output",
        dest="output_path",
        type=Path,
        required=True,
    )
    signature_parser.add_argument(
        "--minimum-markers-per-cell-type",
        type=int,
        default=300,
    )
    signature_parser.add_argument(
        "--maximum-markers-per-cell-type",
        type=int,
        default=500,
    )
    signature_parser.add_argument("--q-value-threshold", type=float, default=0.01)
    add_resource_control_arguments(signature_parser)

    conversion_parser = command_parsers.add_parser(
        "convert-reference",
        help="convert a dense count table to a sparse reference bundle",
    )
    conversion_parser.add_argument("source_path", metavar="source", type=Path)
    conversion_parser.add_argument(
        "destination_path",
        metavar="destination",
        type=Path,
    )
    conversion_parser.add_argument(
        "--matrix-compression",
        choices=("gzip", "bzip2", "xz"),
        default="xz",
    )
    return root_parser.parse_args(arguments)


def main() -> None:
    """Dispatch one BistreRoc command."""

    parsed_arguments = parse_arguments()
    if parsed_arguments.command == "convert-reference":
        dense_integer_tsv_to_sparse_bundle(
            parsed_arguments.source_path,
            parsed_arguments.destination_path,
            matrix_compression=parsed_arguments.matrix_compression,
        )
        print(parsed_arguments.destination_path.resolve())
        return
    if parsed_arguments.command == "deconvolve":
        output_path = estimate_fractions(
            mixture_path=parsed_arguments.mixture_path,
            signature_path=parsed_arguments.signature_path,
            output_path=parsed_arguments.output_path,
            core_count=parsed_arguments.core_count,
            timeout_seconds=parsed_arguments.timeout_seconds,
        )
    else:
        output_path = build_signature(
            reference_path=parsed_arguments.reference_path,
            output_path=parsed_arguments.output_path,
            minimum_markers_per_cell_type=(
                parsed_arguments.minimum_markers_per_cell_type
            ),
            maximum_markers_per_cell_type=(
                parsed_arguments.maximum_markers_per_cell_type
            ),
            q_value_threshold=parsed_arguments.q_value_threshold,
            core_count=parsed_arguments.core_count,
            timeout_seconds=parsed_arguments.timeout_seconds,
        )
    print(output_path)


if __name__ == "__main__":
    main()

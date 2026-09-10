"""Distribution-content tests."""

from pathlib import Path


def test_distribution_contains_only_runtime_source_modules() -> None:
    """Keep the installed source surface focused on runtime modules."""

    package_directory = Path(__file__).parents[1] / "src" / "bistreroc"
    assert sorted(path.name for path in package_directory.glob("*.py")) == [
        "__init__.py",
        "cli.py",
        "compressed_io.py",
        "deconvolution.py",
        "signature_generation.py",
        "sparse_reference.py",
    ]
    assert sorted(path.name for path in (package_directory / "r").glob("*.R")) == [
        "estimate_fractions.R",
        "select_signature_genes.R",
    ]

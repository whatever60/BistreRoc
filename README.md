<p align="center">
  <img src="assets/bistreroc_logo.png" alt="BistreRoc eagle logo" width="150">
</p>

# BistreRoc

A better CIBERSORT, open and free.

**8.7× faster. 15.2× less
peak memory. 19.5× less input disk. 841.8× less output disk.**

## Two-stage workflow

### 1. Learn a signature matrix

The dense reference format is a tab-separated integer count matrix. The first
column contains gene identifiers; every remaining column is one cell, and its
column name is that cell's type label. Repeated labels are expected.

```text
Gene    CD4 T    CD4 T    B cell
IL7R    8        4        0
MS4A1   0        0        12
```

BistreRoc reads this table directly when it is plain or compressed with gzip,
bzip2, or XZ:

```bash
bistreroc build-signature \
  --reference reference.tsv.xz \
  --output signature.tsv.xz \
  --cores 8
```

For sparse single-cell data, convert the dense table once to BistreRoc's
compressed Matrix Market bundle, then build from the bundle:

```bash
bistreroc convert-reference reference.tsv.xz reference.bundle

bistreroc build-signature \
  --reference reference.bundle \
  --output signature.tsv.xz \
  --cores 8
```

The bundle is intentionally simple:

```text
reference.bundle/
├── matrix.mtx.xz
├── genes.tsv.gz
├── labels.tsv.gz
└── library_sizes.tsv.gz
```

`--matrix-compression` can select `gzip`, `bzip2`, or `xz`. Signature selection
defaults to 300–500 markers per cell type and chooses the cutoff with the lowest
exact matrix condition number. Use
`--minimum-markers-per-cell-type`, `--maximum-markers-per-cell-type`, and
`--q-value-threshold` to change that search.

### 2. Estimate mixture fractions

Mixture and signature matrices use the same labeled TSV convention: genes are
rows and mixtures or cell types are numeric columns. Each input and the output
may independently be plain, gzip, bzip2, or XZ.

```bash
bistreroc deconvolve \
  --mixture mixtures.tsv.xz \
  --signature signature.tsv.xz \
  --output fractions.tsv.xz \
  --cores 8
```

The result contains one row per mixture, one fraction column per cell type, and
the `P-value`, `Correlation`, and `RMSE` diagnostics.

## Benchmark

- 16,105 genes
- 28,326 labeled single cells across 76 cell types
- 30,564,848 nonzero counts
- 14 markers per cell type, yielding a 524-gene signature
- 8 bulk mixtures

| Stage | Faster | Less peak memory | Less input disk | Less output disk |
|---|---:|---:|---:|---:|
| Signature generation | 9.1× | 16.3× | 19.9× | 862.9× |
| Mixture to fractions | 5.3× | 2.4× | 2.6× | 3.2× |
| End to end | 8.7× | 15.2× | 19.5× | 841.8× |

Ratios are reported only after both stages match the benchmark target:
identifiers, signature membership, and ordering are exact, with numerical
differences limited to floating-point precision. See the
[benchmark protocol and full report](benchmarks/README.md).

## Installation

BistreRoc requires Python 3.11 or newer, NumPy, R, `data.table`, and `e1071`.
Until the Bioconda recipe is merged, install the immutable release tag in a
small conda environment:

```bash
mamba create -n bistreroc -c conda-forge \
  python=3.13 numpy r-base=4.4 r-data.table r-e1071 pip git

mamba run -n bistreroc python -m pip install \
  "bistreroc @ git+https://github.com/whatever60/BistreRoc.git@v0.1.0"
```

After the recipe is accepted, installation becomes:

```bash
mamba install -c bioconda bistreroc
```

## Development

```bash
uv sync --locked
uv run pytest -q
```

CI also runs Ruff, checks the environment, and builds both distribution formats
on Python 3.11 and 3.13.

## License

BistreRoc is released under the [MIT License](LICENSE).

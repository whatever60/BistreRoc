# Benchmarks

## Comparison fixture

For reproducibility, the comparison target is pinned to
`docker.io/cibersortx/fractions@sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc`.
Maintainers who run it copy these target files verbatim into the ignored local
`dev/reference_image/` fixture:

- `/src/CIBERSORTxFractions`;
- `/src/R_modules/run_SVR.R`;
- `/src/R_modules/run_adjust_w_pseudo.R`;
- `/src/R_modules/run_combatba.R`;
- `/src/R_modules/run_heatmap.R`;
- `/src/R_modules/run_normalize_quantiles.R`;
- `/src/R_modules/run_optimization.R`;
- `/src/filter/CCLE_nonBlood_Avgs.txt`; and
- `/src/filter/XavierGEP.NonImmune.txt`.

This fixture is used only to run the optional comparison. None of its contents
is tracked, packaged, imported, or required by BistreRoc. The two BistreRoc
stages are implemented in `src/bistreroc/r/estimate_fractions.R` and
`src/bistreroc/r/select_signature_genes.R`; compressed and sparse I/O,
orchestration, the CLI, tests, and benchmark runner live in the tracked Python
sources.

## Canonical workload

The recorded canonical workload uses the real sequencing data at realistic
scale. `prepare` accepts any equivalent labeled reference and mixture paths:

- a 16,105-gene, 28,326-cell reference with 76 cell types;
- the sparse coordinate matrix compressed with XZ;
- the first eight bulk mixtures;
- 14 markers per cell type, yielding roughly 500 unique signature genes; and
- separate measurements for signature generation and fraction estimation.

Reference conversion and mixture subsetting are one-time fixture preparation.
They are recorded but excluded from timed stages. Each implementation consumes
its normal stored representation: compressed sparse input for BistreRoc and
uncompressed dense input for the comparison image.

## Run

From a source checkout with its development environment active:

```bash
uv run python -m benchmarks.pipeline_benchmark prepare \
  --reference /path/to/labeled_single_cell_counts.tsv \
  --mixture /path/to/bulk_mixtures.tsv
uv run python -m benchmarks.pipeline_benchmark run \
  --implementation bistreroc
```

The reference run is optional and is the only command accepting credentials:

```bash
uv run python -m benchmarks.pipeline_benchmark run \
  --implementation reference \
  --reference-username "$REFERENCE_USERNAME" \
  --reference-token "$REFERENCE_TOKEN"
```

After both runs finish, validate their numeric tables and build the comparison:

```bash
uv run python -m benchmarks.pipeline_benchmark summarize
```

Working data and logs are written below the ignored
`benchmark_results/realistic_8_mixtures/` directory. The default final report is
`benchmarks/results/realistic_8_mixtures.json`.

## Report contract

The JSON report contains a 3×4 `improvement_matrix`. Its rows are signature
generation, fraction estimation, and end to end. Its columns are:

- `faster`: reference wall time divided by BistreRoc wall time;
- `less_peak_memory`: reference peak aggregate RSS divided by BistreRoc peak;
- `less_input_disk`: reference stored inputs divided by BistreRoc stored inputs;
- `less_output_disk`: reference retained outputs divided by BistreRoc outputs.

Peak memory is sampled from each Linux process tree every 20 ms. For the
container run, the container init process is discovered through Podman and
sampled in addition to the client process. End-to-end wall time is the sum of
the two sequential stages; its peak memory is the larger stage peak. Fractions
stage input includes its learned signature, whereas end-to-end input includes
only the external reference and mixture.

Summarization requires identical row and column labels and finite values. It
uses `rtol=5e-14` and `atol=1e-12` by default, so the report is not produced
when differences exceed floating-point-scale drift.

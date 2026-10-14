# Stochastic-mapping validation

Existing scripts and audits, separate from ML recovery:

- `benchmark-sim-on-tree-classe-td.R`: simulation-on-backbone benchmark.
- `benchmark-conditioned-backbone.R`: standalone conditional-backbone prototype.
- `benchmark-extinction-grid.R`: extinction-grid benchmark.
- `validate-direct-extinction.R`: direct extinction sampler checks.
- The three Markdown documents preserve historical plans and validation evidence, not new validation results from this reorganization.

Run from the repository root, for example:

```sh
Rscript inst/benchmarks/stochastic-mapping/benchmark-sim-on-tree-classe-td.R
```

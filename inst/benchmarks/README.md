# DivExtra benchmarks

This is the single benchmark home, separate from the small checks in `tests/testthat/`.

| Folder | Purpose |
|---|---|
| [01-likelihood-validation](01-likelihood-validation/README.md) | Compare fixed-parameter likelihoods without optimization |
| [02-ml-constant-rates](02-ml-constant-rates/README.md) | Recover three linked rates using episodic ClaSSE |
| [03-ml-one-rate-change](03-ml-one-rate-change/README.md) | Recover four parameters when one rate changes between epochs |
| [04-ml-all-rates-change](04-ml-all-rates-change/README.md) | Recover six epoch-specific rates |
| [stochastic-mapping](stochastic-mapping/README.md) | Validate sampled histories and extinct lineages |
| [archive](archive/README.md) | Preserve earlier pilot workflows and outputs |

## Run organization

Each numbered family contains named `runs/` directories. Keep a run's readable scripts, annotated settings, manifest, datasets, provenance, results, summaries, plots and report together. The publication run is self-contained when copied. Evidence-only historical imports are labelled in their `IMPORT.md`; they are not complete raw-data reproductions. Avoid unnecessary shared wrappers; see [shared](shared/README.md).

## Current publication run

`02-ml-constant-rates/runs/publication-100/` contains the 100-replicate-per-scenario-and-tip-group analysis. From that directory, `bash sync-roihu.sh` retrieves Roihu results and updates status; it does not submit jobs. Final analysis runs only when all expected fit records are present and inference jobs are terminal. Remote paths are unchanged by this reorganization.

## Git and publication data

Track scripts, settings, manifests, compact summaries, reports and selected plots. Large raw trees, optimizer results, logs and source archives are retained locally but ignored by Git. Before publication, deposit the complete run in a data repository, record its identifier and checksums, and verify retrieval. A code-only clone does not contain all raw data.

The entire `inst/benchmarks/` folder is excluded from package builds by `.Rbuildignore`, so these large research benchmarks are not installed with the R package. Tiny package tests remain in `tests/testthat/`.

## Migration

The duplicate top-level `benchmarks/` was consolidated here. Older GeoSSE scripts and outputs moved to `archive/pilots/geosse-parameter-recovery/`; script paths were updated, not scientific settings. Original project folders remain intact. See [migration record](archive/migration.md).

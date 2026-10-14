# Migration record — 14 October 2026

| Original project folder | New GitHub-relative destination | Scope |
|---|---|---|
| benchmarks/publication/constant-rates-100 | inst/benchmarks/02-ml-constant-rates/runs/publication-100 | Complete local run copy plus fresh Roihu retrieval; raw files Git-ignored |
| roihu-constant-size-50 | inst/benchmarks/02-ml-constant-rates/runs/pilot-50-evidence | Original report and selected summary tables |
| roihu-constant-size-50/preparation/likelihood-checks.csv | inst/benchmarks/01-likelihood-validation/runs/constant-50-evidence | Historical fixed-parameter comparisons |
| roihu-speciation-size-50 | inst/benchmarks/03-ml-one-rate-change/runs/speciation-pilot-50 | Original report and selected summaries |
| roihu-dispersal-contraction-size-50 | inst/benchmarks/03-ml-one-rate-change/runs/dispersal-contraction-pilot-50 | Original report and selected summaries |
| roihu-epoch-grid-x010-50 | inst/benchmarks/04-ml-all-rates-change/runs/epoch-grid-x010-pilot-50 | Original report, stored truths and selected summaries |
| GitHub inst/benchmarks/geosse-parameter-recovery | inst/benchmarks/archive/pilots/geosse-parameter-recovery | Historical workflow moved; script paths updated |

Project source: divextra-sunken. Original folders are retained, including the initial publication copy. New downloads are directed to the reorganized publication folder. Remote paths and running scientific jobs are unchanged. The temporary top-level GitHub benchmarks directory was consolidated into inst/benchmarks; its empty directories were removed. Standalone mapping scripts and audits moved into inst/benchmarks/stochastic-mapping.

No Git commit, push, external deposition or raw-data deletion is performed by this migration. Evidence-only folders are not complete standalone reproductions. The publication run's raw files remain available locally but are excluded from ordinary Git tracking.

Current publication jobs: smoke 2300270; inference 2300275; analysis 2300286. Check live status rather than assuming the copied submission note describes completion.

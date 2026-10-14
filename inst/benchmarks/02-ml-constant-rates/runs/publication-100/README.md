# Publication benchmark: constant rates through episodic ClaSSE

## Aim and design

Test maximum-likelihood recovery using **DivExtra `make.classe.td()` with two epochs**, a known boundary at 50% of prescribed simulation stem duration before present, and all rates linked across epochs. There are exactly three estimated parameters: symmetric speciation s, dispersal d and contraction/local extinction x. No classical GeoSSE comparison is run.

| Scenario | s | d | x | Proposal durations |
|---|---:|---:|---:|---|
| Slow | 0.01 | 0.01 | 0.005 | 300, 400, 500, 600 |
| Medium | 0.1 | 0.1 | 0.05 | 30, 40, 50, 60 |
| Fast | 0.3 | 0.3 | 0.15 | 10, 13.3333, 16.6667, 20 |

100 independent selected trees per scenario × tip-size group: 50–99, 100–249, 250–499, 500–999 and 1,000–2,000 extant tips. Total: **1,500 trees and 4,500 fits** (three starts each).

Rates and durations are inversely scaled; these scenarios represent equivalent processes on different time scales. They are not tests of intrinsically faster versus slower diversification. Durations are model time units, or Ma if rates are per Ma.

## Simulation and selection

`simulate.R` makes a fresh fixed pool of 800 independent proposals per scenario-duration combination (9,600 total), using deterministic unique seeds. It uses DivExtra `make.tree.classe.td()` with identical parameter rows, old then recent. One starting stem and a uniformly sampled root state A/B/AB; complete extant sampling. All attempt statuses and state counts are saved.

Retain only full-duration histories with 50–2,000 extant tips. There is **no crossing or minimum-state-count filter** in this constant-rate control. The 2,000-living-lineage safety cap may stop a tree before its endpoint; these histories have unknown final sizes and are excluded. This is path-dependent selection beyond a final-size filter. Trees are not stopped when their desired final size is reached.

Select 100 without replacement in every scenario-size group from eligible proposals; stop with an error if any pool is insufficient. Store full simulation tables and reconstructed crown trees, tip states, truths, seeds and durations. Initial stems are not included in the likelihood tree. Boundary position is based on prescribed stem duration; diagnostics record whether it falls before the reconstructed crown. Do not silently change the boundary or discard such trees.

The simulator copy includes the same narrowly scoped unsplit-living-stem safeguard as the pilots: replace the early `is.null(tb.e1)` return with `n.taxa == 0`. Its exact definition is saved in provenance/simulator-used.R. Installed package files are not modified.

## Likelihood and optimization

The annotated YAML describes both epochs. `model.R` uses the episodic likelihood, explicitly expands each of the three free rates into identical old/recent rates, and uses ROOT.OBS, survival conditioning, complete sampling, tol=1e-10 and eps=1e-8. No biological root state is fixed during inference.

Three fixed jittered starts per tree, seeded independently from simulation; starts do not use generating truth. Log-scale nlminb, bounds 1e-10–100, iter.max=1000, eval.max=3000 and rel.tol=1e-9. Preserve all starts and select the highest finite primary likelihood. Near-zero boundary diagnostics are saved separately, not silently substituted. Each worker has a 20-minute process timeout including boundary diagnostics.

## Local validation

`validate.R` checks 100 trees per group and seed uniqueness, then tests one tree per group at truth and after perturbing each free parameter. Because the rates are linked, moving the bookkeeping boundary from 50% to 25% must leave the likelihood unchanged within numerical tolerance. These are validation evaluations only, not extra fitted epoch scenarios. Nine optimization smoke tests cover all three rate scales. `preflight.R` compares remote runtime function definitions and test likelihoods with local records. No `make.geosse` fit or evaluation is used.

## Reproduce locally

Run from this directory. Do not regenerate over an existing manifest; work in a fresh copy without generated outputs when reproducing simulations.

```bash
Rscript simulate.R
Rscript validate.R
```

Install dependencies `ape`, `diversitree`, `yaml`, `callr` and dependencies required by DivExtra. An archived DivExtra source tree is included in provenance/divextra-source.tar.gz, from clean Git commit 4aa1513ca3c1f9aa493f0ffae779cf67505f2e8e. It can be installed using `R CMD INSTALL provenance/divextra-source.tar.gz` after its dependencies are installed. Package versions, session information and critical runtime definitions are recorded in provenance/. The archived source and runtime checks distinguish reproducibility from a version label alone. Compiled dependencies must be installed for the target operating system; no private keys or credentials are included.

Inference for a single task can also run locally, for example:

```bash
SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=4 Rscript run-task.R
```

## Roihu

Use the existing project library and r-env/461. The array has 150 tasks, ten trees/30 fits per task, four workers, 8 GB, three hours, and at most six concurrent allocations. Three smoke tasks (1, 51, 101) cover Slow/Medium/Fast. Submit the remaining tasks with a dependency on successful smoke completion. A small dependent summary job runs analyse.R after all inference tasks complete.

See submission-status.md for actual job IDs before submitting anything: **do not create duplicate arrays**. Remote directory: /scratch/project_2004996/geosse-benchmark/publication-constant-100.

```bash
Rscript preflight.R
sbatch --array=1,51,101%3 submit.sh
# After recording the returned SMOKE_ID:
sbatch --array=2-50,52-100,102-150%6 --dependency=afterok:SMOKE_ID submit.sh
# After recording MAIN_ID:
sbatch --dependency=afterok:SMOKE_ID:MAIN_ID summarise-job.sh
```

## Outputs and interpretation

After downloading results, run `Rscript analyse.R` to produce full-range and focused recovery boxplots, recovery-versus-tip-size plots, numerical summaries, convergence flags and a short report.md. Focused plots span 0–3 times truth and count estimates above the displayed range; every estimate remains in box calculations. Lower/upper summary bounds are central 95% replicate ranges, not per-fit confidence intervals or coverage.

Simulation tables, trees, manifests, proposals, scripts, annotated YAML, source snapshot and seeds are all within this folder. No pilot data are required. The current results directory contains only completed fits; absence must not be interpreted as successful inference. Convergence code 0 is not proof of a global optimum or reliable recovery.

Code uses explicit statements and loops following the readable-R style. The three-start worker is the only process-management infrastructure; scientific choices remain visible in the scripts.

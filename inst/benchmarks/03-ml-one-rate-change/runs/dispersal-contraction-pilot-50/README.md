# Dispersal and revised contraction recovery

## Design

Two scenarios, each with 750 trees: 50 per boundary (5%, 10%, 25%) and extant-tip group (50–99, 100–249, 250–499, 500–999, 1,000–2,000). Each tree has two epochs, not three. Boundary ages are known and fixed to their generating ages, measured from present relative to simulated stem duration.

- Dispersal: old (s,d,x) = (0.1,0.1,0.05); recent (0.1,0.3,0.05). Estimate shared s, shared x, recent d and old d.
- Revised contraction: old (0.1,0.1,0.005); recent (0.1,0.1,0.05). Estimate shared s, shared d, recent x and old x.

Uniform random root state among A, B and AB, one starting stem, complete extant sampling. Retention: reach full duration, 50–2,000 extant tips and at least 30 reconstructed boundary-crossing lineages. State counts are recorded, not an additional filter, consistent with the immediately preceding screens. A 2,000-living-lineage safety cap introduces additional path-dependent selection; capped histories have unknown final size.

## Generation

Reuse eligible dispersal trees from ../remaining-scenarios-pilot and revised contraction trees from ../contraction-revised-pilot. Fill each group's shortage at its most productive pilot duration among 30,40,50,60, using fresh seeds and a 10,000-proposal limit per group. Retain other eligible trees in the pool. Select 50 without replacement per group. This is an adaptively allocated, size-selected benchmark, not an unconditional tree distribution. Every proposal and selected tree's source are recorded. No earlier trees or results are overwritten.

Use divextra::make.tree.classe.td with the local-only unsplit-living-stem safeguard used in the preceding pilots. No installed package edits. Generating parameter rows are old then recent; likelihood epoch 1 is recent. Save full tables and reconstructed crown trees with state data, truths, links, generating ages and diagnostics.

## Inference and validation

1,500 trees × three fixed jittered starts = 4,500 primary fits. Four estimated parameters per fit. Starts are not centered on generating truth. Same ROOT.OBS treatment, survival conditioning, full sampling, tol=1e-10, eps=1e-8, log-scale nlminb and parameter bounds as the previous benchmark. Save near-zero boundary diagnostics separately from primary results. Classical GeoSSE is not equivalent to these scenarios and is not fitted.

Validation checks counts, seed uniqueness, crossing/node identities and representative likelihoods in all 30 groups. Compare each linked likelihood with explicit six-rate expansion at truth and after perturbing each free parameter (150 comparisons). Run three local starts on one small tree per scenario. Compare remote versus local truth likelihoods before launch.

There are 300 logical batches of five trees/15 fits. Roihu permits at most 200 submitted small-partition jobs, so submit.sh runs only the two smoke batches (1 and 11), and submit-main.sh pairs the other 298 batches into 149 array tasks (ten trees/30 fits each). run-chunk.R maps these tasks to the original batch IDs; results retain their original task-001 through task-300 directory names. Each allocation has four workers, 8 GB and three-hour walltime. Each fit has a 20-minute process timeout including boundary diagnostics. Main tasks depend on successful smoke completion. Smoke tasks fail the dependency gate if any saved fit is erroneous or has an invalid likelihood. Main-array concurrency is six. Do not submit duplicate arrays. Results are checkpointed per fit; failures and timeouts remain explicit.

## Run and interpretation

Locally run `Rscript generate.R`, then `Rscript validate.R`. Generation refuses to overwrite manifest.csv. See preparation/submission-status.md for actual job IDs before any submission. Upload datasets and inference inputs, run preflight.R using the established Roihu project R library, then submit smoke and dependent arrays.

After completion, report all 4,500 attempts, convergence and three-start agreement, parameter recovery by scenario/boundary/size, near-zero estimates, and associations with tip counts, crossings and epoch exposure. These correlated diagnostics do not establish causality. Empirical replicate ranges are not per-fit confidence-interval coverage. Conditional selection and correct-model/full-sampling assumptions limit generalization.

Scripts follow the readable-R preference for explicit steps and loops; established likelihood/worker functions are retained where needed.

# Two-epoch all-rates-change pilot: common percentage boundaries

## Requested design

Two epochs, not three epochs: alternative known boundaries at 5%, 10%, or 25% of prescribed total simulation duration before present.

- S-I: old `(s,d,x)=(0.1,0.1,0.05)`, recent `(0.3,0.3,0.10)`.
- S-II: old `(0.3,0.3,0.10)`, recent `(0.1,0.1,0.05)`.
- All six rates estimated independently; symmetric between regions within an epoch.
- Target 50 trees per scenario x boundary x tip group: 100–249, 250–499, 500–999, 1000–2000.
- Retain complete-duration trees with at least 30 reconstructed boundary crossings. No additional tip-state or recent-node filter. The simulation cap is 2000 living lineages, which introduces path-dependent censoring.
- Random uniform root state during simulation; ROOT.OBS and survival conditioning during inference, complete sampling, known epoch boundary.

## Generation and reuse

`generate.R` uses divextra `make.tree.classe.td` and `pars.geosse2classe`, with old parameters in row 1. It retains the prior local unsplit-living-stem safeguard; package files are not edited. Tables are converted and extinct tips pruned with the same diversitree conversion as earlier pilots.

Existing accepted histories from `all-rates-x010-pilot` are reused. The exact 300 completed S-II trees with 500–2000 tips from `roihu-balanced-rates-50` are preserved and are NOT submitted again. `manifest.csv` records their original indices in `previous.tree.index`; their task column is NA. Their fits remain in the original directory, linked by this explicit mapping.

For each new cell, use the most productive screened duration among existing qualifying trees. If there are none, use the shortest duration with qualifying trees at that boundary. Every new proposal is recorded. Stop after 3000 new attempts for a target cell if 50 qualifying trees cannot be obtained. Never relax the crossing threshold or silently submit fewer than 50 trees in a cell: unfilled cells are omitted and explicitly reported. Retained top-up histories can fill other tip groups at the same boundary. Generation saves pool checkpoints after each cell.

Duration calibration and conditioning differ between groups; comparisons do not isolate the causal effect of N. Smaller S-I/25% cells were very rare or absent in the screening. The achieved design must be read from the manifest, not assumed to equal the nominal 1200 trees.

## Inference and submission

Reuse the validated six-rate likelihood and optimizer from the preceding pilot, three fixed jittered starts per tree (not centered on truth), nlminb on log rates with bounds 1e-10–100, tolerance settings unchanged. Near-zero diagnostic fits remain separate. A 20-minute timeout applies per worker fit including boundary diagnostics.

Five new trees per Slurm task, four worker processes, 8 GB, 3-hour walltime, at most six concurrent tasks. The actual array range is supplied at submission from the manifest. Tasks 1 and 21 are smoke tests (S-I and S-II); remaining tasks are submitted with afterok dependency. Remote preflight must first match local reference likelihoods. Do not rerun submission commands for an already recorded job.

## Reproduce and inspect

Run `Rscript generate.R` in this directory (refuses to overwrite manifest), then `Rscript validate.R`. Run `Rscript preflight.R` on Roihu after transferring scripts, inputs, manifest, and preparation/likelihood-checks.csv. All tasks operate only on rows with nonmissing task numbers.

`preparation/new-proposals.csv`, `preparation/topup-design.csv`, and `preparation/eligible-pool.csv` document simulation effort and duration selection. `preparation/submission-status.md` records the actual retained counts and job IDs once submitted.

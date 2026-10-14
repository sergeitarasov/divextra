# Publication benchmark submission — 14 October 2026

## Prepared

Fresh proposal pool: 9,600 histories. Of these, 7,263 reached the endpoint with survivors, 2,069 went extinct and 268 were safety-cap stopped. There were 4,913 eligible full-duration trees with 50–2,000 extant tips. Selected 1,500 distinct seeds, exactly 100 per rate scenario and size group. No pilot trees were reused. All selected trees happen to have at least five tips per state; no such filter was imposed.

The fixed epoch boundary is half the prescribed stem duration. Nine selected trees have this boundary before their reconstructed crown; they remain included and are identified in tree-diagnostics.csv. No crossing filter is used in this linked-rate control.

## Validation

- All 15 groups and 60 boundary-invariance/parameter-perturbation likelihood comparisons passed.
- Nine local optimization starts converged; within-tree likelihood spreads were below 3e-10.
- Archived clean DivExtra source commit 4aa1513ca3c1f9aa493f0ffae779cf67505f2e8e. All five critical function definitions match the local installed runtime; the documented simulator safeguard is saved separately.
- Roihu runtime definitions matched local definitions. Three reference likelihood differences were below 6e-13.
- Upload checksum verification found no differences before launch.

## Jobs

- **2300270**: three smoke tasks (1, 51, 101). All completed with exit 0. All 90 fits were downloaded and checked: valid likelihoods, three parameters, convergence code 0. Maximum within-tree three-start likelihood spread: 2.11e-7.
- **2300275**: main array, 147 tasks excluding smoke tasks; submitted with afterok:2300270. Six allocations were RUNNING at the post-submission check; the remaining tasks were pending under the six-task concurrency limit.
- **2300286**: dependent analysis job, afterok:2300270:2300275. It will run analyse.R to save numerical summaries, full/focused recovery boxplots, recovery-versus-tip-size plots and report.md. Its successful completion has not yet been verified.

Total: 150 inference tasks, 1,500 trees, 4,500 primary fits. Ten trees and three starts per tree per task; four workers, 8 GB, three-hour allocation, 20-minute timeout per fit. No classical GeoSSE comparison is included.

Remote folder: /scratch/project_2004996/geosse-benchmark/publication-constant-100.

Check with `squeue -j 2300275,2300286`; after jobs leave the queue use `sacct -j 2300270,2300275,2300286 --allocations --format=JobID,State,ExitCode,Elapsed`.

Do not submit duplicates. This is a submission/validation record, not a full recovery result. Main fits, summary outputs and Slurm accounting must be retrieved after completion. No recurring monitoring automation was created.

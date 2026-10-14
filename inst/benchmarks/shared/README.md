# Shared conventions

Keep reusable helpers minimal. Current runs retain their validated likelihood and worker functions beside their inputs for standalone reproducibility. Do not refactor a running job's scientific functions during migration.

Shared conventions: deterministic recorded seeds; known root/sampling conditioning; save all starts and failures; distinguish optimizer agreement from parameter recovery; treat replicate ranges separately from per-fit confidence intervals. Slurm submission templates remain within each versioned run because resource limits and tasks differ.

No new shared scientific implementation is introduced by this reorganization.

# Roihu publication benchmark: completed 14 October 2026

The constant-rate publication run is complete and retrieved into `inst/benchmarks/02-ml-constant-rates/runs/publication-100/` in the DivExtra checkout. All 150 inference array tasks and the summary job completed with exit code 0. No jobs were resubmitted.

## Completeness and convergence

- 1,500 trees: 100 in each of three scenarios by five tip-size groups.
- 4,500/4,500 primary fits valid; 4,498 returned optimizer code 0; no timeouts.
- All 1,500 selected best fits have optimizer code 0. Strict three-start agreement holds for 1,496 trees (99.73%).
- Flagged trees: 161 and 480 (Slow), 1039 and 1071 (Fast). Their best available fits remain in recovery summaries. Tree 480 has a particularly poor third start; multi-start fitting avoided selecting it. A converged fit is not proof of a global optimum.
- Separate zero-boundary diagnostics: 357 fits across 120 trees; 282 returned code 0. Maximum likelihood improvement over the corresponding primary fit was only 2.30e-6. These diagnostics do not replace primary estimates.

## Recovery

Median absolute relative error for 1,000--2,000 tips:

| Scenario | Speciation s | Dispersal d | Contraction x |
|---|---:|---:|---:|
| Slow | 3.4% | 6.0% | 14.4% |
| Medium | 3.7% | 4.9% | 14.5% |
| Fast | 3.3% | 6.2% | 15.5% |

With only 50--99 tips, median relative errors were 11--15% for s, 21--25% for d, and 77--98% for x. In the largest group, 72--74% of x estimates were within 25% of truth, compared with 11--15% in the smallest group. Thus increasing tree size improves recovery substantially, but x remains the least precise parameter. The proportional Slow/Medium/Fast designs show similar recovery; this is not evidence that an absolute rate scale intrinsically improves inference.

The plots and numerical summaries use all best fits. Central 95% replicate ranges are distributions of estimates, not per-fit confidence intervals or coverage tests. The benchmark tests three rates linked across two epochs, not independently estimated epoch rates. Results are conditional on the prescribed simulation durations, size selection and living-lineage cap.

## Runtime and integration

Median primary-fit elapsed time was 2.36 seconds; the 95th percentile was 15.83 seconds and maximum 176.76 seconds. Slurm's largest reported MaxRSS was 1,216,180 KiB (about 1.16 GiB); this is the accounting-reported maximum, not a summed concurrent-worker memory measurement.

`report.md`, `summary/recovery-summary.csv`, `summary/flagged-fits.csv`, the recovery plots, final accounting and the checksum inventory are saved in the run folder. Raw fit records remain local and Git-ignored; compact summaries and plots are Git-visible. Integration means files in the local GitHub checkout; no commit, push or external data deposition was performed.

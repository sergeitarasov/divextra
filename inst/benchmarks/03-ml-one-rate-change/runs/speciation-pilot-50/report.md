# Speciation-change recovery pilot

13 October 2026. Jobs2264359 and2264371. All130 tasks completed with exit code0; all1950 primary fits downloaded. Reproduce tables and plots with `Rscript analyse.R` from this directory.

## Main result

The fitted direction of change (recent s > old s) was correct in649/650 trees (99.85%). This is descriptive recovery of direction, not significance, model-selection power or confidence-interval coverage. Recent speciation and shared dispersal are estimated accurately on larger trees; old speciation is less precise, particularly at the25% boundary. Shared contraction/local extinction x remains the weakest parameter. At1000–2000 tips, median x errors remain28–31%, and all four rates are within25% of truth on only38–46% of trees.

## Design

Old s=.1; recent s=.3; shared d=.1 and x=.05. Boundary ages fixed at their generating ages:5%,10%,25% of total stem duration backward from present.650 unique trees,50 per represented boundary-size group. The25% groups below250 tips were omitted by design. Complete sampling; known model; random starting state. Each tree was fitted from three starts, with best primary likelihood retained. Near-zero-boundary diagnostics are kept separately.

All trees reached full duration, had50–2000 extant tips and at least30 reconstructed boundary crossings. Existing independent pilot trees were supplemented with3042 fresh proposals targeted at shortages. This is a survival/size/crossing-selected, adaptively allocated sample. The2000 living-lineage simulation cap additionally excludes paths that hit that cap before present. Duration and geographic state counts are recorded. No state-count filter or inference of the boundary age was performed.

## Recovery by boundary and size

Median absolute relative errors, abs(estimate/truth - 1), expressed as percentages. Each row represents50 trees.

| Boundary | Tips | Recent s | Old s | Shared d | Shared x | All four within25% |
|---|---|---:|---:|---:|---:|---:|
| 5% | 50–99 | 15.5 | 20.3 | 25.6 | 100.0 | 8% |
| 5% | 100–249 | 6.8 | 19.5 | 17.9 | 92.9 | 6% |
| 5% | 250–499 | 5.9 | 12.2 | 13.7 | 52.8 | 20% |
| 5% | 500–999 | 4.7 | 11.0 | 7.6 | 44.6 | 32% |
| 5% | 1000–2000 | 1.6 | 6.9 | 6.0 | 27.9 | 46% |
| 10% | 50–99 | 11.9 | 13.6 | 15.2 | 100.0 | 4% |
| 10% | 100–249 | 7.4 | 15.9 | 11.9 | 100.0 | 10% |
| 10% | 250–499 | 5.7 | 15.3 | 11.7 | 61.7 | 16% |
| 10% | 500–999 | 4.1 | 11.7 | 7.2 | 41.5 | 24% |
| 10% | 1000–2000 | 3.0 | 9.5 | 4.3 | 30.9 | 46% |
| 25% | 250–499 | 3.6 | 17.6 | 5.8 | 59.1 | 18% |
| 25% | 500–999 | 3.6 | 15.9 | 5.8 | 36.3 | 28% |
| 25% | 1000–2000 | 2.2 | 14.4 | 4.4 | 30.4 | 38% |

Near-zero x estimates (less than1% of true x) occur in40% and56% of50–99-tip trees at5% and10%, respectively. At1000–2000 tips these proportions are0%,2%,4% at5%,10%,25%. A median error close to100% does not necessarily mean the median estimate is zero: both underestimates and overestimates contribute to absolute error.

## Tip counts, crossings and boundary placement

In the1000–2000-tip group, median reconstructed crossings were525 at5%,223 at10%, and43.5 at25%. Corresponding median old-s errors were6.9%,9.5%,14.4%. Thus similar present-day richness can conceal very different amounts of older-tree information. This is consistent with the importance of crossing counts, but is not an isolated causal comparison: durations, selected tree distributions and recent-epoch exposure also differ.

Separate regressions used log1p(relative absolute error) against log2(N) or log2(crossings), controlling for boundary x duration cells, with HC3 robust standard errors. One observation per tree per rate; no pseudoreplication across starts. Both predictors were negatively associated with errors for every parameter (all adjusted p<1e-9). They were fitted separately and should not be interpreted as independent effects.

| Rate | Raw Spearman rho with N | Adjusted N slope | Robust95% interval | Adjusted crossing slope | Robust95% interval |
|---|---:|---:|---|---:|---|
| Recent s | -0.506 | -0.0237 | -0.0287 to -0.0186 | -0.0226 | -0.0279 to -0.0174 |
| Old s | -0.217 | -0.0586 | -0.0758 to -0.0414 | -0.0546 | -0.0716 to -0.0377 |
| Shared d | -0.458 | -0.0428 | -0.0539 to -0.0316 | -0.0425 | -0.0539 to -0.0311 |
| Shared x | -0.520 | -0.1379 | -0.1703 to -0.1055 | -0.1346 | -0.1673 to -0.1020 |

Slopes are per doubling of the predictor on the transformed-error scale, not percentage-point changes in error. Strict-fit sensitivity results are retained in summary/information-regressions.csv. Missing small-tree25% groups prevent a complete factorial comparison across all sizes.

## Numerical performance

All1950 attempts had finite non-penalty likelihoods; no errors or timeouts.1949 returned code0. The single nonzero-code start was tree527,start1 (1258tips,5%); successful alternatives provided the retained best fit. All650 selected best fits have code0.

645/650 trees pass the strict criterion: all three starts converged and likelihood spread<1e-4. The five flagged trees are32,79,176,177,527, with likelihood spreads0.0521,0.5032,0.0589,0.0473,0.00573. Four therefore have multiple solutions despite all starts reporting convergence. These histories were not discarded from the main recovery summaries; best primary likelihood was used and sensitivity results are saved.

352 separate zero-boundary checks:264 code0,88 nonzero. None improved its corresponding primary start by more than1e-4 (maximum gain5.08e-6). These diagnostics did not replace primary fits.

Median fit elapsed time was4.339seconds;95th percentile26.677seconds; maximum619.794seconds, below the20-minute timeout. Task durations ranged9–624seconds. Final Slurm accounting is saved locally.

## Reliability conclusion

For this selected correct-model benchmark, the speciation increase is almost always recovered in direction, and recent s and d have good typical accuracy from a few hundred tips onward. Old s benefits from more crossing lineages and is not equally precise across boundaries at a given N. The≥30-crossing rule improves the design's information floor but does not guarantee precise estimates of all parameters. Shared x remains unreliable on small trees and noticeably uncertain even on large trees. Do not claim universal recovery or a minimum-N guarantee.

The benchmark does not test unknown epoch ages, alternative rate contrasts, incomplete sampling, model misspecification, or extinct-lineage stochastic mapping. All accuracy summaries are conditional on the simulation selection rules. No per-fit confidence intervals were computed, so empirical95% replicate ranges in the tables are not confidence-interval coverage. Regression confidence intervals concern associations only.

## Figures and outputs

- [Parameter recovery](plots/parameter-recovery.png): red dashed lines are true rates; boxes show medians/IQR and1.5IQR whiskers with outliers. Scales differ between panels. Empty25% cells are explicitly labeled omitted.
- [Recovery versus tip count](plots/recovery-versus-tips.png): median error per50-tree group for each boundary; dotted line marks25% error, not a significance threshold.
- summary/recovery-by-group.csv includes means, medians and empirical central95% ranges.
- summary/starts.csv, best-fits.csv, flagged-fits.csv, joint-accuracy.csv, information-regressions.csv and boundary-checks.csv retain diagnostics.

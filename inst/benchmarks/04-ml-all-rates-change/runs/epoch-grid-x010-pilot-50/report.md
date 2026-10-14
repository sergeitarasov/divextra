# GeoSSE two-epoch recovery: 5%, 10%, 25%

Completed 13 October 2026.

## Main result

Speciation is generally recovered much better than contraction. Larger trees help, but tree size alone does not solve the shortage of old-epoch information in S-I. At 25%, old contraction is near zero in 76% of the largest S-I trees, and none of those 50 trees recover all six rates within 50% of truth. S-I at 5% and S-II at 10% offer the most balanced large-tree recovery among these tested percentage settings, without making contraction precise on every tree.

## Design and completeness

S-I old (s,d,x)=(0.1,0.1,0.05), recent=(0.3,0.3,0.10); S-II reverses the epochs. Truths were read from each stored dataset. Six independent rates, known boundary at 5%, 10%, or 25% of prescribed simulation duration before present. Four tip groups, 50 trees per included cell. The S-I/25%/100–249 cell was not simulated for inference because no qualifying trees were obtained; it is not a zero-recovery result.

All 170 new Slurm tasks completed with exit 0. All 2550 new fits were finite. Together with 300 reused S-II trees, there are 1150 trees and 3449/3450 finite fits; the one missing finite result is the previously recorded reused timeout. No fits were rerun. Every tree has a usable estimate. Explicit original-to-current index mapping is retained in preparation/reused-results.csv.

## Optimization and numerical checks

Primary fits: 3447 code 0, 2 code 1, one timeout. All selected highest-likelihood fits returned code 0. Strict three-start agreement (three finite fits, all code 0, likelihood spread <1e-4) holds for 1105/1150 trees (96.1%). All three parameter vectors agree within 5% of truth for every rate in 1112 trees.

There are 1305 separate near-zero diagnostic fits: 656 code 0 and 649 code 1. Maximum likelihood improvement over their primary fits is 0.000808083. Boundary checks are not substituted for primary fits, and their convergence codes are reported separately. The strict-convergence sensitivity changes group median absolute relative errors by at most 12.7 percentage points. Full details are in summary/flagged-fits.csv and boundary-checks.csv.

## Recovery on 1000–2000-tip trees

Median absolute relative error (%), 50 trees in each row. Smaller is better.

setting | s_recent | d_recent | x_recent | s_old | d_old | x_old 
--- | --- | --- | --- | --- | --- | --- 
S-I-b5 | 3.5 | 8 | 35.2 | 10.2 | 18.9 | 42.7 
S-I-b10 | 3 | 5.2 | 31.5 | 17.5 | 40.3 | 88.9 
S-I-b25 | 3.1 | 6.7 | 26.9 | 18.9 | 69.8 | 100 
S-II-b5 | 5 | 21.1 | 77.4 | 5.7 | 9.1 | 31.8 
S-II-b10 | 4.2 | 8.8 | 30.7 | 4.4 | 6.2 | 30.3 
S-II-b25 | 3.9 | 6.4 | 21.5 | 6.4 | 23.7 | 48.3 

S-I old-contraction error rises from 42.7% at 5% to 88.9% at 10% and approximately 100% at 25%. Recent contraction errors instead decline from 35.2% to 31.5% and 26.9%. S-II at 10% gives approximately 30% error for both contraction rates; 5% weakens recent contraction, while 25% weakens old contraction. This is an information-allocation trade-off, not simply optimizer nonconvergence.

## All-six-rate accuracy across sizes

Percent of trees with every rate within 50% of truth. Missing means not simulated, not 0%.

setting | 100-249 | 250-499 | 500-999 | 1000-2000 
--- | --- | --- | --- | --- 
S-I-b5 | 0 | 10 | 32 | 38 
S-I-b10 | 0 | 6 | 6 | 16 
S-I-b25 | NA | 0 | 0 | 0 
S-II-b5 | 0 | 4 | 8 | 26 
S-II-b10 | 2 | 12 | 20 | 36 
S-II-b25 | 2 | 14 | 20 | 44 

Even at the largest size, only 0–16% of trees recover all six rates within 25% of truth (see joint-accuracy.csv). Strong optimizer agreement therefore does not establish reliable recovery of all parameters.

## Tips, crossings, and recent branching

setting | tip.group | median.N | median.crossing | median.recent.nodes 
--- | --- | --- | --- | --- 
S-I-b5 | 100-249 | 158 | 64 | 96 
S-I-b5 | 250-499 | 335.5 | 121.5 | 213 
S-I-b5 | 500-999 | 664.5 | 236 | 445 
S-I-b5 | 1000-2000 | 1359.5 | 441 | 918.5 
S-I-b10 | 100-249 | 172 | 40.5 | 131 
S-I-b10 | 250-499 | 343.5 | 60 | 286 
S-I-b10 | 500-999 | 750 | 105 | 635 
S-I-b10 | 1000-2000 | 1297.5 | 184 | 1121 
S-I-b25 | 250-499 | 446 | 32 | 410 
S-I-b25 | 500-999 | 655 | 37 | 618 
S-I-b25 | 1000-2000 | 1241.5 | 45.5 | 1195 
S-II-b5 | 100-249 | 171.5 | 154.5 | 16 
S-II-b5 | 250-499 | 339 | 299.5 | 37.5 
S-II-b5 | 500-999 | 660.5 | 573 | 87 
S-II-b5 | 1000-2000 | 1291.5 | 1134 | 160 
S-II-b10 | 100-249 | 145.5 | 117 | 26 
S-II-b10 | 250-499 | 348.5 | 279.5 | 70.5 
S-II-b10 | 500-999 | 709.5 | 529 | 170.5 
S-II-b10 | 1000-2000 | 1432 | 1076.5 | 380 
S-II-b25 | 100-249 | 153.5 | 88.5 | 65 
S-II-b25 | 250-499 | 391 | 200 | 194 
S-II-b25 | 500-999 | 716 | 346 | 377.5 
S-II-b25 | 1000-2000 | 1371 | 634 | 730.5 

Exploratory separate-predictor regressions use log(1+absolute relative error) versus log2 information, adjusted for setting and prescribed duration, with HC3 intervals and BH-adjusted p-values. Old-contraction error decreases with tips and crossings in both scenarios. Slopes below describe a doubling of the predictor on the transformed-error scale, not percent changes in raw error.

scenario | predictor | rho | slope | lower | upper | p.BH 
--- | --- | --- | --- | --- | --- | --- 
S-I | N | -0.336 | -0.202 | -0.267 | -0.137 | 6.90511435810681e-09 
S-I | crossing | -0.388 | -0.209 | -0.275 | -0.142 | 5.00439681722207e-09 
S-II | N | -0.412 | -0.151 | -0.202 | -0.099 | 3.10642318089283e-08 
S-II | crossing | -0.449 | -0.148 | -0.2 | -0.097 | 5.09379825429037e-08 

These are associations within a conditionally selected sample, not independent causal effects. Crossings, old nodes (crossings minus one), tips, and branch exposure are strongly correlated. Durations were calibrated differently across cells. The retention rule requires at least 30 reconstructed crossings, and the 2000-living-lineage cap adds path-dependent selection. Regression CSVs include all six rates, six information predictors and strict-convergence sensitivity.

## Runtime and uncertainty

New finite-fit elapsed time: median 9.3 seconds, maximum 561.6 seconds. New task elapsed time: median 1.2 minutes, maximum 13.0 minutes. Maximum Slurm-reported step MaxRSS: 1193.2 MiB; this is not a sum of concurrent worker peaks.

The empirical 2.5–97.5% replicate ranges in recovery-by-group.csv use quantile type 7. They are not per-fit confidence intervals; confidence-interval coverage was not measured. Box plots show medians and quartiles, 1.5-IQR whiskers and explicit outliers. Full-range plots retain every outlier; focused ratio panels show 0–3 and state how many values exceed that limit.

## Plots and reproducibility

Run `Rscript analyse.R`, followed by `Rscript report.R`, from this directory. Scripts follow the readable-R style; original inputs and results are preserved.

- [S-I recovery versus tip groups](plots/S-I-recovery-overview.png)
- [S-II recovery versus tip groups](plots/S-II-recovery-overview.png)
- [Joint accuracy](plots/joint-accuracy.png)
- [Tree information](plots/tree-information.png)

- S-I boundary 5%: [rates](plots/S-I-b5-rates.png), [full ratios](plots/S-I-b5-ratios-full.png), [focused ratios](plots/S-I-b5-ratios-focused.png).
- S-I boundary 10%: [rates](plots/S-I-b10-rates.png), [full ratios](plots/S-I-b10-ratios-full.png), [focused ratios](plots/S-I-b10-ratios-focused.png).
- S-I boundary 25%: [rates](plots/S-I-b25-rates.png), [full ratios](plots/S-I-b25-ratios-full.png), [focused ratios](plots/S-I-b25-ratios-focused.png).
- S-II boundary 5%: [rates](plots/S-II-b5-rates.png), [full ratios](plots/S-II-b5-ratios-full.png), [focused ratios](plots/S-II-b5-ratios-focused.png).
- S-II boundary 10%: [rates](plots/S-II-b10-rates.png), [full ratios](plots/S-II-b10-ratios-full.png), [focused ratios](plots/S-II-b10-ratios-focused.png).
- S-II boundary 25%: [rates](plots/S-II-b25-rates.png), [full ratios](plots/S-II-b25-ratios-full.png), [focused ratios](plots/S-II-b25-ratios-focused.png).
- [S-I error versus N](plots/S-I-error-vs-N.png): descriptive LOWESS, not adjusted regression.
- [S-I error versus crossing](plots/S-I-error-vs-crossing.png): descriptive LOWESS, not adjusted regression.
- [S-I error versus recent.nodes](plots/S-I-error-vs-recent.nodes.png): descriptive LOWESS, not adjusted regression.
- [S-II error versus N](plots/S-II-error-vs-N.png): descriptive LOWESS, not adjusted regression.
- [S-II error versus crossing](plots/S-II-error-vs-crossing.png): descriptive LOWESS, not adjusted regression.
- [S-II error versus recent.nodes](plots/S-II-error-vs-recent.nodes.png): descriptive LOWESS, not adjusted regression.

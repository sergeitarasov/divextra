# Why old contraction is poorly recovered

## Diagnosis

The current evidence points mainly to weak practical information about old x in extant trees, not failed optimization or ODE integration. This is not a proof of structural nonidentifiability of the constrained GeoSSE model.

Across the 30 trees where x_old is independently estimated, 21 estimates are below 1e-5 despite a generating value of 0.02. GeoSSE x is local extinction: it removes a single-area lineage entirely, or contracts AB to A/B. In the symmetric model AB has total contraction hazard 2x. It is therefore not solely a geographic transition rate.

## Likelihood profiles

One median-tip-count tree was chosen in each of the six relevant scenario-boundary cells, independently of recovery error. At each of seven x_old values (including exact zero and truth 0.02), all other independent rates were optimized from two starts: the previous fit and generating nuisance values. All 42 selected profile solutions returned convergence code zero.

| Selected tree | Log-likelihood loss at true x_old = 0.02 |
|---|---:|
| All rates change, 5%, replicate 2 | 0.0115 |
| All rates change, 25%, replicate 4 | 0.3749 |
| All rates change, 50%, replicate 5 | 0.1115 |
| Contraction only, 5%, replicate 4 | 1.1088 |
| Contraction only, 25%, replicate 5 | 0.5422 |
| Contraction only, 50%, replicate 2 | 0.1969 |

The generating rate remains compatible with all six profiles despite five of these fits preferring the zero boundary. A near-zero point estimate is not evidence that extinction was absent. The horizontal 1.92 line in the figure is a conventional likelihood-loss reference, not a coverage-calibrated confidence threshold for a boundary parameter or this selected sample.

No profiled likelihood improved on the original fit by more than 1.1e-6. Changing ODE tolerance across 1e-8, 1e-10 and 1e-12 changed likelihood by at most 1.44e-6. These checks do not support numerical error or a missed substantially better optimum as the main explanation. Some profiles are broad even up to x_old = 0.2.

## What the complete histories reveal

For all-rates-change at the 5% boundary, median old counts are 20 terminal extinctions and 5 AB contractions. Using complete histories, the median count/exposure estimate is 0.01992, versus an extant-tree median near zero. The denominator is time in A/B plus twice time in AB, including the simulated stem. This diagnostic uses information unavailable to ordinary inference and is not an unbiasedness claim for the selected sample.

At 50%, median old AB-contraction counts are zero in both relevant scenarios, and median old terminal extinctions are only one (all-rates-change) or three (contraction-only). Eight of the 30 full histories contain no old AB-contraction event. The current filters ensure branching information but not many events informative about x. Even when events occur, terminal extinct lineages are missing from the observed tree and geographic histories are latent.

## Selection is an additional concern, not yet a quantified explanation

Selection on 50--1000 tips and minimum old/recent branching events changes the sampled histories. The all-rates-change 50% cell needed 500 attempts for five accepted trees. Current inference conditions on survival but does not include a correction for these additional selection criteria. Compare matched filtered and unfiltered ensembles at the same generating rates and durations to measure this contribution; these six profiles alone cannot attribute the zero estimates to filtering.

## Recommended next steps

1. Report profile-based uncertainty and likelihood-supported ranges, not only point-estimate errors. Calibrate coverage by simulation because zero is a boundary. Do not treat 0.02 versus an estimate near 0 as decisive model failure without examining uncertainty.
2. Benchmark linked versus separate x across epochs as explicit alternative models. Linking can increase precision but is misspecified when contraction truly changes; evaluate that bias rather than assuming linking is correct. Partial pooling is another regularized option and requires prior-sensitivity checks.
3. Run matched filtered/unfiltered controls before tightening information thresholds. Increasing tip/crossing minima cannot guarantee contraction information and may amplify selection effects.
4. If biologically justified, combine independent clades under shared or hierarchically related contraction rates. This supplies independent histories; more optimization starts do not supply new information. Historical geographic/fossil observations could add information but require an appropriate likelihood extension.
5. Propagate the remaining parameter uncertainty into extinct-lineage counts and arrival-time summaries. Simulations conditional only on the near-zero MLE omit this source of uncertainty; a likelihood grid is a sensitivity analysis, not automatically a posterior distribution.

The existing simulator and inference functions have not been changed.

## Files and reproduction

- [Profiles](profiles.pdf)
- [Numerical/profile results](profiles.csv)
- [Complete-history diagnostic](full-history-diagnostic.csv)
- Run `R/17-diagnose-old-contraction.R` from the DivExtra repository root. Existing profiles are reused.

## Relevant literature

[Louca and Pennell (2021)](https://pubmed.ncbi.nlm.nih.gov/34019824/) discuss zero extinction estimates from extant timetrees. Their general birth-death result is relevant background, not a direct identifiability theorem for our state-dependent finite-epoch model.

[Pennell, Sarver and Harmon (2012)](https://arxiv.org/abs/1208.3225) show how selecting unusually large clades can distort diversification inference. The effect of our particular filters still needs its own controlled test.

# Random-tree fixed-parameter likelihood validation

Run `Rscript run.R` from this directory. No optimization, cluster submission or parameter estimation is performed. `Rscript summarise.R` regenerates the summary and figure from the saved evaluations.

## Design

- Seed: 20261014; 1,000 independent random ultrametric trees from `ape::rcoal`.
- Tip count: discrete uniform 100--500; crown age: uniform 1--30 time units, obtained by rescaling each tree.
- Tip states: a shuffled approximately balanced allocation of A, B and AB; all three states are present. These numerical test data are not generated conditional on the evaluation rates and do not test biological model adequacy.
- GeoSSE: seven independently log-uniform rates between 0.005 and 0.3 per time unit (sA, sB, sAB, xA, xB, dA, dB). No regional symmetry restrictions.
- To evaluate GeoSSE using DivExtra's episodic implementation, its seven rates are expanded internally with `pars.geosse2classe`. There is no separate GeoSSE-restricted diversitree ClaSSE comparison.
- Traditional unrestricted ClaSSE: 27 independently log-uniform rates between 0.001 and 0.1.
- DivExtra: two identical epochs, with a fixed boundary at half the crown age. This is 55 supplied entries (one boundary plus two 27-entry vectors), not 55 independent rates.
- Complete sampling, ROOT.OBS; test both `condition.surv = FALSE` and `TRUE`. GeoSSE state order is AB,A,B versus ClaSSE A,B,AB; tip codes are converted explicitly.
- ODE controls: tol=1e-10, eps=1e-8. Prespecified pass rule: finite paired log likelihoods and absolute difference <=1e-5.

## Comparisons

1. DivExtra episodic ClaSSE restricted to GeoSSE versus diversitree GeoSSE (7 independent rates).
2. DivExtra episodic versus traditional diversitree unrestricted ClaSSE (27 independent rates).

There is no separate DivExtra `make.geosse` evaluator here: the user's GeoSSE implementation is the GeoSSE restriction of episodic ClaSSE. Each tree supplies two comparisons under both survival settings, giving 4,000 paired checks. This corrected run supersedes the earlier 10--100-tip, three-comparison run; all trees and evaluations were generated afresh.

## Reproducibility and scope

`inputs.rds` contains every tree, tip state, rate vector and boundary. `results/likelihood-checks.csv` retains every comparison, including errors rather than discarding them. `summary.csv`, `report.md` and `plots/likelihood-agreement.pdf` are the outputs. Source files for the episodic implementation and conversion are copied into `provenance/` and used by the script; installed package versions are recorded in session-info.txt. Install divextra and its dependencies before running.

This small validation bundle deliberately includes its inputs and evaluation table in Git, unlike large inference datasets. The reproducible random draws follow the R version and RNG behavior recorded in session-info.txt.

Agreement validates identical-epoch splitting over this tested range. It does not demonstrate correctness of genuinely different rates across epochs, optimization, parameter recovery, or stochastic mapping. No optimization-based GeoSSE comparison is added to the constant-rate publication recovery experiment.

## Numerical failures in the corrected run

Cases 462 (474 tips) and 888 (353 tips) produced `node stack overflow` in both the DivExtra GeoSSE parameterization and traditional diversitree GeoSSE, with and without survival conditioning. A fresh-process rerun using `Rscript recheck-errors.R` reproduced all four paired failures; its outputs are saved separately and do not replace the primary results. Their underlying cause has not been isolated. All 998 other GeoSSE cases and all 1,000 unrestricted ClaSSE cases passed under both settings. Failed cases are retained, not excluded from the denominator.

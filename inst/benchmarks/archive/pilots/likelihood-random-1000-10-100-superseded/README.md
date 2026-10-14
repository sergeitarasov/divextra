# Random-tree fixed-parameter likelihood validation

Run `Rscript run.R` from this directory. No optimization, cluster submission or parameter estimation is performed. `Rscript summarise.R` regenerates the summary and figure from the saved evaluations.

## Design

- Seed: 20261014; 1,000 independent random ultrametric trees from `ape::rcoal`.
- Tip count: discrete uniform 10--100; crown age: uniform 1--30 time units, obtained by rescaling each tree.
- Tip states: a shuffled approximately balanced allocation of A, B and AB; all three states are present. These numerical test data are not generated conditional on the evaluation rates and do not test biological model adequacy.
- GeoSSE: seven independently log-uniform rates between 0.005 and 0.3 per time unit (sA, sB, sAB, xA, xB, dA, dB). No regional symmetry restrictions.
- GeoSSE-as-ClaSSE: the same seven independent rates expanded into 27 entries with `pars.geosse2classe`; constrained entries and structural zeros are intentional, not 27 independent rates.
- Additional unrestricted ClaSSE check: 27 independently log-uniform rates between 0.001 and 0.1.
- DivExtra: two identical epochs, with a fixed boundary at half the crown age. This is 55 supplied entries (one boundary plus two 27-entry vectors), not 55 independent rates.
- Complete sampling, ROOT.OBS; test both `condition.surv = FALSE` and `TRUE`. GeoSSE state order is AB,A,B versus ClaSSE A,B,AB; tip codes are converted explicitly.
- ODE controls: tol=1e-10, eps=1e-8. Prespecified pass rule: finite paired log likelihoods and absolute difference <=1e-5.

## Comparisons

1. DivExtra episodic ClaSSE restricted to GeoSSE versus diversitree GeoSSE (7 independent rates).
2. That same episodic model versus diversitree ClaSSE with the equivalent GeoSSE mapping (27 entries, 7 independent rates).
3. DivExtra episodic versus diversitree unrestricted ClaSSE (27 independent rates).

There is no separate DivExtra `make.geosse` evaluator here: the user's GeoSSE implementation is the GeoSSE restriction of episodic ClaSSE. Each tree supplies all three comparisons under both survival settings, giving 6,000 paired checks. The same reference implementation is not independent evidence twice: comparisons 1 and 2 test two interfaces to the equivalent model.

## Reproducibility and scope

`inputs.rds` contains every tree, tip state, rate vector and boundary. `results/likelihood-checks.csv` retains every comparison, including errors rather than discarding them. `summary.csv`, `report.md` and `plots/likelihood-agreement.pdf` are the outputs. Source files for the episodic implementation and conversion are copied into `provenance/` and used by the script; installed package versions are recorded in session-info.txt. Install divextra and its dependencies before running.

This small validation bundle deliberately includes its inputs and evaluation table in Git, unlike large inference datasets. The reproducible random draws follow the R version and RNG behavior recorded in session-info.txt.

Agreement validates identical-epoch splitting over this tested range. It does not demonstrate correctness of genuinely different rates across epochs, optimization, parameter recovery, or stochastic mapping. No optimization-based GeoSSE comparison is added to the constant-rate publication recovery experiment.

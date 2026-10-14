# Fixed-time GeoSSE recovery pilot

This pilot uses five independently simulated datasets in each of five scenarios.
Each dataset was fitted with the six-rate episodic model from three starting vectors.
Starting vectors use eps = 0.1, 0.4, 0.8 and independent multiplicative perturbations between 0.5 and 2.
Inference uses ROOT.OBS and survival conditioning. Each optimizer has a 90-second wall-time limit; three processes run concurrently.
Optimizer code zero and finite likelihood are required to select a best converged estimate.
Likelihood agreement means all three starts converge and their log-likelihood range is below 0.0001.
Parameter relative spread is also recorded: equal likelihoods do not establish equal parameters or identifiability.

## Duration calibration

For each scenario, 20 independent attempts were simulated for each of ten duration combinations.
Candidates required at least eight usable trees and no taxon-limit truncations in calibration.
The score targeted 30 reconstructed lineages at the boundary and approximately 200 extant tips.
The chosen durations were frozen before independent recovery simulations. No filtering on realized tip count or boundary lineage count was applied.
Trees lacking a reconstructed crown spanning both epochs were excluded; this is additional design conditioning.
A 3000-taxon cap protects computation. This small calibration cannot guarantee caps will never occur.

```
 scenario_id old.time recent.time
         P01       25           6
         P02       35           6
         P03       25           6
         P04       20           6
         P05       30           6
```

## Completed fits

Converged starts: 75 of 75
Datasets with three agreeing converged starts: 18 of 25
Datasets whose best converged likelihood is below truth by more than 0.0001: 0

```
 scenario datasets converged.starts all.three.agree median.N median.cross
      P01        5               15               1      220           75
      P02        5               15               5      716          106
      P03        5               15               4       67           23
      P04        5               15               5       48           22
      P05        5               15               3      212           26
 median.old.length
          411.9117
          751.9997
          145.8062
          129.3343
          264.5687
```

## All three rates change: P05

Old rates: s=0.10, d=0.03, x=0.02. Recent rates: s=0.30, d=0.20, x=0.15.

```
 scenario parameter n truth       median       mean      minimum   maximum
      P05  s_recent 5  0.30 3.413217e-01 0.31695450 2.455417e-01 0.3655495
      P05  d_recent 5  0.20 2.087394e-01 0.23178858 1.244290e-01 0.3315471
      P05  x_recent 5  0.15 2.155644e-01 0.21863675 2.858944e-08 0.4646544
      P05     s_old 5  0.10 6.139350e-02 0.07524728 3.468219e-02 0.1203394
      P05     d_old 5  0.03 5.471818e-01 0.77269979 5.353819e-09 1.9475900
      P05     x_old 5  0.02 7.483299e-08 0.09130196 9.361255e-11 0.2490620
 median.absolute.log10.error within.twofold
                  0.06773351              5
                  0.15664213              5
                  0.25429874              3
                  0.21187758              4
                  1.65766896              0
                  5.42693689              0
```

## Simulation exclusions

```
     
      accepted no_old_crown no_tree
  P01        5            0       0
  P02        5            0       1
  P03        5            0       0
  P04        5            0       0
  P05        5            1       2
```

## Limits

Five trees per scenario support a diagnostic pilot, not precise bias or coverage estimates.
Recovery plots show individual estimates. No 95% confidence interval is claimed from these five replicates.
The scenario durations differ, so comparisons with the earlier taxon-stopped pilot are descriptive rather than paired experiments.
Convergence codes, agreement between starts, parameter recovery, and statistical identifiability are distinct checks.
Simulation cap exclusions and optimizer timeouts must be included in production reports, not silently discarded.

## Pilot interpretation and next steps

Recent speciation and dispersal generally recovered better than old dispersal and contraction/extinction.
P05 old dispersal and old contraction estimates were outside a factor of two of truth in all five datasets.
Increased old-epoch history did not by itself resolve recovery of these parameters.
Likelihood improvement over truth rules out simply failing to reach the generating likelihood; it does not prove the global optimum or validate the entire likelihood implementation.
Parameter spread can be large for rates near zero despite almost identical likelihoods; see the saved fit vectors.
Before a large run, profile old dispersal and contraction for representative P05 datasets and check likelihood trade-offs and numerical accuracy.
Follow with larger independent replicates spanning old-epoch information levels; report all convergence failures and boundary estimates.

## Reproduction

Run R/06-fixed-time-pilot.R, followed by R/07-report-fixed-time-pilot.R, from the DivExtra repository root.
datasets.rds contains the simulated trees and histories; fits.rds includes the initial vectors and full fit objects.

```
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS 27.0.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Europe/Helsinki
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] divextra_0.1.0     kableExtra_1.4.1   knitr_1.51         yaml_2.3.12       
[5] diversitree_0.10-1 ape_5.8-1          ggplot2_4.0.3     

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       dplyr_1.2.1        compiler_4.6.0     tidyselect_1.2.1  
 [5] Rcpp_1.1.1-1.1     xml2_1.5.2         stringr_1.6.0      parallel_4.6.0    
 [9] textshaping_1.0.5  systemfonts_1.3.2  scales_1.4.0       fastmap_1.2.0     
[13] lattice_0.22-9     R6_2.6.1           labeling_0.4.3     generics_0.1.4    
[17] subplex_1.9        tibble_3.3.1       svglite_2.2.2      pillar_1.11.1     
[21] RColorBrewer_1.1-3 rlang_1.2.0        stringi_1.8.7      deSolve_1.42      
[25] xfun_0.57          S7_0.2.2           otel_0.2.0         viridisLite_0.4.3 
[29] cli_3.6.6          withr_3.0.2        magrittr_2.0.5     digest_0.6.39     
[33] grid_4.6.0         rstudioapi_0.18.0  lifecycle_1.0.5    nlme_3.1-169      
[37] vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1         farver_2.1.2      
[41] ragg_1.5.2         codetools_0.2-20   rmarkdown_2.31     tools_4.6.0       
[45] pkgconfig_2.0.3    htmltools_0.5.9   
```

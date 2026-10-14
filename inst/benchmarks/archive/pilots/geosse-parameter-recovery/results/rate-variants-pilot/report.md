# Rate variants and classical GeoSSE control

Five new fixed-time trees per scenario. Each episodic and constant ClaSSE model has three starts.
The five constant control trees are additionally fitted using classical diversitree::make.geosse with sA=sB=sAB, dA=dB, xA=xB.
Parameters use ROOT.OBS, condition.surv=TRUE, complete sampling and tolerance 1e-10.
ClaSSE state order is A, B, AB; GeoSSE uses AB=0, A=1, B=2.
All models use the same nlminb optimizer, log rates, independently perturbed starting vectors and rate bounds [1e-10,100].
Three processes run concurrently; each start has a 120-second wall-time limit. Failed starts remain in starts.csv.
The P05-family durations are held fixed at old=30 and recent=6; constant-control durations are 25 and 6.
No selection on tip count or realized old-epoch information. Trees lacking a reconstructed crown across both epochs are excluded; all exclusions are logged.

## Design

```
   scenario_id s_old d_old x_old s_recent d_recent x_recent old.time
      baseline  0.10  0.03  0.02     0.30      0.2     0.15       30
 slow_recent_s  0.10  0.03  0.02     0.15      0.2     0.15       30
  higher_old_d  0.10  0.10  0.02     0.30      0.2     0.15       30
  higher_old_x  0.10  0.03  0.05     0.30      0.2     0.15       30
      constant  0.15  0.10  0.05     0.15      0.1     0.05       25
 recent.time
           6
           6
           6
           6
           6
```

## Convergence and model comparison

Converged starts: 155 of 165
Returned fits at the upper rate bound: 1

Incomplete or nonconverged starts:
```
             id           model start convergence.code
     baseline-2 constant_classe     2               NA
     baseline-3        episodic     1               NA
     baseline-4 constant_classe     2               NA
     baseline-5 constant_classe     2               NA
 higher_old_d-1 constant_classe     1               NA
 higher_old_d-3 constant_classe     1               NA
 higher_old_d-5 constant_classe     1               NA
 higher_old_x-1        episodic     3                1
 higher_old_x-4 constant_classe     1               NA
 higher_old_x-5        episodic     1               NA
                                                                                              message
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
                                                                                false convergence (8)
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
```
Nested likelihood violations (episodic below constant by >0.0001): 0
Best episodic fits below the generating-parameter likelihood by >0.0001: 0

```
      scenario median.N median.boundary.lineages episodic.converged.starts
      baseline       75                        8                        14
 slow_recent_s       49                       15                        15
  higher_old_d      454                       54                        15
  higher_old_x       40                        9                        13
      constant      133                       49                        15
 all.starts.agree episodic.AIC.support
                2                    4
                3                    0
                4                    5
                3                    0
                4                    0
```
Agreement requires all three starts to converge and a log-likelihood range <0.0001. This does not imply accurate parameter recovery.
convergence.csv also records maximum relative parameter spread across converged starts, with a 1e-6 denominator floor for rates near zero.
Episodic AIC support counts datasets where constant AIC minus episodic AIC exceeds 2.

## Classical GeoSSE equivalence

Maximum absolute matched-parameter log-likelihood difference: 3.13999e-06
Maximum best-fit log-likelihood difference: 1.18885e-08
Maximum absolute parameter-estimate difference: 7.2633e-10
Matched parameter tests include true, low, high and both fitted parameter vectors, with and without survival conditioning.
The GeoSSE and constrained ClaSSE fits describe the same model; they are not counted as separate alternatives in AIC selection.

## Recovery

```
      scenario parameter n truth       median        mean within.twofold
      baseline  s_recent 5  0.30 2.912820e-01  0.27714051              5
      baseline  d_recent 5  0.20 2.248964e-01  0.21069815              5
      baseline  x_recent 5  0.15 5.655506e-02  0.10767609              2
      baseline     s_old 5  0.10 7.093408e-02  0.05580556              3
      baseline     d_old 5  0.03 1.275339e-07 20.02748251              0
      baseline     x_old 5  0.02 2.452509e-05  0.51598001              1
 slow_recent_s  s_recent 5  0.15 1.611956e-01  0.16651397              5
 slow_recent_s  d_recent 5  0.20 2.581742e-01  0.21669913              4
 slow_recent_s  x_recent 5  0.15 2.607530e-01  0.31281187              3
 slow_recent_s     s_old 5  0.10 7.864704e-02  0.58273346              4
 slow_recent_s     d_old 5  0.03 8.963455e-02  0.11964169              0
 slow_recent_s     x_old 5  0.02 2.170549e-08  0.58259813              0
  higher_old_d  s_recent 5  0.30 3.064512e-01  0.29854756              5
  higher_old_d  d_recent 5  0.20 1.955279e-01  0.19780211              5
  higher_old_d  x_recent 5  0.15 1.487541e-01  0.14795056              5
  higher_old_d     s_old 5  0.10 1.085831e-01  0.10614065              4
  higher_old_d     d_old 5  0.10 1.040075e-01  0.25416822              2
  higher_old_d     x_old 5  0.02 3.757683e-07  0.03808141              0
  higher_old_x  s_recent 5  0.30 3.925318e-01  0.38579624              5
  higher_old_x  d_recent 5  0.20 2.258890e-01  0.29193200              3
  higher_old_x  x_recent 5  0.15 3.427919e-01  0.44024649              1
  higher_old_x     s_old 5  0.10 7.281164e-02  0.09425865              4
  higher_old_x     d_old 5  0.03 5.127635e-02  1.14344803              1
  higher_old_x     x_old 5  0.05 5.164790e-08  0.22240260              0
      constant  s_recent 5  0.15 1.251042e-01  0.12819599              5
      constant  d_recent 5  0.10 9.888139e-02  0.08435090              4
      constant  x_recent 5  0.05 2.353081e-02  0.03566027              2
      constant     s_old 5  0.15 1.609621e-01  0.13770952              4
      constant     d_old 5  0.10 1.325245e-01  0.43886131              4
      constant     x_old 5  0.05 8.639586e-02  0.09054006              1
 median.absolute.log10.error
                  0.06350118
                  0.05510689
                  0.42361976
                  0.14914508
                  5.37149561
                  2.91141933
                  0.05188366
                  0.12047548
                  0.24013812
                  0.24571892
                  0.96596920
                  5.96446032
                  0.01469450
                  0.01749623
                  0.04429403
                  0.04000545
                  0.31094214
                  4.72610989
                  0.11675364
                  0.09223851
                  0.71371662
                  0.20216938
                  2.23139285
                  5.98591737
                  0.07881929
                  0.10987148
                  0.32733308
                  0.05243486
                  0.17215389
                  0.56691818
```

## Interpretation of this pilot

The constant-rate implementation agrees closely with classical GeoSSE, both at matched parameter vectors and after optimization. This is strong evidence for the constant-rate special case, not proof of parameter identifiability or of every epoch-boundary calculation.
Recent speciation and dispersal are recovered much more consistently than old dispersal and contraction. Several old-rate estimates approach zero; one reaches the upper bound. Optimizer convergence alone therefore cannot establish reliable recovery.
Reducing recent speciation did not consistently resolve old-rate recovery. Increasing old dispersal produced the strongest recent-rate recovery, but also much larger trees (median 454 tips versus 75 in the baseline), so this cannot be attributed to the rate change alone.
Increasing old contraction did not resolve its recovery. Even in the higher-old-dispersal scenario, none of the five old-contraction estimates was within a factor of two of truth.
For the five constant controls, AIC preferred the constant model over the six-parameter episodic model. Five controls are too few to estimate a false-positive rate precisely.
For the next benchmark, retain both difficult and informative scenarios, increase replication, report recovery against realized tree size and old-epoch branch exposure, and extend/restart failed or disagreeing optimization runs. Do not treat these pilot results as confidence-interval coverage validation.

### Constant-model estimates (not the overparameterized episodic fits)

```
         id parameter truth       classe       geosse lnLik.difference
 constant-1         s  0.15 1.276539e-01 1.276539e-01     5.684910e-10
 constant-1         d  0.10 8.921911e-02 8.921911e-02     5.684910e-10
 constant-1         x  0.05 5.226535e-02 5.226535e-02     5.684910e-10
 constant-2         s  0.15 1.170870e-01 1.170870e-01     1.188852e-08
 constant-2         d  0.10 1.324545e-01 1.324545e-01     1.188852e-08
 constant-2         x  0.05 3.205655e-02 3.205655e-02     1.188852e-08
 constant-3         s  0.15 8.200082e-02 8.200082e-02    -1.858915e-09
 constant-3         d  0.10 1.262678e-01 1.262678e-01    -1.858915e-09
 constant-3         x  0.05 2.344628e-09 2.344628e-09    -1.858915e-09
 constant-4         s  0.15 1.457884e-01 1.457884e-01    -5.070433e-09
 constant-4         d  0.10 8.755087e-02 8.755087e-02    -5.070433e-09
 constant-4         x  0.05 4.543363e-02 4.543363e-02    -5.070433e-09
 constant-5         s  0.15 1.685640e-01 1.685640e-01    -8.709776e-09
 constant-5         d  0.10 1.201864e-01 1.201864e-01    -8.709776e-09
 constant-5         x  0.05 7.391618e-02 7.391618e-02    -8.709776e-09
```

## Figures

- [Parameter recovery](recovery.pdf): estimate divided by truth; the dashed line is exact recovery. Each point is one tree.
- [Classical GeoSSE comparison](geosse-comparison.pdf): matched constant-model parameter estimates; the dashed line is equality, not the generating truth.
- [Starting-point agreement](convergence.pdf): likelihood loss relative to the best returned start for each tree/model. Consult starts.csv for optimizer status; a finite likelihood is not necessarily converged.

## Limits

Five replicates per scenario are diagnostic only. Unequal realized tree sizes and old-epoch histories limit comparisons.
These are independent trees across scenarios, not paired common-history simulations. No confidence-interval coverage claim is made.
Changing old dispersal or contraction changes diversification and state occupancy as well as the number of transitions.
Slower recent speciation does not guarantee better recovery. Do not choose a final benchmark only because a pilot gives favorable estimates.
Only the constant control is checked against classical GeoSSE; this does not independently validate every episodic likelihood component.

## Reproduction

Run R/10-rate-variants-pilot.R and then R/11-report-rate-variants.R from the DivExtra repository root.
The runner reuses saved fits. To change the design use a separate results directory; do not reuse an old fits.rds.

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
[1] ggplot2_4.0.3      divextra_0.1.0     kableExtra_1.4.1   knitr_1.51        
[5] yaml_2.3.12        diversitree_0.10-1 ape_5.8-1         

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       dplyr_1.2.1        compiler_4.6.0     tidyselect_1.2.1  
 [5] Rcpp_1.1.1-1.1     xml2_1.5.2         stringr_1.6.0      parallel_4.6.0    
 [9] systemfonts_1.3.2  scales_1.4.0       textshaping_1.0.5  fastmap_1.2.0     
[13] lattice_0.22-9     R6_2.6.1           labeling_0.4.3     generics_0.1.4    
[17] subplex_1.9        tibble_3.3.1       svglite_2.2.2      pillar_1.11.1     
[21] RColorBrewer_1.1-3 rlang_1.2.0        deSolve_1.42       stringi_1.8.7     
[25] xfun_0.57          S7_0.2.2           otel_0.2.0         viridisLite_0.4.3 
[29] cli_3.6.6          withr_3.0.2        magrittr_2.0.5     digest_0.6.39     
[33] grid_4.6.0         rstudioapi_0.18.0  lifecycle_1.0.5    nlme_3.1-169      
[37] vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1         farver_2.1.2      
[41] ragg_1.5.2         codetools_0.2-20   rmarkdown_2.31     pkgconfig_2.0.3   
[45] tools_4.6.0        htmltools_0.5.9   
```

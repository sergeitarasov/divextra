# Calibrated epoch-grid recovery pilot

## Main findings

All 15 scenario-boundary combinations are represented by five qualifying trees. Of 270 optimization starts, 269 converged and one timed out; every tree has a converged fit. All three ClaSSE starts agree in likelihood for 71 of 75 trees. No returned estimate hits the upper parameter bound.
Recent speciation and dispersal medians are generally close to truth. Old contraction remains difficult, often approaching zero. In the all-rates-change scenario its median is near zero at all three boundaries. The information rules therefore do not guarantee recovery.
The cell-adjusted regressions suggest smaller recent-speciation errors with more tips, crossing lineages and recent branch exposure under both error metrics. Shared speciation improves with tip count, and shared dispersal with crossing lineages/branch exposure. Contraction slopes remain uncertain. These are exploratory patterns across 54 related fits, not multiplicity-adjusted discoveries.
At 50%, old speciation medians are elevated in the all-rates-change and speciation-only scenarios (about 0.169 and 0.168 versus truth 0.10). Qualifying trees are rare in these cells, so results must be interpreted conditional on selection. The pilot cannot separate selection effects from other sources of error.

## Procedure

Five independent trees per scenario-boundary cell, 15 cells, 75 trees. Boundaries are 5%, 25%, 50% of simulation duration measured backwards from present, NOT reconstructed crown age. Total durations are calibrated per cell using 30 independent attempts at each of 18,24,30,36,42,48 units; highest acceptance wins, ties favor median size nearest 200.
Simulation uses unchanged make.tree.classe.td(), uniformly sampled initial states A/B/AB, full sampling, and a 1000-taxon safety cap. Final trees require 50--1000 tips, at least 10 crossing lineages and 10 reconstructed branching events per epoch. All attempts are logged. Calibration never uses recovery quality; final trees use fresh seeds.
Inference imposes the scenario-specific links: 6/4/4/4/3 free parameters. Shared rates are listed once, not duplicated. ROOT.OBS, survival conditioning, ODE tolerance 1e-10, log-rate nlminb bounds 1e-10 to 100, three independently perturbed starts, 180-second limit per start.
Convergence requires optimizer code zero and a finite likelihood. All-start agreement requires three converged starts within 1e-4 log-likelihood. It is not proof of identifiable parameters.

## Design

s = speciation; d = range expansion; x = local extinction/range contraction.

~~~
          scenario_id s_old d_old x_old s_recent d_recent x_recent link_s
 all_rates_change_b05   0.1  0.03  0.02      0.3     0.20     0.05  FALSE
 all_rates_change_b25   0.1  0.03  0.02      0.3     0.20     0.05  FALSE
 all_rates_change_b50   0.1  0.03  0.02      0.3     0.20     0.05  FALSE
  speciation_only_b05   0.1  0.03  0.02      0.3     0.03     0.02  FALSE
  speciation_only_b25   0.1  0.03  0.02      0.3     0.03     0.02  FALSE
  speciation_only_b50   0.1  0.03  0.02      0.3     0.03     0.02  FALSE
   dispersal_only_b05   0.1  0.03  0.02      0.1     0.20     0.02   TRUE
   dispersal_only_b25   0.1  0.03  0.02      0.1     0.20     0.02   TRUE
   dispersal_only_b50   0.1  0.03  0.02      0.1     0.20     0.02   TRUE
 contraction_only_b05   0.1  0.03  0.02      0.1     0.03     0.05   TRUE
 contraction_only_b25   0.1  0.03  0.02      0.1     0.03     0.05   TRUE
 contraction_only_b50   0.1  0.03  0.02      0.1     0.03     0.05   TRUE
         constant_b05   0.1  0.03  0.02      0.1     0.03     0.02   TRUE
         constant_b25   0.1  0.03  0.02      0.1     0.03     0.02   TRUE
         constant_b50   0.1  0.03  0.02      0.1     0.03     0.02   TRUE
 link_d link_x n_parameters    rate.scenario fraction total.time old.time
  FALSE  FALSE            6 all_rates_change     0.05         48     45.6
  FALSE  FALSE            6 all_rates_change     0.25         30     22.5
  FALSE  FALSE            6 all_rates_change     0.50         18      9.0
   TRUE   TRUE            4  speciation_only     0.05         48     45.6
   TRUE   TRUE            4  speciation_only     0.25         36     27.0
   TRUE   TRUE            4  speciation_only     0.50         24     12.0
  FALSE   TRUE            4   dispersal_only     0.05         48     45.6
  FALSE   TRUE            4   dispersal_only     0.25         42     31.5
  FALSE   TRUE            4   dispersal_only     0.50         42     21.0
   TRUE  FALSE            4 contraction_only     0.05         48     45.6
   TRUE  FALSE            4 contraction_only     0.25         48     36.0
   TRUE  FALSE            4 contraction_only     0.50         48     24.0
   TRUE   TRUE            3         constant     0.05         48     45.6
   TRUE   TRUE            3         constant     0.25         42     31.5
   TRUE   TRUE            3         constant     0.50         42     21.0
 recent.time calibration.accepted
         2.4                   25
         7.5                   13
         9.0                    0
         2.4                   24
         9.0                   13
        12.0                    3
         2.4                   21
        10.5                   21
        21.0                   10
         2.4                   23
        12.0                   17
        24.0                   10
         2.4                   21
        10.5                   23
        21.0                   15
~~~

## Convergence

Converged starts: 269 / 270
Returned starts at upper parameter bound: 0

~~~
             scenario parameters median.N median.crossing converged.starts
 all_rates_change_b05          6      258             104               15
 all_rates_change_b25          6      191              14               15
 all_rates_change_b50          6      432              12               15
  speciation_only_b05          4      448             219               15
  speciation_only_b25          4      374              20               15
  speciation_only_b50          4      477              12               15
   dispersal_only_b05          4       82              66               15
   dispersal_only_b25          4      183              33               15
   dispersal_only_b50          4      689              19               15
 contraction_only_b05          4       80              63               15
 contraction_only_b25          4      153              50               15
 contraction_only_b50          4      177              22               14
         constant_b05          3      234             182               15
         constant_b25          3       85              26               15
         constant_b50          3      112              12               15
 datasets.with.fit all.three.agree attempts acceptance total.duration
                 5               5        6 0.83333333             48
                 5               5       14 0.35714286             30
                 5               4      500 0.01000000             18
                 5               5        6 0.83333333             48
                 5               5       15 0.33333333             36
                 5               5      106 0.04716981             24
                 5               4        8 0.62500000             48
                 5               5        5 1.00000000             42
                 5               4       27 0.18518519             42
                 5               5       11 0.45454545             48
                 5               5       12 0.41666667             48
                 5               4       14 0.35714286             48
                 5               5        7 0.71428571             48
                 5               5        6 0.83333333             42
                 5               5        9 0.55555556             42
~~~

Nonconverged starts:
~~~
                     id         model start
 contraction_only_b50-2 linked_classe     1
                                                                                              message
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
~~~

## Recovery

~~~
             scenario parameter n truth       median         mean
 all_rates_change_b05  s_recent 5  0.30 3.090846e-01 3.023860e-01
 all_rates_change_b05  d_recent 5  0.20 1.791584e-01 1.991439e-01
 all_rates_change_b05  x_recent 5  0.05 5.095952e-02 7.466282e-02
 all_rates_change_b05     s_old 5  0.10 8.683522e-02 8.962176e-02
 all_rates_change_b05     d_old 5  0.03 4.001963e-02 4.659255e-02
 all_rates_change_b05     x_old 5  0.02 1.902611e-07 3.732445e-03
 all_rates_change_b25  s_recent 5  0.30 2.994697e-01 2.928832e-01
 all_rates_change_b25  d_recent 5  0.20 2.105209e-01 2.119553e-01
 all_rates_change_b25  x_recent 5  0.05 8.409366e-02 6.042484e-02
 all_rates_change_b25     s_old 5  0.10 1.138759e-01 1.389256e-01
 all_rates_change_b25     d_old 5  0.03 1.730201e-08 1.319162e-02
 all_rates_change_b25     x_old 5  0.02 1.170835e-07 6.181784e-02
 all_rates_change_b50  s_recent 5  0.30 3.113004e-01 3.036776e-01
 all_rates_change_b50  d_recent 5  0.20 2.114597e-01 2.118574e-01
 all_rates_change_b50  x_recent 5  0.05 5.774976e-02 6.959700e-02
 all_rates_change_b50     s_old 5  0.10 1.692272e-01 1.668048e-01
 all_rates_change_b50     d_old 5  0.03 5.968660e-02 1.195386e-01
 all_rates_change_b50     x_old 5  0.02 1.769972e-07 2.738470e-07
  speciation_only_b05  s_recent 5  0.30 3.067127e-01 3.130578e-01
  speciation_only_b05     s_old 5  0.10 9.583051e-02 1.045642e-01
  speciation_only_b05  d_shared 5  0.03 3.044894e-02 2.984003e-02
  speciation_only_b05  x_shared 5  0.02 1.335203e-02 3.327430e-02
  speciation_only_b25  s_recent 5  0.30 3.030283e-01 3.219568e-01
  speciation_only_b25     s_old 5  0.10 1.091001e-01 1.319798e-01
  speciation_only_b25  d_shared 5  0.03 2.689170e-02 2.910619e-02
  speciation_only_b25  x_shared 5  0.02 2.933099e-02 5.780828e-02
  speciation_only_b50  s_recent 5  0.30 2.959263e-01 3.041600e-01
  speciation_only_b50     s_old 5  0.10 1.682802e-01 1.853575e-01
  speciation_only_b50  d_shared 5  0.03 3.659154e-02 3.451869e-02
  speciation_only_b50  x_shared 5  0.02 3.242752e-02 3.394450e-02
   dispersal_only_b05  d_recent 5  0.20 2.081625e-01 2.029217e-01
   dispersal_only_b05     d_old 5  0.03 2.479304e-02 2.420469e-02
   dispersal_only_b05  x_shared 5  0.02 1.145359e-02 1.975150e-02
   dispersal_only_b05  s_shared 5  0.10 9.217933e-02 9.852903e-02
   dispersal_only_b25  d_recent 5  0.20 2.027105e-01 2.364446e-01
   dispersal_only_b25     d_old 5  0.03 1.423950e-01 1.384644e-01
   dispersal_only_b25  x_shared 5  0.02 5.946813e-02 6.949648e-02
   dispersal_only_b25  s_shared 5  0.10 1.119589e-01 1.149390e-01
   dispersal_only_b50  d_recent 5  0.20 2.060431e-01 2.133732e-01
   dispersal_only_b50     d_old 5  0.03 4.154462e-02 3.643918e-02
   dispersal_only_b50  x_shared 5  0.02 3.196998e-02 2.543689e-02
   dispersal_only_b50  s_shared 5  0.10 1.016314e-01 9.950453e-02
 contraction_only_b05  x_recent 5  0.05 6.393694e-02 8.603573e-02
 contraction_only_b05     x_old 5  0.02 8.809052e-03 1.015366e-02
 contraction_only_b05  d_shared 5  0.03 3.075820e-02 2.762138e-02
 contraction_only_b05  s_shared 5  0.10 9.714085e-02 9.564098e-02
 contraction_only_b25  x_recent 5  0.05 3.760774e-02 4.274778e-02
 contraction_only_b25     x_old 5  0.02 3.978762e-08 3.242444e-03
 contraction_only_b25  d_shared 5  0.03 2.991672e-02 3.077122e-02
 contraction_only_b25  s_shared 5  0.10 9.399460e-02 9.488316e-02
 contraction_only_b50  x_recent 5  0.05 5.163742e-02 5.265114e-02
 contraction_only_b50     x_old 5  0.02 5.144908e-08 4.082348e-03
 contraction_only_b50  d_shared 5  0.03 3.259907e-02 3.437487e-02
 contraction_only_b50  s_shared 5  0.10 1.022279e-01 1.029349e-01
         constant_b05  d_shared 5  0.03 2.944442e-02 2.661152e-02
         constant_b05  x_shared 5  0.02 2.072995e-02 2.012259e-02
         constant_b05  s_shared 5  0.10 1.018392e-01 9.872072e-02
         constant_b25  d_shared 5  0.03 3.011438e-02 3.493160e-02
         constant_b25  x_shared 5  0.02 1.888367e-02 3.073735e-02
         constant_b25  s_shared 5  0.10 1.093947e-01 1.090437e-01
         constant_b50  d_shared 5  0.03 2.943821e-02 3.011151e-02
         constant_b50  x_shared 5  0.02 1.651303e-02 3.544557e-02
         constant_b50  s_shared 5  0.10 1.025010e-01 1.089987e-01
 within.twofold
              5
              5
              1
              5
              4
              1
              5
              5
              2
              4
              2
              0
              5
              5
              4
              4
              2
              0
              5
              5
              5
              2
              5
              5
              5
              2
              5
              3
              5
              2
              5
              5
              2
              5
              4
              0
              1
              5
              5
              1
              2
              5
              3
              2
              5
              5
              3
              1
              5
              5
              5
              1
              5
              5
              5
              4
              5
              5
              2
              5
              5
              4
              5
~~~

## Classical GeoSSE comparison

Maximum matched-parameter log-likelihood difference: 8.53527762956219e-06
Maximum fitted log-likelihood difference: 5.51062839804217e-09
Maximum fitted parameter difference: 1.12418040588286e-08
Maximum linked-versus-full ClaSSE likelihood difference: 4.54747350886464e-13
Linked fits below generating-parameter likelihood by >1e-4: 0

## Information and limitations

tree-information.csv records extant tip counts by state, crown age, boundary age, crossing lineages, branch exposure and reconstructed branching events by epoch. It separately records full-history simulated speciation, terminal extinction, dispersal and contraction events; simulation tables and transition histories are retained in datasets.rds.
The reconstructed root is counted as a branching event. A node exactly on the boundary belongs to the recent epoch.
Five replicates per cell cannot establish bias or coverage reliably. Sizes are outcomes, not controlled treatments. Plots of error against information are descriptive, not causal comparisons. No parameter-recovery confidence intervals or coverage are estimated in this pilot; regression slope intervals are reported separately.
The all-rates-change 50% cell had zero accepted trees in its 30-attempt calibration batches. The deterministic fallback selected the shortest candidate duration (18); the independent final sample required 500 attempts to obtain five qualifying trees. This very low acceptance must not be hidden when reporting recovery.
Models match the generating links; this is parameter recovery, not a test of model selection. No stochastic-mapping benchmark is performed here.

## Information regressions

Separate regressions for each independent parameter and each predictor: error ~ log2(information) + scenario-boundary cell fixed effects. Branch exposure is epoch-specific for epoch-specific rates and total exposure for shared rates. This controls duration and rate settings through cell effects. Crossing counts are not fitted jointly with tip counts or branch exposure because they are correlated.
Outcomes: absolute log10 estimate/truth error, and bounded symmetric error abs(estimate-truth)/(estimate+truth). Negative slopes mean better recovery with more information. HC3 robust 95% intervals are exploratory, not adjusted for multiple comparisons. Sparse cells, selected trees, boundary estimates and failed fits limit inference; no causal conclusion.

~~~
 parameter         predictor              outcome  n cells slope.per.doubling
  s_recent                 N absolute.log10.error 30     6       -0.026335879
  s_recent                 N        bounded.error 30     6       -0.030152334
  s_recent boundary.lineages absolute.log10.error 30     6       -0.020641365
  s_recent boundary.lineages        bounded.error 30     6       -0.023638152
  s_recent   branch.exposure absolute.log10.error 30     6       -0.024833451
  s_recent   branch.exposure        bounded.error 30     6       -0.028438357
  d_recent                 N absolute.log10.error 30     6       -0.009672641
  d_recent                 N        bounded.error 30     6       -0.010863952
  d_recent boundary.lineages absolute.log10.error 30     6       -0.001739228
  d_recent boundary.lineages        bounded.error 30     6       -0.001985995
  d_recent   branch.exposure absolute.log10.error 30     6       -0.010324994
  d_recent   branch.exposure        bounded.error 30     6       -0.011555703
  x_recent                 N absolute.log10.error 30     6       -0.422478391
  x_recent                 N        bounded.error 30     6       -0.093233981
  x_recent boundary.lineages absolute.log10.error 30     6       -0.400531252
  x_recent boundary.lineages        bounded.error 30     6       -0.093014840
  x_recent   branch.exposure absolute.log10.error 30     6       -0.431374372
  x_recent   branch.exposure        bounded.error 30     6       -0.095118112
     s_old                 N absolute.log10.error 30     6       -0.025868012
     s_old                 N        bounded.error 30     6       -0.030039778
     s_old boundary.lineages absolute.log10.error 30     6       -0.029589677
     s_old boundary.lineages        bounded.error 30     6       -0.033173371
     s_old   branch.exposure absolute.log10.error 30     6       -0.045689887
     s_old   branch.exposure        bounded.error 30     6       -0.050241044
     d_old                 N absolute.log10.error 30     6       -0.533458273
     d_old                 N        bounded.error 30     6       -0.116189770
     d_old boundary.lineages absolute.log10.error 30     6       -0.463260804
     d_old boundary.lineages        bounded.error 30     6       -0.113074718
     d_old   branch.exposure absolute.log10.error 30     6       -0.080016775
     d_old   branch.exposure        bounded.error 30     6       -0.066613078
     x_old                 N absolute.log10.error 30     6       -0.633765791
     x_old                 N        bounded.error 30     6       -0.060735927
     x_old boundary.lineages absolute.log10.error 30     6       -0.622172916
     x_old boundary.lineages        bounded.error 30     6       -0.054896906
     x_old   branch.exposure absolute.log10.error 30     6       -0.549740379
     x_old   branch.exposure        bounded.error 30     6       -0.039013911
  s_shared                 N absolute.log10.error 45     9       -0.012794749
  s_shared                 N        bounded.error 45     9       -0.014626972
  s_shared boundary.lineages absolute.log10.error 45     9       -0.011027281
  s_shared boundary.lineages        bounded.error 45     9       -0.012613753
  s_shared   branch.exposure absolute.log10.error 45     9       -0.012324780
  s_shared   branch.exposure        bounded.error 45     9       -0.014084196
  d_shared                 N absolute.log10.error 45     9       -0.030671844
  d_shared                 N        bounded.error 45     9       -0.034356906
  d_shared boundary.lineages absolute.log10.error 45     9       -0.033560119
  d_shared boundary.lineages        bounded.error 45     9       -0.037587015
  d_shared   branch.exposure absolute.log10.error 45     9       -0.031847991
  d_shared   branch.exposure        bounded.error 45     9       -0.035680196
  x_shared                 N absolute.log10.error 45     9       -0.140827604
  x_shared                 N        bounded.error 45     9       -0.107383853
  x_shared boundary.lineages absolute.log10.error 45     9       -0.124125837
  x_shared boundary.lineages        bounded.error 45     9       -0.100849512
  x_shared   branch.exposure absolute.log10.error 45     9       -0.155188146
  x_shared   branch.exposure        bounded.error 45     9       -0.104596727
       lower         upper
 -0.04500795 -0.0076638105
 -0.05155093 -0.0087537385
 -0.03529285 -0.0059898770
 -0.04043667 -0.0068396317
 -0.04214857 -0.0075183328
 -0.04829033 -0.0085863805
 -0.03456903  0.0152237532
 -0.03900417  0.0172762601
 -0.02835023  0.0248717734
 -0.03210327  0.0281312798
 -0.03624932  0.0155993295
 -0.04078716  0.0176757503
 -1.13309728  0.2881404959
 -0.19621990  0.0097519335
 -1.11365078  0.3125882731
 -0.20219781  0.0161681315
 -1.16856145  0.3058127050
 -0.20005512  0.0098188951
 -0.06952705  0.0177910236
 -0.07750089  0.0174213366
 -0.06381442  0.0046350645
 -0.07036217  0.0040154255
 -0.09231093  0.0009311619
 -0.09963735 -0.0008447374
 -1.60687571  0.5399591617
 -0.27363099  0.0412514480
 -1.42394585  0.4974242391
 -0.27524037  0.0490909315
 -0.76438999  0.6043564424
 -0.19123035  0.0580041918
 -2.07750486  0.8099732756
 -0.23567990  0.1142080418
 -2.21792455  0.9735787220
 -0.25132086  0.1415270430
 -2.19985777  1.1003770089
 -0.25104836  0.1730205350
 -0.02467067 -0.0009188329
 -0.02820037 -0.0010535681
 -0.02337330  0.0013187391
 -0.02673616  0.0015086511
 -0.02517704  0.0005274776
 -0.02877096  0.0006025721
 -0.06150897  0.0001652795
 -0.06897274  0.0002589258
 -0.06484533 -0.0022749121
 -0.07271132 -0.0024627128
 -0.06349761 -0.0001983681
 -0.07121979 -0.0001405987
 -0.92704171  0.6453865022
 -0.22066826  0.0059005499
 -1.02055948  0.7723078026
 -0.22217696  0.0204779336
 -0.99096898  0.6805926897
 -0.22399633  0.0148028737
~~~

## Figures

- [Information regressions](information-regressions.pdf)
- [Tree sizes and crossing lineages](tree-information.pdf)
- [Recovery](recovery.pdf)
- [Recovery versus tip number](information-recovery.pdf)
- [Old-rate recovery versus crossing lineages](old-information.pdf)
- [GeoSSE comparison](geosse-comparison.pdf)

## Reproduce

Run R/15-epoch-grid-pilot.R followed by R/16-report-epoch-grid-pilot.R from the DivExtra repository root. Use a fresh output directory if changing settings; existing fits are reused.

~~~
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
~~~

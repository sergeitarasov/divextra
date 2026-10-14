# Linked-rate recovery pilot

Five independent trees per scenario; three optimization starts per fit. All scenarios use 30 old + 6 recent time units. The boundary is fixed at age 6; this pilot does not vary boundary position.
Simulation uses unchanged make.tree.classe.td(), uniformly sampled initial states A/B/AB, full sampling, and a 4000-taxon safety cap. Trees that stop early or lack a crown spanning both epochs are excluded; all attempts are logged.
Inference imposes the scenario-specific links: 6/4/4/4/3 free parameters. Shared rates are listed once, not duplicated. ROOT.OBS, survival conditioning, ODE tolerance 1e-10, log-rate nlminb bounds 1e-10 to 100, three independently perturbed starts, 180-second limit per start.
Convergence requires optimizer code zero and a finite likelihood. All-start agreement requires three converged starts within 1e-4 log-likelihood. It is not proof of identifiable parameters.

## Design

s = speciation; d = range expansion; x = local extinction/range contraction.

~~~
      scenario_id s_old d_old x_old s_recent d_recent x_recent link_s link_d
 all_rates_change   0.1  0.03  0.02      0.3     0.20     0.05  FALSE  FALSE
  speciation_only   0.1  0.03  0.02      0.3     0.03     0.02  FALSE   TRUE
   dispersal_only   0.1  0.03  0.02      0.1     0.20     0.02   TRUE  FALSE
 contraction_only   0.1  0.03  0.02      0.1     0.03     0.05   TRUE   TRUE
         constant   0.1  0.03  0.02      0.1     0.03     0.02   TRUE   TRUE
 link_x n_parameters old.time recent.time
  FALSE            6       30           6
   TRUE            4       30           6
   TRUE            4       30           6
  FALSE            4       30           6
   TRUE            3       30           6
~~~

## Convergence

Converged starts: 87 / 90

~~~
         scenario parameters median.N median.crossing converged.starts
 all_rates_change          6      152              14               12
  speciation_only          4      136              19               15
   dispersal_only          4      135              47               15
 contraction_only          4       27              16               15
         constant          3       71              36               15
 datasets.with.fit all.three.agree
                 4               4
                 5               5
                 5               5
                 5               5
                 5               5
~~~

Nonconverged starts:
~~~
                 id         model start               message
 all_rates_change-1 linked_classe     1 false convergence (8)
 all_rates_change-1 linked_classe     2 false convergence (8)
 all_rates_change-1 linked_classe     3 false convergence (8)
~~~

## Recovery

~~~
         scenario parameter n truth       median         mean within.twofold
 all_rates_change  s_recent 4  0.30 2.899183e-01  0.289112495              4
 all_rates_change  d_recent 4  0.20 2.172539e-01  0.212094932              4
 all_rates_change  x_recent 4  0.05 7.679122e-02  0.078883206              2
 all_rates_change     s_old 4  0.10 9.887153e-02  0.101735643              4
 all_rates_change     d_old 4  0.03 1.323737e+00 25.667575310              2
 all_rates_change     x_old 4  0.02 3.587971e-01  1.094836547              0
  speciation_only  s_recent 5  0.30 2.985201e-01  0.299809950              5
  speciation_only     s_old 5  0.10 1.038898e-01  0.104324986              5
  speciation_only  d_shared 5  0.03 2.482380e-02  0.030169522              5
  speciation_only  x_shared 5  0.02 3.169679e-08  0.015898702              1
   dispersal_only  d_recent 5  0.20 1.964087e-01  0.193416462              5
   dispersal_only     d_old 5  0.03 3.235419e-02  0.038802900              3
   dispersal_only  x_shared 5  0.02 3.231783e-08  0.004616498              1
   dispersal_only  s_shared 5  0.10 1.043882e-01  0.105599547              5
 contraction_only  x_recent 5  0.05 1.596013e-08  0.012734056              1
 contraction_only     x_old 5  0.02 1.262247e-10  0.003613764              1
 contraction_only  d_shared 5  0.03 2.743068e-02  0.029141893              5
 contraction_only  s_shared 5  0.10 8.320045e-02  0.079266924              4
         constant  d_shared 5  0.03 2.887807e-02  0.035358065              5
         constant  x_shared 5  0.02 1.392564e-02  0.032421318              1
         constant  s_shared 5  0.10 1.159001e-01  0.115583388              5
~~~

## Classical GeoSSE comparison

Maximum matched-parameter log-likelihood difference: 1.40363286504908e-06
Maximum fitted log-likelihood difference: 1.71701231010957e-09
Maximum fitted parameter difference: 8.05515751323149e-10
Maximum linked-versus-full ClaSSE likelihood difference: 0
Linked fits below generating-parameter likelihood by >1e-4: 0

## Information and limitations

tree-information.csv records extant tip counts by state, crown age, boundary age, crossing lineages, branch exposure and reconstructed branching events by epoch. It separately records full-history simulated speciation, terminal extinction, dispersal and contraction events; simulation tables and transition histories are retained in datasets.rds.
The reconstructed root is counted as a branching event. A node exactly on the boundary belongs to the recent epoch.
Five replicates per scenario cannot establish bias or coverage reliably. Sizes are outcomes, not controlled treatments. Plots of error against information are descriptive, not causal comparisons. No confidence intervals or coverage are estimated in this pilot.
Models match the generating links; this is parameter recovery, not a test of model selection. No stochastic-mapping benchmark is performed here.

## Figures

- [Recovery](recovery.pdf)
- [Recovery versus tip number](information-recovery.pdf)
- [Old-rate recovery versus crossing lineages](old-information.pdf)
- [GeoSSE comparison](geosse-comparison.pdf)

## Reproduce

Run R/12-linked-rate-pilot.R followed by R/13-report-linked-rate-pilot.R from the DivExtra repository root. Use a fresh output directory if changing settings; existing fits are reused.

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

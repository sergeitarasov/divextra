# Tree-information screening

50 fresh attempts per cell: five rate scenarios times three boundary positions. Initial state sampled uniformly. No inference.
Total stem simulation duration is 36; boundary ages are 9, 18, 27 before present. These are NOT fixed fractions of reconstructed crown age.
Acceptance: at least 50 extant tips, 10 reconstructed crossing lineages, and 10 reconstructed branching events per epoch.
For a binary reconstructed tree spanning the boundary, old.nodes = crossing - 1 and recent.nodes = tips - crossing. Therefore the old-node rule effectively requires at least 11 crossing lineages.
All distributions below are conditional on acceptance and avoiding the 4000-taxon safety cap. Capped runs are not small trees or extinction failures: they are censored and their final tip counts are unknown.
Ranges and quantiles describe this screening sample, not confidence intervals. Small or zero retained samples cannot characterize a distribution.

~~~
         scenario boundary.fraction accepted capped N.median N.min N.max
 all_rates_change              0.25       28      0    727.5   194  2340
 all_rates_change              0.50        0     33       NA    NA    NA
 all_rates_change              0.75        0     30       NA    NA    NA
  speciation_only              0.25       22      0    407.5   106  1545
  speciation_only              0.50        7      9   3052.0  2243  3707
  speciation_only              0.75        0     28       NA    NA    NA
   dispersal_only              0.25       25      0    116.0    56   516
   dispersal_only              0.50       10      0    288.0   174   737
   dispersal_only              0.75        1      0   1178.0  1178  1178
 contraction_only              0.25       12      0     76.0    50   103
 contraction_only              0.50       10      0     84.0    64   184
 contraction_only              0.75        0      0       NA    NA    NA
         constant              0.25       26      0     81.5    51   334
         constant              0.50       14      0    109.0    59   285
         constant              0.75        2      0    342.5   276   409
 crossing.median crossing.min crossing.max
            24.0           11           69
              NA           NA           NA
              NA           NA           NA
            30.5           11           92
            12.0           11           18
              NA           NA           NA
            27.0           11          114
            14.0           11           31
            11.0           11           11
            29.0           17           48
            15.5           13           32
              NA           NA           NA
            32.0           17          108
            16.0           11           39
            13.0           13           13
~~~

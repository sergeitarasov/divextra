# Tip-state representation and recovery

## Main findings

All 75 trees already contain A, B and AB, and 72 meet the five-tips-per-state rule. Missing-state trees cannot explain the existing recovery failures.
For old contraction, 21/30 estimates are below 1e-5. Applying the five-tip rule leaves 29 trees and all 21 near-zero estimates; recovery within a factor of two changes only from 5/30 to 5/29. Thus this rule does not resolve the problem in the existing sample.
After controlling tree size and scenario-boundary cell, old-contraction association intervals include zero for both state-balance predictors and both error definitions. The sample does not establish absence of an effect, but supplies no clear evidence of improved old-contraction recovery with more balanced present-day states.
Across all 36 exploratory regressions, two nominal intervals exclude zero, both in the opposite direction (more AB associated with worse log-scale error for recent/shared contraction). Neither is reproduced with bounded error. These are not multiplicity-adjusted findings and do not support choosing a state threshold.

## Procedure

Existing 75-tree epoch-grid pilot only. No new simulations, rejections or inference fits. State counts were verified to sum to tip number.
All three states present: 75 of 75
At least five tips per state: 72 of 75

## Median state counts by cell

These are separate column medians and need not sum to median tree size.
~~~
             scenario n.A n.B n.AB min.count all.present five.each
 all_rates_change_b05  77 106   44        44           5         5
 all_rates_change_b25  75  76   45        39           5         5
 all_rates_change_b50 178 164   92        92           5         5
         constant_b05 105 109   20        20           5         5
         constant_b25  36  42   10        10           5         4
         constant_b50  56  72   13        13           5         5
 contraction_only_b05  53  32    8         8           5         4
 contraction_only_b25  65  74   21        21           5         5
 contraction_only_b50  91  70   22        22           5         5
   dispersal_only_b05  32  17   33        15           5         5
   dispersal_only_b25  65  54   69        54           5         5
   dispersal_only_b50 221 205  263       205           5         5
  speciation_only_b05 269 143   36        36           5         4
  speciation_only_b25 242 116   15        15           5         5
  speciation_only_b50 153 176   25        25           5         5
~~~

## Trees failing five tips per state

~~~
                     id   N n.A n.B n.AB
  speciation_only_b05-5 101  57  43    1
 contraction_only_b05-3  64  37  23    4
         constant_b25-5  72  36  33    3
~~~

## Recovery with and without the proposed five-tip rule

Summaries are by independently estimated parameter. The two groups overlap; they are not randomized treatments.
~~~
 parameter      rule  n within.twofold near.zero median.bounded.error
  s_recent       all 30             30         0           0.03380222
  s_recent five.each 29             29         0           0.02941415
  d_recent       all 30             29         0           0.04832418
  d_recent five.each 30             29         0           0.04832418
  x_recent       all 30             18         3           0.24481044
  x_recent five.each 29             18         3           0.23536828
     s_old       all 30             26         0           0.11976754
     s_old five.each 29             25         0           0.11718793
     d_old       all 30             14         6           0.36160840
     d_old five.each 30             14         6           0.36160840
     x_old       all 30              5        21           0.99998256
     x_old five.each 29              5        21           0.99998281
  s_shared       all 45             45         0           0.03644342
  s_shared five.each 43             43         0           0.03095653
  d_shared       all 45             45         0           0.06562574
  d_shared five.each 42             42         0           0.06448667
  x_shared       all 45             21         7           0.35464657
  x_shared five.each 43             21         7           0.33819941
~~~

## Continuous association checks

Separate models for each parameter and state predictor: error ~ state proportion + log2(tip number) + scenario-boundary cell. Slopes are per 10 percentage-point increase; HC3 95% intervals. Outcomes are absolute log10 estimate/truth error and bounded symmetric error. Negative slopes mean better recovery.
No missing-state trees exist here, and only three trees fail the five-tip rule. This dataset cannot estimate the benefit of requiring all states or give a reliable threshold-treatment comparison. Continuous results are exploratory, with no multiple-testing adjustment. Five trees per cell and correlated predictors limit precision.
~~~
 parameter    predictor              outcome  n slope.per.ten.percentage.points
  s_recent min.fraction absolute.log10.error 30                     0.056793860
  s_recent min.fraction        bounded.error 30                     0.065086170
  s_recent  AB.fraction absolute.log10.error 30                     0.069325026
  s_recent  AB.fraction        bounded.error 30                     0.079181718
  d_recent min.fraction absolute.log10.error 30                    -0.067820116
  d_recent min.fraction        bounded.error 30                    -0.075409399
  d_recent  AB.fraction absolute.log10.error 30                     0.037900480
  d_recent  AB.fraction        bounded.error 30                     0.042930221
  x_recent min.fraction absolute.log10.error 30                     2.195329425
  x_recent min.fraction        bounded.error 30                     0.137816650
  x_recent  AB.fraction absolute.log10.error 30                     2.257111228
  x_recent  AB.fraction        bounded.error 30                     0.268401441
     s_old min.fraction absolute.log10.error 30                     0.008792606
     s_old min.fraction        bounded.error 30                     0.011748829
     s_old  AB.fraction absolute.log10.error 30                    -0.044010375
     s_old  AB.fraction        bounded.error 30                    -0.045736712
     d_old min.fraction absolute.log10.error 30                     2.250607978
     d_old min.fraction        bounded.error 30                     0.316955289
     d_old  AB.fraction absolute.log10.error 30                     1.115705420
     d_old  AB.fraction        bounded.error 30                    -0.043785773
     x_old min.fraction absolute.log10.error 30                    -4.035771697
     x_old min.fraction        bounded.error 30                    -0.269460504
     x_old  AB.fraction absolute.log10.error 30                    -1.706998921
     x_old  AB.fraction        bounded.error 30                    -0.106257550
  s_shared min.fraction absolute.log10.error 45                    -0.013489071
  s_shared min.fraction        bounded.error 45                    -0.015389436
  s_shared  AB.fraction absolute.log10.error 45                    -0.028323034
  s_shared  AB.fraction        bounded.error 45                    -0.032395142
  d_shared min.fraction absolute.log10.error 45                     0.094754259
  d_shared min.fraction        bounded.error 45                     0.107252173
  d_shared  AB.fraction absolute.log10.error 45                     0.094754259
  d_shared  AB.fraction        bounded.error 45                     0.107252173
  x_shared min.fraction absolute.log10.error 45                     0.945946985
  x_shared min.fraction        bounded.error 45                     0.238468557
  x_shared  AB.fraction absolute.log10.error 45                     3.818093497
  x_shared  AB.fraction        bounded.error 45                     0.311981699
        lower       upper
  -0.02922663 0.142814349
  -0.03326763 0.163439971
  -0.01209230 0.150742353
  -0.01330109 0.171664530
  -0.21399105 0.078350819
  -0.23743856 0.086619764
  -0.02384403 0.099644989
  -0.02659564 0.112456081
  -2.04831045 6.438969298
  -0.63285698 0.908490275
   0.21304629 4.301176165
  -0.27940329 0.816206174
  -0.23828502 0.255870234
  -0.25270112 0.276198777
  -0.25058958 0.162568833
  -0.26756759 0.176094164
  -1.77885419 6.280070145
  -0.09456851 0.728479087
  -3.02355964 5.254970479
  -0.48117755 0.393606007
 -12.56994024 4.498396846
  -1.53658378 0.997662777
  -7.82319414 4.409196304
  -0.69146532 0.478950216
  -0.05458925 0.027611105
  -0.06238003 0.031601154
  -0.06430410 0.007658036
  -0.07355195 0.008761669
  -0.02624488 0.215753403
  -0.02781526 0.242319609
  -0.02624488 0.215753403
  -0.02781526 0.242319609
  -4.48013982 6.372033789
  -0.45617724 0.933114355
   1.14541898 6.490768018
  -0.11618766 0.740151054
~~~

[Old contraction versus smallest state count](old-contraction-tip-states.pdf)

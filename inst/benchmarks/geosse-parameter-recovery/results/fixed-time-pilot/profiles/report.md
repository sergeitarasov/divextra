# P05 likelihood profiles and numerical audit

The evidence favors weak practical identifiability over integration error in these three datasets.
For all six profiles the true focal rate is less than 1.92 log-likelihood units below the best value, after optimizing nuisance rates.
P05-3 dispersal stays within 1.45 units across the full evaluated range from approximately 1e-8 to 10.
Some extreme tails timed out or encountered a nuisance bound and remain unresolved. No formal confidence interval is claimed.

Representative datasets: P05-1 (499 tips), P05-3 (216 tips), P05-4 (68 tips).
Each grid value fixes one old-epoch rate while fitting the other five rates.
Two nuisance starts are used: the previous best fit and the true parameters.
Optimization uses log rates, nlminb, and GSL rkck tolerance 1e-10. Free rates are bounded between 1e-10 and 100.
Both optimizer objects and nuisance vectors are saved. Each grid task has a 120-second process limit.

Completed grid tasks: 52 of 57
Nonzero chosen optimizer codes: 0
Tasks hitting the nuisance upper bound: 1

Failed or timed-out points (not interpreted as evidence against their rate values):
```
    id  rate value
 P05-1 mu1.2   2.0
 P05-3 mu1.2   0.6
 P05-4 q13.2  10.0
 P05-4 mu1.2   0.6
 P05-4 mu1.2   2.0
                                                                                                error
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
 ! callr subprocess failed: could not start R, exited with non-zero status, has crashed or was killed
```

## Profiles

```
    id  rate best.grid.value true.value loss.at.truth   grid.lower grid.upper
 P05-1 q13.2    5.471818e-01       0.03    1.55628642 1.000000e-02      10.00
 P05-1 mu1.2    2.074477e-01       0.02    1.01358683 1.000000e-08       0.60
 P05-3 q13.2    1.363924e+00       0.03    1.16771823 1.000000e-08      10.00
 P05-3 mu1.2    1.000000e-08       0.02    0.05260663 1.000000e-08       0.20
 P05-4 q13.2    5.353819e-09       0.03    0.76317957 5.353819e-09       0.03
 P05-4 mu1.2    9.361255e-11       0.02    0.64347894 9.361255e-11       0.06
 lower.edge.supported upper.edge.supported
                FALSE                 TRUE
                 TRUE                 TRUE
                 TRUE                 TRUE
                 TRUE                FALSE
                 TRUE                FALSE
                 TRUE                FALSE
```
The reported grid limits are sampled values with likelihood loss <=1.92. They are not interpolated confidence limits.
A supported grid edge means that direction is unresolved within the tested range.
Nuisance optima at the upper bound 100 are constrained by the chosen search range; those points cannot establish an unconstrained tail endpoint.
The 1.92 threshold is a conventional one-parameter likelihood reference. Boundary parameters and this sampling design require care; coverage was not calibrated.

## Numerical checks

Maximum likelihood range across rkck tolerances 1e-8, 1e-10, 1e-12: 7.21559e-06
Maximum absolute difference between rk8pd and rkck at 1e-12: 8.89941e-10
The deSolve backend was attempted but returns "Not yet possible" for this likelihood; it does not provide an independent validation here.
Agreement between GSL methods and tolerances tests integration stability. It does not rule out a shared equation or implementation error.

## Slower rates

```
 time.multiplier rate.multiplier likelihood.contrast
               1             1.0            5.729677
               2             0.5            5.729677
               5             0.2            5.729677
```
Uniformly slowing every rate while lengthening all times by the reciprocal factor preserves likelihood contrasts.
Slowing every rate at fixed duration usually produces fewer events and is not a general remedy for poor recovery.
A useful next design changes relative rates: slower RECENT speciation to reduce recent dominance, with an independently calibrated old interval.
For a controlled follow-up, compare recent s=0.15 versus 0.30 while keeping old s=0.10 and the other P05 rates fixed; both epochs still differ in all three rates.
Also vary old d independently (e.g. 0.03, 0.10) and old x (e.g. 0.02, 0.05) to test whether rare old events drive the weakness.
Changing these rates is a proposed follow-up, not an analysis performed in this report.

## Interpretation boundaries

Broad nuisance-optimized profiles with stable numerical likelihoods support weak practical identifiability for these observed datasets.
They do not prove structural non-identifiability, global optimality, or validity of every implementation component.
Only three pilot datasets are profiled. No general recovery guarantee follows from these cases.

## Reproduction

Run R/08-profile-p05.R then R/09-report-p05-profiles.R from the DivExtra repository root.

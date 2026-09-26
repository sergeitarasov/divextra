# Direct extinction sampler: reliability audit

## Conclusion

The event sampler has strong mathematical and simulation support for its stated
target: a side clade with a fixed founding age and state, conditional on complete
extinction before present. At ordinary test rates, event times, extinct-tip
counts, branch lengths, state occupancy, and transitions agree with independent
answers. This is empirical validation, not proof of all numerical cases.

Two limitations remain:

1. The complete mapping procedure is ASR-weighted augmentation, not a posterior
   conditional on the observed reconstructed tree. Conditioning each attached
   clade on extinction does not condition the birth locations on the backbone.
2. The fixed absolute interpolation step can be inaccurate at very fast rates.
   Rate-aware/adaptive refinement should precede claims of general numerical
   reliability. No production simulator changes were made during this audit.

Reproduce with the current package loaded:

```r
devtools::load_all("/Users/taravser/GitHub/divextra")
source("/Users/taravser/GitHub/divextra/inst/benchmarks/validate-direct-extinction.R")
```

## Mathematical target and derivation

Let E_i(t) be extinction probability for a lineage in state i at age t before
present, and r_i(t) its original total event rate. E_i(0)=0. For unordered
daughter pairs, the backward equation is

    E_i' = mu_i - r_i E_i + sum_j q_ij E_j
           + sum_{j<=k} lambda_ijk E_j E_k.

There is no extra factor of two in the unordered-pair birth sum. For a lineage
starting at age a, the probability of no event until younger age u, conditional
on extinction, is

    S_i(u | a, extinct) = exp(-integral_u^a r_i(v) dv) E_i(u)/E_i(a).

The code inverts this expression. At the sampled event age, event weights are
mu_i, q_ij E_j, and lambda_ijk E_j E_k. These are the correct conditional hazards
after canceling their common denominator E_i. After birth the two daughter
clades are independent conditional on both becoming extinct. E is continuous
across epochs but the original rates change at their absolute-age boundaries.
The numerical implementation approximates E by linear interpolation.

## Analytical one-state checks

For constant lambda and mu, founder age T=3:

    E(t) = mu (1-exp(-(lambda-mu)t)) /
           (lambda-mu exp(-(lambda-mu)t)),   lambda != mu;
    E(t) = lambda t/(1+lambda t),            lambda == mu.

The conditional CDF of time to extinction of the entire clade is E(s)/E(T).
For lambda=0, the founder's lifetime is a truncated exponential.

Writing p0=E(s), beta=lambda*p0/mu and z=E(T-s), the mean number of lineages
alive at elapsed time s, conditional on extinction by T, is

    m(s) = (1-p0)(1-beta) z / ((1-beta*z)^2 E(T)).

Expected total branch length is integral m(s) ds. Expected births are
integral lambda E(T-s) m(s) ds; expected extinct tips are one plus this value.
The probability of no births is

    [mu/(lambda+mu)] [1-exp(-(lambda+mu)T)] / E(T).

The audit sampled 2,000 independent side trees for each of pure death,
subcritical, supercritical, and critical cases. It checked all these quantities,
not merely the absence of surviving tips.

| Rates (lambda, mu) | Mean tips simulated | Mathematical mean |
|---|---:|---:|
| (0, 0.7) | 1.0000 | 1.0000 |
| (0.3, 0.7) | 1.2330 | 1.2297 |
| (0.7, 0.3) | 1.2400 | 1.2297 |
| (0.5, 0.5), initial 2,000 | 1.2530 | 1.3000 |
| (0.5, 0.5), independent 20,000 | 1.30315 | 1.3000 |

The initial critical result was 3.43 estimated Monte Carlo standard errors below
the expectation and was retained, not discarded. An independent 20,000-tree
check gave tips 1.30315 +/- 0.00490 MC SE, branch length 1.20446 versus 1.2,
midpoint richness 0.35065 versus 0.35, and no-birth probability 0.78905 versus
0.79184. All four means were within one MC SE of their analytical expectations.
Its extinction-time probability-integral-transform ECDF distance was 0.00450.

## Independent multistate, epoch-specific oracle

A separate ordinary forward Gillespie simulator was implemented in the audit,
without calling the package's event-rate parser, extinction cache, or event
generator. It accepts only trees that become extinct. The two-state model has
anagenetic changes and both same-state and different-state daughter pairs, with
a rate change at age 1 and founder age 3.

Compared 2,000 direct trees with 2,000 accepted reference trees (2,373 attempts).

| Summary | Direct | Independent rejection |
|---|---:|---:|
| Extinct tips | 1.1950 | 1.1855 |
| Anagenetic transitions | 0.1060 | 0.1025 |
| Total branch length | 0.89371 | 0.87946 |
| Midpoint richness | 0.2005 | 0.2005 |
| Midpoint richness in state 2 | 0.0515 | 0.0520 |
| Extinct tips in state 2 | 0.1040 | 0.1035 |
| Elapsed time to final extinction | 0.80707 | 0.80227 |

All nine compared summaries, including births and no-birth probability, differed
by less than 0.59 combined Monte Carlo SE. This is agreement at the resolution
of these experiments, not a guarantee of equality in every distributional tail.

## Important backbone counterexample

The public function was tested on a two-tip tree with both branches length 2,
one state, lambda=0.5, and mu=0.05. Across 1,000 maps the mean hidden birth count
was 2.018 (MC SE 0.0440), agreeing with its designed value lambda*total_length=2.

For a fully sampled reconstructed tree under the actual birth-death process,
the conditional hidden-birth intensity along a retained branch is 2 lambda E(t).
Locally, either daughter can carry the observed descendants while the other
becomes extinct: lambda dt (D E + E D)/D = 2 lambda E dt.
The corresponding expected total on this two-tip tree is

    2 * integral_0^2 [2 lambda E(t)] dt = 0.1445842.

Thus the current whole-map procedure does NOT sample that posterior. This is
not a new implementation bug: it is a consequence of the previously agreed
ASR-weighted backbone with unrestricted hidden births. Correct extinction-
conditioned side trees do not remove this discrepancy. General multistate
posterior mapping requires descendant likelihood messages as well as E, rather
than replacing lambda alone or relying on marginal ASR weights.

## Numerical stress test

With the default max.step=0.01, pure-death probability errors were:

| mu | Maximum absolute E error |
|---|---:|
| 0.7 | 0.00000606 |
| 100 | 0.0286731 |
| 10,000 | Explicit ODE probability-check failure |

For mu=100 and founder age 0.025, deterministic comparison of the conditional
event-time CDF against its exact formula gave errors 0.0080019, 0.00011184, and
0.00000112 for max.step 0.01, 0.001, and 0.0001, respectively. This shows grid
convergence, but also that ODE tolerances alone do not control interpolation
error. The current settings are not invariant to arbitrary changes of time units.

## Interpretation and next step

The side-tree sampler is supported for the tested ordinary-rate scenarios.
Do not equate that with reliable historical diversity inference for the full
pipeline. Before presenting complete maps as posterior historical reconstructions,
replace the backbone augmentation law with a jointly conditioned law and validate
against the one-state 2 lambda E result. Independently, make interpolation
accuracy rate-aware or adaptive. Keep the present sampler as a clearly labeled
ASR-weighted exploratory option.

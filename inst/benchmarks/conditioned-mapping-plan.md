# Plan: tree-conditioned ClaSSE mapping and numerical safeguards

## Scope and scientific target

Add a distinct tree-conditioned mapping mode; preserve the existing ASR-weighted
mode for reproducibility. First target: fixed parameters, rooted binary observed
tree, state likelihoods at tips, explicit fixed root prior, and **no additional
surviving lineages at present**. This is an extinction-only/full-sampling target.
Keep absolute epoch boundaries and biological extinction E(0)=0.

If the original likelihood used sampling fractions below one, its posterior is a
different target: absence of sampled descendants permits unsampled survivors.
Do not call the extinction-only mode the original fitted-model posterior in that
case. Supporting that sampling process is a separate extension; it is not part
of this first implementation. No extinction-age argument is needed.

Implementation update: `make.simmap.classe.td()` now provides the distinct
`backbone.method="tree-conditioned"` mode. The historical benchmarks below are
bounded prototypes; production validation is described in
`tree-conditioned-validation.md` and `tests/testthat/test-tree-conditioned.R`.

## 1. Shared E and branch-specific D caches

- E_i(t): probability of complete extinction, computed once per parameter set.
- D_e,i(t): likelihood of observed descendants below age t on backbone edge e.
- Compute D in postorder, starting from tip state likelihoods and combining both
  daughter messages at every observed speciation node.
- Integrate the following equations within each epoch, retaining continuity at
  boundaries (r_i is total original event rate, daughter pairs are unordered):

      E_i' = mu_i-r_i E_i+sum_j q_ij E_j+sum_jk lambda_ijk E_j E_k
      D_i' = -r_i D_i+sum_j q_ij D_j
             +sum_jk lambda_ijk(E_j D_k+D_j E_k)

- At an observed binary node, combine ordered daughter messages with
  lambda_ijk (D_left,j D_right,k + D_left,k D_right,j)/2.
  For j=k this reduces to lambda_ijj D_left,j D_right,j.
- Store coherent log scaling for D. Normalizing each grid row independently
  without retaining its scale destroys the time ratios needed for sampling.
- Validate both E and D against the existing compiled likelihood on tiny trees.

## 2. Conditional backbone traversal

- Sample root from root.prior_i * D_root,i, normalized.
- An explicitly supplied, compatible root **posterior** is a distinct input and
  is used directly, not multiplied by D a second time. Do not silently reinterpret
  existing root-ASR inputs as priors. First benchmarks use fixed root priors.
- Sample branches forward in chronological time using no-event probability

      exp(-integral_u^a r_i(v)dv) * D_e,i(u)/D_e,i(a).

  Unlike extinct side trees, there is a nonzero probability of reaching the
  child endpoint without another event; implement this mass explicitly.
- At events use weights q_ij D_j for state transitions and
  lambda_ijk D_j E_k / lambda_ijk D_k E_j for the two hidden-birth orientations.
  Divide by D_i for hazards. There is **no factor 1/2** in these hidden-event
  weights; identical daughter states still have two indistinguishable
  orientations and yield 2 lambda E in the one-state limit.
- At observed nodes sample daughters using the joint lambda/D-left/D-right
  weights. No ASR endpoint rejection or redraw-both-branches loop is needed.
- Simulate each selected hidden side clade with the existing E-conditioned
  direct sampler. Reuse the existing graph, history, and simmap construction.
- Keep internal ASR marginals as validation outputs, not endpoint weights in
  this mode. They do not substitute for descendant likelihood messages.

## 3. Correct empirical input before focus selection

Use the full intended observed tree or complete extracted clade. In particular,
the five-tip Sisyphini tree with inherited full-tree ASRs is not equivalent to
the original 22-tip observed clade: its deleted observed species are not extinct.
Compute D on the complete observed clade/tree, expand side trees only along the
four selected Nesosisyphus paths, and focus-prune after sampling.

Outside-clade information can be carried by a compatible full-tree root posterior
for a complete extracted clade. Internal marginal ASRs cannot restore omitted
descendant data. Use an identical model, sampling target, and root convention
when claiming equivalence between full-tree and extracted-clade runs.

## 4. Numerical improvement

- Seed each epoch's grid using a dimensionless local rate scale, not an absolute
  step of 0.01 regardless of units/rates.
- Adaptively refine intervals using midpoint/quarter-point integration checks.
  Control both absolute E error and log/relative error when E is small and positive.
- Apply equivalent log-likelihood/time-ratio checks to D.
- Check monotonicity of conditional no-event probabilities, not monotonicity of
  E or D themselves. State/epoch changes can make E decrease toward older ages.
- Keep true zeros; distinguish impossible histories from numerical underflow.
- Retry integration with stricter controls on probability-domain violations.
  Roundoff tests near one must reflect relative and absolute solver tolerances.
- Limit grid size and emit convergence diagnostics instead of silently clipping
  or discarding difficult maps.
- Initial analytical-case target: conditional-event CDF error around 1e-6 or
  less, verified under tolerance halving and equivalent time-unit rescaling.

## Benchmarks already executed

### Corrected backbone prototype

Run `inst/benchmarks/benchmark-conditioned-backbone.R` (standalone base R).
Each case has 5,000 conditional prototype draws and 5,000 accepted independent
full forward birth-death trees. The reference retains trees with exactly one
surviving tip, then traces that survivor's ancestry to recover hidden births.

These counts are per single stem, not a whole two-tip crown tree.

| Model | Exact mean hidden births | Conditional prototype | Independent forward rejection |
|---|---:|---:|---:|
| lambda=.5, mu=.05, duration2 | .072292 | .0704 | .0744 |
| lambda=.5, mu=1, duration2 | 1.020240 | 1.0154 | 1.0154 |
| Two epochs, duration3 | 1.261535 | 1.2662 | 1.2480 |

The prototype/reference differences were below 0.82 combined MC standard errors.
Constant-rate count-CDF discrepancies from the predicted Poisson distribution
were at most .00669. Conditional event-age CDF discrepancies were .0089–.0432
for the prototype and .0157–.0273 for the reference (low-extinction runs contain
only ~350 birth events, so age-CDF precision is lower there).

In the epoch test the expected fraction of births younger than the age1 boundary
was .09495, versus .08924 (prototype) and .08974 (reference).

Reference attempts were 14,081, 35,836 and 19,165 respectively. Prototype times
were approximately .03, .19 and .09 seconds; reference times .08, .13 and .26
seconds. There is **no universal speedup claim** from these tiny one-state tests.
They establish the conditional distribution, not full ClaSSE performance.

### Numerical-grid prototypes

Run `inst/benchmarks/benchmark-extinction-grid.R` with current divextra loaded.
The prototype changes only a local function clone, never the production function.

For pure death with mu=100 and founder age .025, the default grid's maximum
conditional-event CDF error was .00803. A rate-scaled step .0001 reduced it to
1.12e-6. For mu=10,000 and age .00025, a step1e-6 gave 1.11e-6 error.

A locally stretched-grid **pure-death proof of concept**, with horizon3 and
tighter ODE tolerance, achieved CDF errors <=1.50e-6 across mu=.7,100,10000 using
540–3,604 grid points. This shows feasibility without a globally tiny uniform
step, but is NOT a validated general multistate adaptive grid.

## Implementation and release gates

1. Implement/test numerically controlled E/D caches and explicit root semantics.
2. Implement exact conditional branch/node sampling as a separate selectable mode.
3. Pass one-state Poisson/count/age tests, including multiple epochs, and the
   mu=0 case with no hidden extinct births.
4. Pass CTMC bridge tests and exact daughter-pair frequency tests.
5. Match marginal root/node probabilities to independent likelihood/ASR results
   under identical assumptions; test joint states on tiny enumeratable examples.
6. Verify mapping lengths, event records, focus selection, complete-clade
   extraction equivalence, and unchanged ASR-weighted-mode seeded outputs.
7. Only then benchmark Nesovinsonia and the complete Sisyphini clade under gr6/gr7:
   initially 10 independent fixed seeds each, expanding only selected paths.
   Report cache time, per-map time, failure reasons, hidden founders, extinct
   tips, state-specific LTT, and numerical-refinement stability. Never replace
   failed seeds or keep only fast/small maps.

Until these gates pass, the current empirical outputs remain ASR-weighted
exploratory simulations rather than tree-conditioned posterior reconstructions.

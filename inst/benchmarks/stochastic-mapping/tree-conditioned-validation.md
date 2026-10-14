# Tree-conditioned ClaSSE mapping: implementation and validation

Implemented in `make.simmap.classe.td()` as the opt-in
`backbone.method="tree-conditioned"`. The existing `"asr-weighted"` default is
unchanged. Core helpers are in `R/conditioned-classe-td.R`; shared extinction
curves and direct side-tree simulation are in `R/extinction-classe-td.R`.

## Target and scope

The target is the supplied reconstructed tree and tip-state information at fixed
parameters, with a user-specified root prior and no additional living descendants
at present. E(0)=0. All epochs retain their absolute ages. This is a complete-census,
extinction-only target, not the original incomplete-sampling posterior when fitted
parameters used sampling fractions below one. Parameter/model uncertainty is not
included by a single call. Extracting a clade with a flat root prior does not retain
ancestral information from the rest of the original tree.

Marginal ASRs are not used as endpoint acceptance weights in this mode. Root
states are weighted by the prior times descendant likelihood. Observed-node
daughter states use both daughter likelihoods. Hidden births use both orientations
lambda D_j E_k and lambda D_k E_j, including two terms when j=k. Side clades are
sampled directly conditional on extinction, not generated and subsequently pruned.

## Production tests

`tests/testthat/test-tree-conditioned.R` includes:

- 1,000 one-state maps on a two-tip age-2 tree, lambda=0.5, mu=0.05.
  Hidden births must follow the Poisson process with rate 2 lambda E(t), not
  lambda. Expected total hidden births per map: 0.1445842; observed: 0.135
  (Monte Carlo SE about 0.0114). The zero-event frequency is also tested.
- 800 two-state maps on a three-tip tree with two epochs. Root likelihoods and
  root-state probabilities agree with the independent compiled ClaSSE likelihood
  and marginal ASR calculation. Sampled internal-node frequencies agree within
  Monte Carlo tolerance. Observed tip states are retained.
- 500 zero-extinction two-state maps: root probabilities match matrix-exponential
  CTMC bridge calculations and no hidden extinct clades are generated.
- Selective side-tree expansion leaves sampled backbone histories unchanged for
  the matched random seed.
- Ambiguous input combinations (ASR probabilities used as priors, incompatible
  side-tree modes) fail explicitly.

The full package suite passed 2,178 checks, with zero failures, warnings, or skips.
A previously saved unrestricted-mode fixture retains identical tree mappings and
identical final RNG state. These checks support the implemented distribution;
they do not constitute a universal proof of numerical accuracy for every parameter
set.

## Numerical safeguards

E and D are cached once per batch, with adaptive grids, epoch splits, and coherent
log scaling for D. Near-zero geometric points and rate-dependent initial spacing
resolve fast dynamics. True zero probabilities are not replaced by positive floors.
Refinement/resource failures stop rather than discard difficult histories.

Default interpolation tolerance is 1e-6, ODE rtol=1e-10, atol=1e-13. D relative
checks additionally budget solver uncertainty (10 rtol + 10 atol/p), because
interpolation refinement cannot resolve ODE noise in tiny probabilities. Absolute
checks remain active. Returned numerical diagnostics report errors and
solver-limited entries. These local checks are not rigorous global error bounds.
For sensitive analyses, tighten interpolation.tol, rtol, and atol and compare
summaries; max.step alone does not control the adaptive grid's accuracy.

## Empirical benchmark

The project script `asr/check-tree-conditioned.R` runs ten maps for each gr6/gr7
fit and complete extracted clade, seed 10 and explicit uniform 16-state root prior.
Nesovinsonia retains both observed tips; Sisyphini retains all 22 observed tips.
Only Nesovinsonia or the four Nesosisyphus paths receive expanded side trees.
Generated extinct-tip counts therefore describe those selected paths, not the
whole clade. Times include cache construction and all ten maps, on this machine;
they are not per-map timings or matched speed comparisons against the old target.

| Clade/model | Maps | Seconds | Extinct tips min / median / max | Log-likelihood difference from compiled ClaSSE |
|---|---:|---:|---:|---:|
| Nesovinsonia gr6 | 10 | 16.017 | 4 / 57.5 / 81 | 4.01e-7 |
| Nesovinsonia gr7 | 10 | 15.146 | 111 / 124.5 / 221 | 4.39e-7 |
| Sisyphini gr6 | 10 | 56.014 | 6 / 13.5 / 19 | -1.59e-6 |
| Sisyphini gr7 | 10 | 57.785 | 14 / 37.5 / 75 | -1.56e-6 |

The script verifies mapped lengths, preservation of all observed tips, no extra
survivors, and no endpoint rejection attempts. Reference likelihoods use complete
sampling, the same root prior, and condition.surv=FALSE. The empirical fixture
files belong to the analysis project and are not distributed with the package.

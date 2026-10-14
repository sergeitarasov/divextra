# Fixed-parameter likelihood validation

Compare log likelihoods at identical parameters and conditioning, without optimization.

- [Random-tree validation: 1,000 cases](runs/random-trees-1000/README.md): comparisons with traditional diversitree GeoSSE and unrestricted ClaSSE on fresh 100--500-tip trees. Includes inputs, source snapshots, results, report and plot. Epoch rates are identical, separated halfway through each tree. The earlier three-comparison, 10--100-tip design is superseded and preserved in `../archive/pilots/likelihood-random-1000-10-100-superseded/`.
- [Historical constant-50 evidence](runs/constant-50-evidence/IMPORT.md): earlier fixed-parameter comparisons at 5%, 25% and 50% boundaries; not new fits.

Keep optimized parameter comparisons under ML benchmarks. Constant-rate references cannot validate genuinely different epoch rates; these require separate mathematical or independently implemented checks.

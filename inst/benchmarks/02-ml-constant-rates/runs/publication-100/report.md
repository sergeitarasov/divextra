# Constant-rate recovery: publication benchmark

Valid fits: 4500/4500; optimizer code 0: 4498; timeouts: 0.
Trees with a valid fit: 1500/1500; strict three-start likelihood agreement: 1496/1500.

Median absolute relative error (%) in the largest size group:

| Scenario | s | d | x |
|---|---:|---:|---:|
| Slow | 3.4 | 6.0 | 14.4 |
| Medium | 3.7 | 4.9 | 14.5 |
| Fast | 3.3 | 6.2 | 15.5 |

The highest finite primary likelihood is used per tree; boundary diagnostics are separate. Strict agreement requires three converged starts and a likelihood spread below 0.0001. Convergence is not proof of accurate recovery or a global optimum.

Both epochs share s, d and x, with a known boundary at half the prescribed stem duration. This tests the episodic implementation under constant rates, not independent epoch-rate recovery. Slow/Medium/Fast are proportional time-scaled processes. Results are conditional on size selection and the living-lineage cap; no crossing or tip-state filter was imposed.

The lower/upper columns are central 95% replicate ranges, NOT confidence intervals or coverage. Focused plots retain all data in calculations and count estimates above the displayed range. See README.md for methods and summary/ for all numerical results.

# Exact-tip GeoSSE simulator pilot

Pilot replicates per scenario: 3

All accepted trees had exact target size: TRUE

The CSV files in this directory contain the complete numerical results.
The realized boundary fraction is measured on the reconstructed crown tree.

## Scenario summary

| Scenario | N simulations | Target tips | Exact | Median attempts | Median extinct | Median boundary fraction | Median runtime (s) |
|---|---:|---:|:---:|---:|---:|---:|---:|
| P01 | 3 | 200 | TRUE | 1 | 41 | 0.487 | 0.012 |
| P02 | 3 | 200 | TRUE | 2 | 30 | 0.610 | 0.010 |
| P03 | 3 | 200 | TRUE | 2 | 35 | 0.415 | 0.010 |
| P04 | 3 | 200 | TRUE | 1 | 100 | 0.204 | 0.014 |
| P05 | 3 | 200 | TRUE | 3 | 90 | 0.698 | 0.016 |
| E25 | 3 | 200 | TRUE | 2 | 30 | 0.179 | 0.012 |
| E75 | 3 | 200 | TRUE | 1 | 30 | 0.892 | 0.008 |
| N050 | 3 | 50 | TRUE | 1 | 2 | 0.535 | 0.003 |
| N500 | 3 | 500 | TRUE | 1 | 64 | 0.672 | 0.023 |

## Interpretation

- Exact target size was recovered in 27 of 27 accepted simulations.
- Accepted trees contained 1 to 383 extinct terminal lineages.
- Accepted simulations required 1 to 4 attempts.
- The realized epoch fraction varies because final tree age is stochastic.
  Boundary calibration should be assessed separately from exact tip count.

# Episodic GeoSSE parameter-recovery benchmark

This benchmark tests whether `divextra::make.classe.td()` recovers known
GeoSSE parameters when rates differ between two time epochs. It is intended
for the supplementary material of the Sunken Islands paper.

The benchmark focuses on the episodic extension. The underlying likelihood is
the ClaSSE likelihood implemented in `diversitree`.

## Questions

The analysis asks whether recovery depends on:

1. which GeoSSE rates differ between epochs;
2. the position of the epoch boundary;
3. the approximate number of extant tips.

The benchmark uses the three GeoSSE states `A`, `B`, and `AB`. Within an epoch,
the model is symmetric between areas:

- `sA = sB = sAB = s`;
- `dA = dB = d`;
- `xA = xB = x`.

The fitted episodic model therefore has six rate parameters:
`s_recent`, `d_recent`, `x_recent`, `s_old`, `d_old`, and `x_old`. The epoch
boundary is fixed to its simulated value.

## Simulation and inference workflow

```mermaid
flowchart TD
    A[Read scenario<br/>scenarios.csv] --> B[Build old and recent<br/>GeoSSE parameter vectors]
    B --> C[Convert GeoSSE to ClaSSE<br/>pars.geosse2classe]
    C --> D[Randomly sample root state<br/>A, B, or AB]
    D --> E[Simulate complete two-epoch process<br/>make.tree.classe.td]

    E --> F{Accept simulation?}
    F -->|Entire process extinct| D
    F -->|Maximum-taxa safety limit reached| D
    F -->|Reconstructed root younger<br/>than epoch boundary| D
    F -->|Accepted| G[Complete simulation table<br/>extant and extinct lineages]

    G --> H[Prune extinct lineages<br/>table2tree]
    H --> I[Observed reconstructed tree<br/>and tip states]

    I --> J[Construct episodic likelihood<br/>make.classe.td]
    I --> K[Construct constant likelihood<br/>make.classe.td]

    L[Read YAML constraints<br/>read_yaml_pars_td] --> M[Create parameter constraints<br/>make_constraints_sse_td]
    M --> J
    M --> K

    J --> N[Create three starting points<br/>init.pars.classe_td]
    K --> N
    N --> O[Maximum-likelihood inference<br/>find.mle]
    O --> P[Root treatment<br/>ROOT.OBS]
    P --> Q[Keep highest likelihood<br/>across EPS 0.1, 0.4, 0.8]
    Q --> R[Compare estimates with<br/>known simulated parameters]
    R --> S[Bias, RMSE, direction,<br/>AIC support and convergence]
```

Every accepted replicate therefore contains a newly simulated tree. The
analysis does not fit repeatedly to one fixed tree.

## Epoch orientation

The simulator and likelihood number epochs in opposite directions.

- In `make.tree.classe.td()`, row 1 is the old, rootward epoch and row 2 is the
  recent, tipward epoch.
- In `make.classe.td()`, parameter suffix `.1` is the recent epoch and suffix
  `.2` is the old epoch.

The scripts preserve this conversion explicitly. Simulation truth from row 1
is compared with fitted parameters ending in `.2`, and simulation truth from
row 2 is compared with parameters ending in `.1`.

## Scenarios

`config/scenarios.csv` contains nine unique conditions.

### Parameter scenarios

Five parameter scenarios use a nominal 50% epoch boundary and approximately
200 extant tips:

1. constant rates;
2. speciation shift;
3. dispersal shift;
4. contraction shift;
5. all three rates differ between epochs.

### Epoch-boundary scenarios

The speciation-shift scenario is also simulated with nominal boundaries at
25% and 75% of the total simulated duration. The 50% condition is already the
baseline parameter scenario.

Because `make.tree.classe.td()` begins with a stem lineage and `table2tree()`
returns the reconstructed crown tree, the realized boundary fraction can
differ slightly from the nominal fraction. Both the absolute boundary age and
the realized fraction of reconstructed tree height are saved.

### Tree-size scenarios

The speciation-shift scenario with a 50% boundary is simulated at approximate
median sizes of 50, 200, and 500 extant tips. Simulation duration, rather than
post hoc selection of an exact tip count, controls tree size. This avoids
conditioning parameter recovery on a narrow observed-tip interval.

## Replication and starting points

The local pilot uses five independently simulated datasets per scenario. The
production analysis uses 100 datasets per scenario, or 900 datasets in total.

Each model is fitted from three starting points, following the strategy used
in the Sunken Islands analyses:

```r
EPS <- c(0.1, 0.4, 0.8)
```

All fits are retained. The fit with the highest log likelihood is used for
parameter-recovery summaries. Agreement among starting points is reported as
an optimization diagnostic.

By default, each dataset receives both:

- the six-rate episodic model;
- a three-rate constant model constrained across epochs.

The constant comparison is represented as two epochs with all rates shared.
This gives the same likelihood as ordinary time-homogeneous ClaSSE while
keeping the same epoch boundary and likelihood machinery. The current
`make.classe.td(..., n.epoch = 1)` path is not used because it fails inside the
time-dependent branch evaluator; this is a separate compatibility issue to
repair in DivExtra.

## Root and survival conditioning

At the start of every simulation attempt, the root state is sampled uniformly
from `A`, `B`, and `AB`. If that process goes entirely extinct, a new attempt
draws a new root state and a new tree.

Inference does not receive the simulated root state. It uses the
Maddison-FitzJohn treatment implemented in `diversitree` as `ROOT.OBS`. Under
this treatment, the root weights are proportional to the conditional
likelihood of the observed data in each root state. Fits use
`condition.surv = TRUE`, matching the retention of surviving simulated trees.

Trees whose reconstructed root is younger than the epoch boundary contain no
information from the old epoch and are resimulated. Entirely extinct
simulations are also resimulated. No filtering on the final number of tips is
used in the main benchmark.

## Directory contents

```text
geosse-parameter-recovery/
├── README.md
├── config/
│   ├── model-constant.yml
│   ├── model-episodic.yml
│   └── scenarios.csv
├── R/
│   ├── benchmark-functions.R
│   ├── 00-calibrate-tree-sizes.R
│   ├── 01-run-one.R
│   ├── 02-run-local.R
│   └── 03-summarize.R
├── roihu/
│   ├── run-array.sh
│   └── submit.sh
└── results/
    ├── calibrated/
    ├── raw/
    ├── summaries/
    └── logs/
```

The scripts create the `results` subdirectories when needed. Full raw results
should not be committed to the package or submitted to CRAN. Keep them in the
Sunken Islands analysis repository or archive them separately. Only scripts,
configuration files, and compact final summary tables belong in the package.

## Software setup

Run the scripts from the root of the DivExtra repository. Install the current
checkout before starting:

```sh
R CMD INSTALL .
```

The required R packages are `divextra`, `diversitree`, and `ape`.

## Step 1: calibrate tree sizes

Calibration estimates a fixed simulation duration for each scenario so that
the median number of extant tips is close to its target. It uses a small
stochastic bisection search and does not create datasets used in parameter
recovery.

```sh
Rscript inst/benchmarks/geosse-parameter-recovery/R/00-calibrate-tree-sizes.R
```

The resulting file is:

```text
inst/benchmarks/geosse-parameter-recovery/results/calibrated/scenarios-calibrated.csv
```

Inspect the calibrated durations and median tip counts before continuing.

## Step 2: smoke test one replicate per scenario

```sh
Rscript inst/benchmarks/geosse-parameter-recovery/R/02-run-local.R 1
```

This creates nine simulated datasets and fits each from three starting points.

## Step 3: local pilot with five replicates

```sh
Rscript inst/benchmarks/geosse-parameter-recovery/R/02-run-local.R 5
```

Existing successful result files are skipped, so the command can be resumed.
To rerun a dataset, remove its individual RDS file first.

## Step 4: summarize the pilot

```sh
Rscript inst/benchmarks/geosse-parameter-recovery/R/03-summarize.R
```

The script writes dataset-level estimates, start-level optimization results,
and scenario-level recovery summaries to `results/summaries/`.

Before the Roihu run, inspect:

- the distribution of actual tip counts;
- simulation and optimization failure rates;
- agreement among the three starting points;
- estimates at numerical boundaries;
- runtime as a function of tree size;
- whether the old and recent epochs both contain enough tree history.

## Step 5: send the production analysis to Roihu

Copy or clone the DivExtra repository on Roihu, install the same package
checkout, and run calibration once if the calibrated CSV was not copied.

The production design contains 900 array tasks:

```text
9 scenarios x 100 replicates = 900 datasets
```

Review the Slurm account, partition, memory, and time settings in
`roihu/run-array.sh`, then submit:

```sh
bash inst/benchmarks/geosse-parameter-recovery/roihu/submit.sh
```

One array task simulates and fits one dataset. Each task writes one RDS file,
so interrupted arrays can be submitted again safely. The array script limits
simultaneous jobs with `%25`; adjust this only after checking Roihu policy and
the pilot runtime.

## Outputs

Each raw RDS file stores:

- scenario and replicate identifiers;
- random seed;
- simulated stem-root state and reconstructed crown-root state;
- true GeoSSE rates;
- simulated tree and simulation metadata;
- actual number of extant tips;
- tree height and epoch boundary;
- all fits from all starting points;
- best episodic and constant fits;
- parameter estimates, likelihood, AIC, convergence code, and runtime;
- R and package versions.

The summary script reports parameter bias, relative bias, RMSE, empirical 95%
quantiles, recovery of the direction of epoch change, AIC support for the
episodic model, and optimization agreement across starting points.

## Interpretation

The empirical 95% interval across recovered estimates describes the sampling
distribution of the estimator. It is not a confidence interval calculated for
an individual dataset.

The constant scenario measures false detection of an epoch effect. The
nonconstant scenarios measure recovery of effect direction and magnitude.
Poor recovery in small trees or near extreme epoch boundaries should be
reported as a limitation rather than removed from the analysis.


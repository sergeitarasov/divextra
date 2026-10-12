#!/bin/bash
#SBATCH --job-name=geosse-recovery
#SBATCH --array=1-900%25
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=inst/benchmarks/geosse-parameter-recovery/results/logs/%A_%a.out
#SBATCH --error=inst/benchmarks/geosse-parameter-recovery/results/logs/%A_%a.err

set -euo pipefail

# Add the appropriate Roihu account and partition above if required.
# Submit from the root of the DivExtra repository.

REPOSITORY_DIR="${REPOSITORY_DIR:-$PWD}"
cd "$REPOSITORY_DIR"

SCENARIO_NUMBER=$(( (SLURM_ARRAY_TASK_ID - 1) / 100 + 1 ))
REPLICATE_NUMBER=$(( (SLURM_ARRAY_TASK_ID - 1) % 100 + 1 ))

SCENARIO_FILE="inst/benchmarks/geosse-parameter-recovery/results/calibrated/scenarios-calibrated.csv"
SCENARIO_ROW=$((SCENARIO_NUMBER + 1))
SCENARIO_ID=$(awk -F, -v row="$SCENARIO_ROW" 'NR == row {print $1}' "$SCENARIO_FILE")

if [ -z "$SCENARIO_ID" ]; then
  echo "Could not resolve scenario for task $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

Rscript \
  inst/benchmarks/geosse-parameter-recovery/R/01-run-one.R \
  "$SCENARIO_ID" \
  "$REPLICATE_NUMBER"


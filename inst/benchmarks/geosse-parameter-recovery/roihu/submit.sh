#!/bin/bash

set -euo pipefail

# Run this script from the root of the DivExtra repository.

BENCHMARK_DIR="inst/benchmarks/geosse-parameter-recovery"
CALIBRATION_FILE="$BENCHMARK_DIR/results/calibrated/scenarios-calibrated.csv"

if [ ! -f "$CALIBRATION_FILE" ]; then
  echo "Missing calibration file: $CALIBRATION_FILE" >&2
  echo "Run 00-calibrate-tree-sizes.R or copy the reviewed local calibration." >&2
  exit 1
fi

mkdir -p "$BENCHMARK_DIR/results/logs"
mkdir -p "$BENCHMARK_DIR/results/raw"
mkdir -p "$BENCHMARK_DIR/results/summaries"

sbatch "$BENCHMARK_DIR/roihu/run-array.sh"


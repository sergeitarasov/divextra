#!/bin/bash
#SBATCH --job-name=geosse-pub-constant
#SBATCH --account=project_2004996
#SBATCH --partition=small
#SBATCH --array=1-150%6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%A_%a.out
set -eo pipefail
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load r-env/461
set -u
export R_LIBS_USER=/projappl/project_2004996/geosse-benchmark/Rlib/4.6.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
srun Rscript --vanilla run-task.R

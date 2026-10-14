#!/bin/bash
#SBATCH --job-name=geosse-pub-summary
#SBATCH --account=project_2004996
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=summary-%j.out
set -eo pipefail
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load r-env/461
srun Rscript --vanilla analyse.R

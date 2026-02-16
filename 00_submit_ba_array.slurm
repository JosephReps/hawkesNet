#!/bin/bash
#SBATCH --job-name=hawkesNet_ba350
#SBATCH --account=uoa04249
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-350
#SBATCH --output=logs/ba_%A_%a.out
#SBATCH --error=logs/ba_%A_%a.err

set -euo pipefail

module purge
module load R

export R_LIBS_USER="/home/jrep150/R/foss-2023a/4.3"
export R_LIBS="/home/jrep150/R/foss-2023a/4.3"

OUTDIR="results/ba_run_${SLURM_ARRAY_JOB_ID}"
mkdir -p logs "${OUTDIR}"

# If your package is not installed on the compute nodes, you can run from a checked-out repo:
# export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3
# Or use renv, etc.

Rscript simstudy_ba/run_ba_rep.R "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}"

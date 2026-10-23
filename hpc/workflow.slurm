#!/usr/bin/env bash
#SBATCH --job-name=synthetic_qc
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/qc-%j.out
#SBATCH --error=logs/qc-%j.err
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
printf 'Host: '; hostname
printf 'Slurm job: %s\n' "${SLURM_JOB_ID:?}"
# Activate site-approved modules/environment BEFORE submission; sbatch inherits PATH.
# Optional first argument selects summary, reads, genotype, or all (default).
bash run.sh --stage "${1:-all}" --threads "${SLURM_CPUS_PER_TASK:-2}" --out "runs/slurm-${SLURM_JOB_ID}"

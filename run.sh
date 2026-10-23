#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec python3 scripts/run_workflow.py "$@"

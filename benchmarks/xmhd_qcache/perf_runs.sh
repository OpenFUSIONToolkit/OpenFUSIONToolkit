#!/bin/bash
# Perf sampling runs used for phase0_baseline.md (serial process, all cores via OpenMP unless noted).
# Usage: bash perf_runs.sh <build_dir> [perf_exe]
BUILD=${1:?build dir}
PERF=${2:-$(command -v perf)}
cd "$(dirname "$0")"
run() { python run_bench.py --build-dir "$BUILD" --perf --perf-exe "$PERF" "$@" | tail -1; }
# Jacobian builds only (prof_reps=5,0,0) and NL residual only (1,100,0)
run --tag perf_build --size 0.5 --prof-reps 5,0,0 --cases 'tet_p2,tet_p3_mg,tet_p4_mg,hex_p3,curved_p3_mg'
run --tag perf_nlop --size 0.5 --prof-reps 1,100,0 --cases 'tet_p2,tet_p3,tet_p4,hex_p3,curved_p3_mg'
# Single thread (no lock contention) for comparison
run --tag perf_build_omp1 --size 0.5 --prof-reps 3,0,0 --cases 'tet_p3_mg_omp1'
run --tag perf_nlop_omp1 --size 0.5 --prof-reps 1,50,0 --cases 'tet_p3_mg_omp1'
# Full time steps
run --tag perf_run --cases 'tet_p3_mg_run,tet_p3_mg_mfnk_run'

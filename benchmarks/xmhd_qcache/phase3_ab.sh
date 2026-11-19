#!/bin/bash
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Phase 3 A/B runs: xmhd_qcache=0 (direct) vs 2 (REFERENCE) vs 1 (FULL), interleaved repetitions.
#   bash phase3_ab.sh <build_dir> [nrep]
# Results: results/p3ab_t1_q<mode>_r<rep>.json (1 thread, NL residual only, profile cases) and
#          results/p3ab_t4_q<mode>_r<rep>.json (4 threads, default matrix incl. run/MFNK cases).
# Summarize with: python phase3_tables.py
BUILD=${1:?build dir}
NREP=${2:-2}
HERE=$(cd $(dirname $0) && pwd)
cd $HERE
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
T1CASES='tet_p2$,tet_p3$,tet_p4$,hex_p3$,curved_p3_mg$'
for rep in $(seq 1 $NREP); do
  for q in 0 2; do
    python3 run_bench.py --build-dir $BUILD --tag p3ab_t1_q${q}_r${rep} --cases "$T1CASES" --threads 1 \
      --prof-reps 1,10,0 --xmhd-opt xmhd_qcache=$q
  done
done
python3 run_bench.py --build-dir $BUILD --tag p3ab_t1_q1_r1 --cases "$T1CASES" --threads 1 \
  --prof-reps 1,10,0 --xmhd-opt xmhd_qcache=1
for rep in $(seq 1 $NREP); do
  for q in 0 2; do
    python3 run_bench.py --build-dir $BUILD --tag p3ab_t4_q${q}_r${rep} --threads 4 --xmhd-opt xmhd_qcache=$q
  done
done
python3 run_bench.py --build-dir $BUILD --tag p3ab_t4_q1_r1 --threads 4 --xmhd-opt xmhd_qcache=1

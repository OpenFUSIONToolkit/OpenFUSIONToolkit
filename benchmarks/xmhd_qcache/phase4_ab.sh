#!/bin/bash
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Phase 4 A/B runs (Jacobian on the quadrature cache + CSR position map), interleaved repetitions.
#   bash phase4_ab.sh <build_dir> <new_xmhd_bench> <phase3_xmhd_bench> [nrep] [part]
# <new_xmhd_bench> is a copy of the branch's xmhd_bench (copy it so rebuilding does not disturb the runs), and
# <phase3_xmhd_bench> an xmhd_bench built from the Phase 3 branch (original Jacobian). <build_dir> is only used
# for the build metadata.
# Configurations ("p3" = Phase 3 binary with xmhd_qcache=2, ie. residual on the cache, original Jacobian):
#   p3      : original Jacobian
#   q0pm0   : new direct path (shared kernel + batched contractions), CSR search assembly
#   q0pm1   : new direct path with the CSR position map
#   q2pm0   : cache path, CSR search assembly
#   q2pm1   : cache path with the CSR position map (new default)
#   q2pm1h  : as q2pm1 with the hand-written GEMM loop instead of BLAS dgemm
# part=t1  : 1 thread, profile cases, builds only (prof_reps 2,1,0)
# part=t4  : 4 threads, default matrix (profile, MPI, run/MFNK cases), p3 vs q2pm1
# Results: results/p4ab_<part>_<cfg>_r<rep>.json. Tables: python phase4_tables.py
BUILD=${1:?build dir}
NEWEXE=$(realpath ${2:?new xmhd_bench})
P3EXE=$(realpath ${3:?phase 3 xmhd_bench})
NREP=${4:-2}
PART=${5:-all}
HERE=$(cd $(dirname $0) && pwd)
cd $HERE
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
T1CASES='tet_p2$,tet_p3$,tet_p4$,hex_p3$,curved_p3_mg$,tet_p3_mg$'
run_cfg() { # tag threads cases prof_reps cfg
  local tag=$1 thr=$2 cases=$3 reps=$4 cfg=$5
  local repsopt=""
  [ -n "$reps" ] && repsopt="--prof-reps $reps"
  case $cfg in
    p3)     python3 run_bench.py --build-dir $BUILD --exe $P3EXE --tag $tag --cases "$cases" --threads $thr $repsopt --xmhd-opt xmhd_qcache=2 ;;
    q0pm0)  python3 run_bench.py --build-dir $BUILD --exe $NEWEXE --tag $tag --cases "$cases" --threads $thr $repsopt --xmhd-opt xmhd_qcache=0 --xmhd-opt xmhd_jac_posmap=0 ;;
    q0pm1)  python3 run_bench.py --build-dir $BUILD --exe $NEWEXE --tag $tag --cases "$cases" --threads $thr $repsopt --xmhd-opt xmhd_qcache=0 --xmhd-opt xmhd_jac_posmap=1 ;;
    q2pm0)  python3 run_bench.py --build-dir $BUILD --exe $NEWEXE --tag $tag --cases "$cases" --threads $thr $repsopt --xmhd-opt xmhd_qcache=2 --xmhd-opt xmhd_jac_posmap=0 ;;
    q2pm1)  python3 run_bench.py --build-dir $BUILD --exe $NEWEXE --tag $tag --cases "$cases" --threads $thr $repsopt --xmhd-opt xmhd_qcache=2 --xmhd-opt xmhd_jac_posmap=1 ;;
    q2pm1h) python3 run_bench.py --build-dir $BUILD --exe $NEWEXE --tag $tag --cases "$cases" --threads $thr $repsopt --xmhd-opt xmhd_qcache=2 --xmhd-opt xmhd_jac_posmap=1 --bench-opt jac_blas=F ;;
  esac
}
if [ "$PART" = "all" ] || [ "$PART" = "t1" ]; then
  for rep in $(seq 1 $NREP); do
    for cfg in p3 q0pm0 q0pm1 q2pm0 q2pm1 q2pm1h; do
      run_cfg p4ab_t1_${cfg}_r${rep} 1 "$T1CASES" 2,1,0 $cfg
    done
  done
fi
if [ "$PART" = "all" ] || [ "$PART" = "t4" ]; then
  for rep in $(seq 1 $NREP); do
    for cfg in p3 q2pm1; do
      run_cfg p4ab_t4_${cfg}_r${rep} 4 '.*' "" $cfg
    done
  done
fi

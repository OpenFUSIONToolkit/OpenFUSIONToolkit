#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Bucket a flat `perf report --no-children --stdio` self-time profile into cost categories.

Usage: python perf_summary.py <perf_report.txt> [more reports ...]

Categories are matched in order against the symbol name (first match wins). OFT is built
with LTO, so some small routines are inlined into their callers; their cost then shows
up under the caller (eg. under xmhd_build_ops or xmhd_errmatrix_apply).
'''
import re
import sys

CATEGORIES = [
    ('log2phys', r'log2phys'),
    ('mesh jacobian', r'mesh_jacobian|_MOD_\w*jacobian|cell_is_curved|hcurl_get_cgops'),
    ('basis: lagrange', r'lag_basis|oft_lag_eval|oft_lag_geval|lag_eval_all|lag_geval_all'),
    ('basis: h1 grad', r'h1_basis|oft_h1_geval|oft_h1_eval|h1_geval_all'),
    ('basis: hcurl', r'hcurl_basis|oft_hcurl_eval|oft_hcurl_ceval|hcurl_eval_all|hcurl_ceval_all'),
    ('omp lock/critical/barrier wait', r'gomp_mutex|GOMP_critical|omp_set_lock|gomp_\w*lock|futex|gomp_barrier|'
                                       r'gomp_team_barrier|do_wait|do_spin|gomp_sem|GOMP_barrier'),
    ('local->global scatter (mat_add + CSR column search)', r'mat_add_local|add_values|fem_mat_add|oft_sort_MOD_search'),
    ('xmhd_interp_apply (field reconstruction)', r'xmhd_interp_apply|cc_interp'),
    ('xmhd_build_ops (self)', r'xmhd_build_ops'),
    ('xmhd_errmatrix_apply (self)', r'xmhd_errmatrix_apply'),
    ('other xmhd', r'xmhd'),
    ('linear solver / factorization', r'ilu|superlu|lu_|dgstrf|dgstrs|sp_|gmres|native_matrix_apply|mat_apply|'
                                      r'csr|dgemv|dgemm|daxpy|ddot|blas|lapack|openblas|vec_|_MOD_\w*vector'),
    ('math lib (pow/exp/sqrt)', r'^(pow|exp|log|sqrt|__ieee|__pow|__exp|__log)'),
    ('malloc/memcpy/kernel', r'malloc|free|memcpy|memset|memmove|\[k\]|page|clear_page|_int_'),
]
LINE = re.compile(r'^\s*([\d.]+)%\s+(\S+)\s+(\S+)\s+\[([.k])\]\s+(.*)$')


def summarize(path):
    totals = {}
    other = {}
    with open(path) as fid:
        for line in fid:
            m = LINE.match(line)
            if not m:
                continue
            pct, dso, sym = float(m.group(1)), m.group(3), m.group(5).strip()
            if m.group(4) == 'k':
                sym = '[k] ' + sym
            if 'libgomp' in dso:
                # libgomp has no symbols here: spin/lock/barrier waits (critical sections, locks, implicit barriers)
                sym = 'GOMP_critical (libgomp)'
            for cat, pat in CATEGORIES:
                if re.search(pat, sym):
                    totals[cat] = totals.get(cat, 0.) + pct
                    break
            else:
                other[sym] = other.get(sym, 0.) + pct
    print('== {0}'.format(path))
    for cat, _ in CATEGORIES:
        if cat in totals:
            print('  {0:45s} {1:6.1f}%'.format(cat, totals[cat]))
    print('  {0:45s} {1:6.1f}%'.format('other (listed >= 0.3%)', sum(other.values())))
    for sym, pct in sorted(other.items(), key=lambda x: -x[1])[:8]:
        print('      {0:6.2f}%  {1}'.format(pct, sym[:90]))
    return totals


if __name__ == '__main__':
    for p in sys.argv[1:]:
        summarize(p)

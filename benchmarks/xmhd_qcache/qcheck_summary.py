#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Summarize the xMHD residual equivalence checks (Phase 3) from a verbose pytest log

Usage: OMP_NUM_THREADS=2 ./run_test.sh physics/MUG/test_MUG_qcache.py -v -s > qcheck.log
       python qcheck_summary.py qcheck.log

For every test case and mode pair (0:2 direct vs REFERENCE, 0:1 direct vs FULL, 2:1 REFERENCE vs FULL)
prints the largest difference over the three residual evaluations (dt>0, dt>0 with back_step, dt<0)
relative to the largest entry of the same equation (B, V, n, T, Te, n2, j2), the largest difference of the
Grad(H^1) block relative to its own maximum (a near-cancellation), and the largest relative difference of
the diagnostic integrals.
'''
import re
import sys

EQ = {'Bc': 'B', 'Bg': 'B', 'Vx': 'V', 'Vy': 'V', 'Vz': 'V'}


def summarize(lines):
    tests = []
    name = None
    for line in lines:
        m = re.match(r'^\S*test_MUG_qcache\.py::(\S+)', line)
        if m:
            name = m.group(1)
            tests.append({'name': name, 'blocks': [], 'diags': []})
            continue
        if (name is None) or (not line.startswith('XMHD_QCHECK')):
            continue
        vals = line.split()
        if vals[1] == 'field':
            case = re.search(r'case=(\S+)', line).group(1)
            modes = re.search(r'modes=(\S+)', line).group(1)
            block = re.search(r'block=(\S+)', line).group(1)
            tests[-1]['blocks'].append((case, modes, block) + tuple(float(v) for v in vals[-3:]))
        elif vals[1] == 'diag':
            tests[-1]['diags'].append((vals[2].split('=')[1], float(vals[-1])))
    rows = []
    for t in tests:
        if not t['blocks']:
            continue
        row = {'name': t['name']}
        for pair in ('0:2', '0:1', '2:1'):
            rel_eq, rel_bg = 0., 0.
            for case in set(b[0] for b in t['blocks']):
                group = [b for b in t['blocks'] if b[0] == case and b[1] == pair]
                scale = {}
                for b in group:
                    eq = EQ.get(b[2], b[2])
                    scale[eq] = max(scale.get(eq, 0.), b[3])
                for b in group:
                    rel_eq = max(rel_eq, b[4]/max(scale[EQ.get(b[2], b[2])], 1.E-300))
                    if b[2] == 'Bg':
                        rel_bg = max(rel_bg, b[5])
            diag = max([d[1] for d in t['diags'] if d[0] == pair] + [0.])
            row[pair] = (rel_eq, rel_bg, diag)
        rows.append(row)
    return rows


if __name__ == '__main__':
    with open(sys.argv[1]) as fid:
        rows = summarize(fid.readlines())
    print('| case | 0 vs 2: blocks | 0 vs 2: Bg (own max) | 0 vs 2: diag | 0 vs 1: blocks | 2 vs 1: blocks |')
    print('|---|---|---|---|---|---|')
    worst = [0., 0., 0., 0., 0.]
    for r in rows:
        vals = [r['0:2'][0], r['0:2'][1], r['0:2'][2], r['0:1'][0], r['2:1'][0]]
        worst = [max(a, b) for a, b in zip(worst, vals)]
        print('| `{0}` | {1:.1e} | {2:.1e} | {3:.1e} | {4:.1e} | {5:.1e} |'.format(r['name'], *vals))
    print('| **max** | {0:.1e} | {1:.1e} | {2:.1e} | {3:.1e} | {4:.1e} |'.format(*worst))

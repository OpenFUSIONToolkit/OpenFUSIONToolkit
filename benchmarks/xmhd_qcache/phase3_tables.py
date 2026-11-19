#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Median A/B tables for Phase 3 (xmhd_qcache=0 direct, 2 REFERENCE, 1 FULL) from results/p3ab_*.json

Run after phase3_ab.sh. Prints markdown tables: NL residual time per apply (profile counter and loop timer),
cache setup time and memory, Jacobian build time and RSS (unchanged by the cache; noise check), and for the
run-mode cases the time per step, iteration counts and final-field checksums.
'''
import glob
import json
import os
import re
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))


def load(prefix):
    '''Return {mode: {case: [records over repetitions]}} for results/<prefix>_q<mode>_r<rep>.json'''
    out = {}
    for path in sorted(glob.glob(os.path.join(HERE, 'results', prefix + '_q*_r*.json'))):
        m = re.search(r'_q(\d)_r(\d+)\.json$', path)
        mode = int(m.group(1))
        with open(path) as fid:
            for rec in json.load(fid)['results']:
                if rec.get('returncode') != 0:
                    continue
                out.setdefault(mode, {}).setdefault(rec['name'], []).append(rec)
    return out


def med(recs, key, per_call=False):
    vals = []
    for r in recs:
        v = r.get(key)
        if v is None:
            continue
        if per_call:
            n = r.get('prof_nlop_n')
            if not n:
                continue
            v = v/n
        vals.append(v)
    return statistics.median(vals) if vals else None


def fmt(v, f='{0:.3f}'):
    return '-' if v is None else f.format(v)


def ratio(a, b):
    return '-' if (a is None or b is None or b == 0) else '{0:.2f}'.format(a/b)


def rel(a, b):
    if a is None or b is None:
        return '-'
    return '{0:.1e}'.format(abs(a-b)/max(abs(a), 1.e-300))


def profile_table(data, title):
    print('\n### ' + title + '\n')
    print('| case | NL apply q0 [s] | q2 [s] | q0/q2 | q1 [s] | q0/q1 | loop/apply q0 [s] | q2 [s] | q0/q2 | '
          'setup q2 [s] | mem q2 [MB] | setup q1 [s] | mem q1 [MB] | build q0 [s] | build q2 [s] | RSS q0/q2/q1 [MB] | nlop_norm rel q2, q1 |')
    print('|' + '---|'*17)
    cases = [c for c in data.get(0, {}) if not c.endswith('_run')]
    for c in cases:
        r0, r2, r1 = data.get(0, {}).get(c, []), data.get(2, {}).get(c, []), data.get(1, {}).get(c, [])
        n0, n2, n1 = med(r0, 'nlop_s'), med(r2, 'nlop_s'), med(r1, 'nlop_s')
        l0, l2 = med(r0, 'prof_nlop_loop_total_s', True), med(r2, 'prof_nlop_loop_total_s', True)
        print('| `{0}` | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} | {9} | {10} | {11} | {12} | {13} | {14} | {15}/{16}/{17} | {18}, {19} |'.format(
            c, fmt(n0, '{0:.4f}'), fmt(n2, '{0:.4f}'), ratio(n0, n2), fmt(n1, '{0:.4f}'), ratio(n0, n1),
            fmt(l0, '{0:.4f}'), fmt(l2, '{0:.4f}'), ratio(l0, l2),
            fmt(med(r2, 'qcache_setup_s'), '{0:.4f}'), fmt(med(r2, 'qcache_mem_mb'), '{0:.1f}'),
            fmt(med(r1, 'qcache_setup_s'), '{0:.3f}'), fmt(med(r1, 'qcache_mem_mb'), '{0:.0f}'),
            fmt(med(r0, 'build_s'), '{0:.2f}'), fmt(med(r2, 'build_s'), '{0:.2f}'),
            fmt(med(r0, 'max_rss_mb'), '{0:.0f}'), fmt(med(r2, 'max_rss_mb'), '{0:.0f}'), fmt(med(r1, 'max_rss_mb'), '{0:.0f}'),
            rel(med(r0, 'nlop_norm'), med(r2, 'nlop_norm')), rel(med(r0, 'nlop_norm'), med(r1, 'nlop_norm'))))
    print('\nRepetitions (q0/q2/q1): ' + ', '.join('{0}: {1}/{2}/{3}'.format(
        c, len(data.get(0, {}).get(c, [])), len(data.get(2, {}).get(c, [])), len(data.get(1, {}).get(c, []))) for c in cases))


def run_table(data, title):
    cases = [c for c in data.get(0, {}) if c.endswith('_run')]
    if not cases:
        return
    print('\n### ' + title + '\n')
    print('| case | mode | step time (excl. 1st) [s] | run time [s] | NL residual total [s] | NL applies | lits | nlits | '
          'per-step its = q0 | final B, V, Ti rel. diff vs q0 |')
    print('|' + '---|'*10)
    for c in cases:
        base = data[0][c][0]
        its0 = [(s['lits'], s['nlits']) for s in base.get('steps', [])]
        for mode in (0, 2, 1):
            recs = data.get(mode, {}).get(c, [])
            if not recs:
                continue
            its = [[(s['lits'], s['nlits']) for s in r.get('steps', [])] for r in recs]
            same = all(i == its0 for i in its)
            diffs = ', '.join(rel(base.get(k), med(recs, k)) for k in ('final_norm_B', 'final_norm_V', 'final_norm_Ti'))
            print('| `{0}` | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} | {9} |'.format(
                c, mode, fmt(med(recs, 'step_time_avg_excl1_s')), fmt(med(recs, 'run_time_s')),
                fmt(med(recs, 'prof_nlop_total_s')), fmt(med(recs, 'prof_nlop_n'), '{0:.0f}'),
                fmt(med(recs, 'lits_total'), '{0:.0f}'), fmt(med(recs, 'nlits_total'), '{0:.0f}'), same, diffs))


if __name__ == '__main__':
    profile_table(load('p3ab_t1'), '1 thread (NL residual only: prof_reps=1,10,0)')
    t4 = load('p3ab_t4')
    profile_table(t4, '4 threads (default matrix, profile cases)')
    run_table(t4, '4 threads, run mode (5 time steps; NK and MFNK)')

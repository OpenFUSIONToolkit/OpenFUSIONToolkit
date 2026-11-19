#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Median A/B tables for Phase 4 (Jacobian on the quadrature cache, CSR position map) from results/p4ab_*.json

Run after phase4_ab.sh. Configurations: p3 (Phase 3 binary: original Jacobian), q0pm0/q0pm1 (new direct path
without/with the CSR position map), q2pm0/q2pm1 (cache path without/with the position map; q2pm1 is the
default), q2pm1h (default with the hand-written GEMM loop). Prints markdown tables of the Jacobian build time
per level, the setup/memory of the caches and position maps, RSS, and for the run-mode cases the step times,
per-step iteration counts and final-field checksums.
'''
import glob
import json
import os
import re
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
CFGS = ['p3', 'q0pm0', 'q0pm1', 'q2pm0', 'q2pm1', 'q2pm1h']


def load(prefix):
    '''Return {cfg: {case: [records over repetitions]}} for results/<prefix>_<cfg>_r<rep>.json'''
    out = {}
    for path in sorted(glob.glob(os.path.join(HERE, 'results', prefix + '_*_r*.json'))):
        m = re.search(prefix + r'_(\w+?)_r(\d+)\.json$', path)
        if m is None:
            continue
        cfg = m.group(1)
        with open(path) as fid:
            for rec in json.load(fid)['results']:
                if rec.get('returncode') != 0:
                    continue
                out.setdefault(cfg, {}).setdefault(rec['name'], []).append(rec)
    return out


def med(vals):
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def medk(recs, key):
    return med([r.get(key) for r in recs])


def level_avg(recs, which):
    '''Median build time per call of the finest level ('fine') or the sum over coarser levels ('coarse')'''
    vals = []
    for r in recs:
        lv = r.get('levels', [])
        if not lv:
            continue
        top = max(x['level'] for x in lv)
        if which == 'fine':
            vals.append(sum(x['avg_s'] for x in lv if x['level'] == top))
        else:
            c = [x['avg_s'] for x in lv if x['level'] != top and x['finemesh']]
            vals.append(sum(c) if c else None)
    return med(vals)


def fmt(v, f='{0:.3f}'):
    return '-' if v is None else f.format(v)


def ratio(a, b):
    return '-' if (a is None or b is None or b == 0) else '{0:.2f}'.format(a/b)


def rel(a, b):
    if a is None or b is None:
        return '-'
    return '{0:.1e}'.format(abs(a-b)/max(abs(a), 1.e-300))


def build_table(data, title, cfgs):
    print('\n### ' + title + '\n')
    cases = [c for c in data.get('q2pm1', {}) if not c.endswith('_run')]
    hdr = '| case | level | ' + ' | '.join('{0} [s]'.format(c) for c in cfgs) + ' | p3/q2pm1 | q2pm0/q2pm1 (posmap) | p3/q2pm0 (cache+GEMM) | q0pm1/q2pm1 |'
    print(hdr)
    print('|' + '---|'*(len(cfgs)+6))
    for c in cases:
        for which in ('fine', 'coarse'):
            vals = {cfg: level_avg(data.get(cfg, {}).get(c, []), which) for cfg in cfgs}
            if vals.get('q2pm1') is None:
                continue
            print('| `{0}` | {1} | {2} | {3} | {4} | {5} | {6} |'.format(
                c, which, ' | '.join(fmt(vals[cfg]) for cfg in cfgs),
                ratio(vals.get('p3'), vals.get('q2pm1')), ratio(vals.get('q2pm0'), vals.get('q2pm1')),
                ratio(vals.get('p3'), vals.get('q2pm0')), ratio(vals.get('q0pm1'), vals.get('q2pm1'))))
    print('\nRepetitions: ' + ', '.join('{0}: {1}'.format(cfg, '/'.join(str(len(data.get(cfg, {}).get(c, []))) for c in cases))
                                      for cfg in cfgs))


def mem_table(data, title):
    print('\n### ' + title + '\n')
    print('| case | cells | qcache [MB] | coarse-level caches [MB] | posmap (all levels) [MB] | posmap KB/cell (finest) | '
          'posmap setup [s] | RSS p3 [MB] | RSS q2pm0 [MB] | RSS q2pm1 [MB] | matvec_norm rel q2pm1 vs p3 |')
    print('|' + '---|'*11)
    for c, recs in data.get('q2pm1', {}).items():
        if c.endswith('_run'):
            continue
        r = recs[0]
        lv = r.get('levels', [])
        top = max(x['level'] for x in lv) if lv else None
        fine_pm = sum(x.get('posmap_mb', 0.) for x in lv if x['level'] == top)
        ncell = r.get('ncells_local')
        kb = 1.e3*fine_pm/ncell if ncell else None
        p3 = data.get('p3', {}).get(c, [])
        print('| `{0}` | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} | {9} | {10} |'.format(
            c, ncell, fmt(r.get('qcache_mem_mb'), '{0:.1f}'), fmt(r.get('jcache_coarse_mem_mb'), '{0:.1f}'),
            fmt(r.get('posmap_mem_mb'), '{0:.0f}'), fmt(kb, '{0:.1f}'), fmt(r.get('posmap_setup_s'), '{0:.2f}'),
            fmt(medk(p3, 'max_rss_mb'), '{0:.0f}'), fmt(medk(data.get('q2pm0', {}).get(c, []), 'max_rss_mb'), '{0:.0f}'),
            fmt(medk(recs, 'max_rss_mb'), '{0:.0f}'), rel(medk(p3, 'matvec_norm'), medk(recs, 'matvec_norm'))))


def run_table(data, title):
    cases = [c for c in data.get('q2pm1', {}) if c.endswith('_run')]
    if not cases:
        return
    print('\n### ' + title + '\n')
    print('| case | cfg | step time (excl. 1st) [s] | run time [s] | Jacobian builds total [s] | builds | precond. update total [s] | '
          'NL residual total [s] | lits | nlits | per-step its = p3 | final B, V, Ti rel. diff vs p3 |')
    print('|' + '---|'*12)
    for c in cases:
        base = data.get('p3', {}).get(c, [None])[0]
        its0 = [(s['lits'], s['nlits']) for s in base.get('steps', [])] if base else None
        for cfg in ('p3', 'q0pm1', 'q2pm1'):
            recs = data.get(cfg, {}).get(c, [])
            if not recs:
                continue
            its = [[(s['lits'], s['nlits']) for s in r.get('steps', [])] for r in recs]
            same = all(i == its0 for i in its) if its0 is not None else '-'
            diffs = ', '.join(rel(base.get(k) if base else None, medk(recs, k))
                              for k in ('final_norm_B', 'final_norm_V', 'final_norm_Ti'))
            print('| `{0}` | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} | {9} | {10} | {11} |'.format(
                c, cfg, fmt(medk(recs, 'step_time_avg_excl1_s')), fmt(medk(recs, 'run_time_s')),
                fmt(medk(recs, 'prof_setops_total_s')), fmt(medk(recs, 'prof_setops_n'), '{0:.0f}'),
                fmt(medk(recs, 'prof_preupdate_total_s')), fmt(medk(recs, 'prof_nlop_total_s')),
                fmt(medk(recs, 'lits_total'), '{0:.0f}'), fmt(medk(recs, 'nlits_total'), '{0:.0f}'), same, diffs))


if __name__ == '__main__':
    t1 = load('p4ab_t1')
    build_table(t1, '1 thread: Jacobian build time per call (prof_reps=2,1,0)', CFGS)
    mem_table(t1, '1 thread: memory')
    t4 = load('p4ab_t4')
    build_table(t4, '4 threads: Jacobian build time per call (default matrix)', ['p3', 'q2pm1'])
    mem_table(t4, '4 threads: memory')
    run_table(t4, '4 threads, run mode (5 time steps; NK and MFNK)')

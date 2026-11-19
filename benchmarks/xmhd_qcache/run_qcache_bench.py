#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Phase 2 kernel benchmark driver for the quadrature-point cache (src/fem/fem_quad_cache.F90).

Runs `fem_qcache_bench` (src/tests/fem/fem_qcache_bench.F90, built with OFT_BUILD_TESTS) on the
meshes/orders of the Phase 0 harness (`run_bench.py`) and writes results/<tag>.json plus Markdown tables.

Examples:
  python run_qcache_bench.py --build-dir ../../builds/build_release --tag phase2_kernels
  python run_qcache_bench.py --build-dir ... --cases tet_p3,hex_p3 --threads 1 --tag x
  python run_qcache_bench.py --report results/phase2_kernels.json
'''
import argparse
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_bench import cube_mesh, cyl_mesh, CUBE_NI  # noqa: E402

OFT_IN = """&runtime_options
 ppn=1
 debug=0
/

&mesh_options
{mesh_options}/

{extra_mesh_groups}
&qcache_bench_options
 order={order}
 nreps={nreps}
/
"""

VARIANTS = ['direct', 'full', 'ref', 'refc', 'fullc', 'refcel']
KERNELS = ['basis', 'recon', 'resid']
VARIANT_DESC = {
    'direct': 'direct evaluators (current xMHD)',
    'full': 'FULL, per-qp kernels',
    'ref': 'REFERENCE, per-qp kernels, gathered rows',
    'refc': 'REFERENCE, per-qp kernels, canonical frame where available',
    'fullc': 'FULL, cell-batched kernels',
    'refcel': 'REFERENCE, cell-batched kernels',
}


def cases(size=1.0):
    def ni(n):
        return max(1, int(round(n*size)))
    return [
        {'name': 'tet_p2', 'mesh': cube_mesh(ni(CUBE_NI[2])), 'order': 2},
        {'name': 'tet_p3', 'mesh': cube_mesh(ni(CUBE_NI[3])), 'order': 3},
        {'name': 'tet_p4', 'mesh': cube_mesh(ni(CUBE_NI[4])), 'order': 4},
        {'name': 'hex_p3', 'mesh': cube_mesh(ni(CUBE_NI['hex3']), hex_mesh=True), 'order': 3},
        {'name': 'curved_p3', 'mesh': cyl_mesh(), 'order': 3},
    ]


def parse(out):
    rec = {'times': {}, 'sums': {}}
    m = re.search(r'QCACHE_BENCH info \S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', out)
    if m:
        keys = ['ncells', 'nqp', 'nthreads', 'nce_lag', 'nce_h1', 'nce_hcurl', 'nvar']
        rec.update({k: int(v) for k, v in zip(keys, m.groups())})
    m = re.search(r'QCACHE_BENCH geom setup_s,bytes\s+(\S+)\s+(\d+)', out)
    if m:
        rec['geom_setup_s'] = float(m.group(1))
        rec['geom_bytes'] = int(m.group(2))
    for m in re.finditer(r'QCACHE_BENCH cache (\w+)\s+' + r'\s+'.join([r'(\S+)']*6), out):
        vals = [float(v) for v in m.groups()[1:]]
        rec['cache_' + m.group(1)] = {'setup_s': vals[0:3], 'bytes': vals[3:6]}
    for m in re.finditer(r'QCACHE_BENCH xform space,nident,nperm,ndense,nnone,resid\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)', out):
        rec.setdefault('xform', {})[int(m.group(1))] = [int(v) for v in m.groups()[1:5]] + [float(m.group(6))]
    for m in re.finditer(r'QCACHE_BENCH time (\w+) (\w+)[ \t]+(\S+)((?:[ \t]+\S+)+)', out):
        rec['times'].setdefault(m.group(1), {})[m.group(2)] = [float(v) for v in m.group(4).split()]
        rec['sums'].setdefault(m.group(1), {})[m.group(2)] = float(m.group(3))
    return rec


def run_case(case, args, exe):
    tag = '{0}_t{1}'.format(case['name'], args.threads)
    rdir = os.path.join(HERE, 'runs', 'qcache_' + args.tag, tag)
    os.makedirs(rdir, exist_ok=True)
    src_root = os.path.abspath(os.path.join(HERE, '..', '..', 'src'))
    for f in case['mesh']['mesh_files']:
        shutil.copy(os.path.join(src_root, f), rdir)
    with open(os.path.join(rdir, 'oft.in'), 'w') as fid:
        fid.write(OFT_IN.format(mesh_options=case['mesh']['mesh_options'],
                                extra_mesh_groups=case['mesh']['extra_mesh_groups'],
                                order=case['order'], nreps=args.nreps))
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    cmd = [exe, 'oft.in']
    if shutil.which('/usr/bin/time'):
        cmd = ['/usr/bin/time', '-v'] + cmd
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=rdir, env=env, capture_output=True, text=True, timeout=args.timeout)
    with open(os.path.join(rdir, 'output.log'), 'w') as fid:
        fid.write(proc.stdout)
        fid.write(proc.stderr)
    rec = parse(proc.stdout)
    rec.update({'name': case['name'], 'order': case['order'], 'threads': args.threads,
                'returncode': proc.returncode, 'wall_s': time.time()-t0})
    m = re.search(r'Maximum resident set size \(kbytes\): (\d+)', proc.stderr)
    if m:
        rec['max_rss_mb'] = int(m.group(1))/1024.
    return rec


def fmt_ms(vals):
    return '{0:.2f} ({1:.2f})'.format(1.E3*statistics.median(vals), 1.E3*min(vals))


def report(results):
    lines = []
    for rec in results:
        if rec.get('returncode', 1) != 0 or not rec['times']:
            lines.append('{0}: FAILED (returncode {1})'.format(rec['name'], rec.get('returncode')))
            continue
        lines.append('### {0} ({1} cells, {2} qp/cell, {3} thread(s); nce lag/h1/hcurl = {4}/{5}/{6})'.format(
            rec['name'], rec['ncells'], rec['nqp'], rec['threads'], rec['nce_lag'], rec['nce_h1'], rec['nce_hcurl']))
        lines.append('')
        lines.append('Time per pass over all cells and quadrature points, median (min) of {0} interleaved repetitions, ms. '
                     'Ratio = median / direct median.'.format(len(rec['times']['direct']['basis'])))
        lines.append('')
        lines.append('| variant | basis | ratio | recon | ratio | resid | ratio |')
        lines.append('|---|---|---|---|---|---|---|')
        for v in VARIANTS:
            if v not in rec['times']:
                continue
            row = [v]
            for k in KERNELS:
                vals = rec['times'][v][k]
                ref = statistics.median(rec['times']['direct'][k])
                row += [fmt_ms(vals), '{0:.2f}'.format(statistics.median(vals)/ref)]
            lines.append('| ' + ' | '.join(row) + ' |')
        # Checksums
        dev = 0.
        for k in ('recon', 'resid'):
            base = rec['sums']['direct'][k]
            for v in VARIANTS:
                if v in rec['sums']:
                    dev = max(dev, abs(rec['sums'][v][k]-base)/max(abs(base), 1.E-300))
        lines.append('')
        lines.append('Max relative checksum difference vs direct (recon, resid): {0:.1e}'.format(dev))
        lines.append('')
    return '\n'.join(lines)


def memory_table(results):
    lines = ['| case | cells | qp | geom setup [s] | geom [MB] | FULL setup lag/h1/hcurl [s] | FULL [MB] | FULL [KB/cell] '
             '| REF setup lag/h1/hcurl [s] | REF tables [MB] | max RSS [MB] |',
             '|---|---|---|---|---|---|---|---|---|---|---|']
    for rec in results:
        if rec.get('returncode', 1) != 0 or 'cache_full' not in rec:
            continue
        full = rec['cache_full']
        refc = rec['cache_ref']
        fb = sum(full['bytes'])
        lines.append('| {0} | {1} | {2} | {3:.4f} | {4:.2f} | {5} | {6:.1f} | {7:.1f} | {8} | {9:.2f} | {10} |'.format(
            rec['name'], rec['ncells'], rec['nqp'], rec['geom_setup_s'], rec['geom_bytes']/1.E6,
            '/'.join('{0:.3f}'.format(v) for v in full['setup_s']), fb/1.E6, fb/rec['ncells']/1.E3,
            '/'.join('{0:.4f}'.format(v) for v in refc['setup_s']), sum(refc['bytes'])/1.E6,
            '{0:.0f}'.format(rec['max_rss_mb']) if 'max_rss_mb' in rec else '-'))
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--build-dir', help='OFT build directory (with tests/fem/fem_qcache_bench)')
    parser.add_argument('--tag', default='phase2_kernels')
    parser.add_argument('--cases', default='', help='Comma-separated case names (default: all)')
    parser.add_argument('--threads', type=int, default=os.cpu_count())
    parser.add_argument('--nreps', type=int, default=7)
    parser.add_argument('--size', type=float, default=1.0)
    parser.add_argument('--timeout', type=int, default=3600)
    parser.add_argument('--report', help='Print tables for an existing result file and exit')
    args = parser.parse_args()
    if args.report:
        with open(args.report) as fid:
            data = json.load(fid)
        print(report(data['results']))
        print(memory_table(data['results']))
        return
    exe = os.path.join(os.path.abspath(args.build_dir), 'tests', 'fem', 'fem_qcache_bench')
    sel = [c for c in cases(args.size) if (not args.cases) or (c['name'] in args.cases.split(','))]
    # Results are merged into results/<tag>.json after every case (same name+threads replaced)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    out_file = os.path.join(HERE, 'results', args.tag + '.json')
    data = {'meta': {}, 'results': []}
    if os.path.exists(out_file):
        with open(out_file) as fid:
            data = json.load(fid)
    data['meta'] = {'host': platform.node(), 'cpu': platform.processor(), 'ncores': os.cpu_count(),
                    'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'nreps': args.nreps}
    for case in sel:
        print('Running {0} ({1} threads)...'.format(case['name'], args.threads), flush=True)
        rec = run_case(case, args, exe)
        data['results'] = [r for r in data['results']
                           if not (r['name'] == rec['name'] and r['threads'] == rec['threads'])] + [rec]
        with open(out_file, 'w') as fid:
            json.dump(data, fid, indent=1)
        print(report([rec]), flush=True)
    print(memory_table(data['results']))


if __name__ == '__main__':
    main()

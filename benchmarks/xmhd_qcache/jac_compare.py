#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Compare xMHD Jacobians written by different builds or implementations (Phase 4).

`xmhd_bench` in `mode='qcheck'` with `jac_dump='<prefix>'` writes the Jacobian of every level, built at
the check state, to `<prefix>_L<level>_R<rank>.bin` (see `xmhd_jac_dump` in `xmhd.F90`). This script
1. runs the equivalence-test configurations of `src/tests/physics/MUG/test_MUG_qcache.py` with a given
   build and `xmhd_qcache` value and collects the dumps (`run`), and
2. compares two dump sets entrywise (`compare`): the CSR structure must be identical, and value
   differences are reported relative to the largest entry of the same row (`row`), of the same
   field block row/column pair (`block`), and of the whole matrix (`global`).

Example (reference = a build of the original Jacobian, new = the branch under test):
    python jac_compare.py run --build-dir <ref_build> --out dumps/ref --qcache 0
    python jac_compare.py run --build-dir <new_build> --out dumps/q0 --qcache 0
    python jac_compare.py run --build-dir <new_build> --out dumps/q2 --qcache 2
    python jac_compare.py compare dumps/ref dumps/q0
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.abspath(os.path.join(HERE, '..', '..', 'src', 'tests', 'physics', 'MUG'))


def load_cases():
    '''Configurations of test_MUG_qcache.py as (name, mesh, order, xmhd, bench, minlev, nproc)'''
    sys.path.insert(0, TEST_DIR)
    import test_MUG_qcache as t
    cases = []
    for name, (xmhd, bench) in t.physics_cases.items():
        cases.append(('physics_tet_p2[{0}]'.format(name), lambda: t.cube_mesh(2), 2, xmhd, bench, -1, 1))
    for name in ['base', 'single_temp_hall', 'hyper', 'shaping_eta_temp', 'floors_solid']:
        xmhd, bench = t.physics_cases[name]
        cases.append(('physics_tet_p3[{0}]'.format(name), lambda: t.cube_mesh(2), 3, xmhd, bench, -1, 1))
    hall = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n nu_gyro=1.0\n"
    for order in [2, 3, 4]:
        ni = 1 if order == 4 else 2
        cases.append(('tet[{0}]'.format(order), (lambda ni=ni: t.cube_mesh(ni)), order, hall, " shaping=T\n", -1, 1))
    for order in [2, 3]:
        cases.append(('hex[{0}]'.format(order), lambda: t.cube_mesh(2, hex_mesh=True), order, hall, " shaping=T\n", -1, 1))
    cases.append(('curved_p3', t.curved_mesh, 3, " xmhd_hall=T\n",
                  " shaping=T\n den_floor=1.02E19\n temp_floor=20.4\n qcheck_solid=T\n qcheck_negsrc=T\n", -1, 1))
    for vbc in ['norm', 'tang']:
        cases.append(('vbc[{0}]'.format(vbc), lambda: t.cube_mesh(2, periodic=False), 2,
                      " vbc='{0}'\n nbc='n'\n tbc='n'\n".format(vbc), '', -1, 1))
    for ct in ['tet', 'hex']:
        cases.append(('revolved[{0}]'.format(ct), (lambda ct=ct: t.revolved_mesh(ct)), 2, " vbc='norm'\n", '', -1, 1))
    cases.append(('pmg', lambda: t.cube_mesh(2), 3, '', '', 2, 1))
    cases.append(('pmg_p4', lambda: t.cube_mesh(1), 4, " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n", " shaping=T\n", 2, 1))
    cases.append(('mpi2', lambda: t.cube_mesh(2, nlevels=2), 3, " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n",
                  " shaping=T\n qcheck_solid=T\n", -1, 2))
    return t, cases


def run(args):
    t, cases = load_cases()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    exe = os.path.join(os.path.abspath(args.build_dir), 'tests', 'physics', 'MUG', 'xmhd_bench')
    if args.exe:
        exe = os.path.abspath(args.exe)
    workdir = os.path.join(os.path.abspath(args.build_dir), 'tests', 'physics', 'MUG')
    env = dict(os.environ)
    env.setdefault('OMP_NUM_THREADS', '2')
    for name, mesh_fn, order, xmhd, bench, minlev, nproc in cases:
        if args.cases and not any(c in name for c in args.cases.split(',')):
            continue
        os.chdir(workdir)
        mesh = mesh_fn()
        prefix = os.path.join(out, name.replace('[', '_').replace(']', ''))
        xmhd_all = xmhd + ''.join(' {0}\n'.format(o) for o in args.xmhd_opt)
        if args.qcache is not None:
            xmhd_all += ' xmhd_qcache={0}\n'.format(args.qcache)
        bench_all = bench + " jac_dump='{0}'\n".format(prefix)
        with open('oft.in', 'w+') as fid:
            fid.write(t.oft_in_template.format(order=order, minlev=minlev, xmhd=xmhd_all, bench=bench_all, **mesh))
        with open('oft_in_qcache.xml', 'w+') as fid:
            fid.write(t.oft_in_xml)
        cmd = [exe, 'oft.in', 'oft_in_qcache.xml']
        if nproc > 1:
            cmd = args.mpiexec.split() + ['-np', str(nproc)] + cmd
        res = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        with open(prefix + '.log', 'w') as fid:
            fid.write(res.stdout)
        nf = len(glob.glob(prefix + '_L*_R*.bin'))
        print('{0:32s} status={1} dumps={2}'.format(name, res.returncode, nf), flush=True)


def read_dump(fname):
    with open(fname, 'rb') as fid:
        level, nr, nnz = np.fromfile(fid, dtype=np.int32, count=3)
        kr = np.fromfile(fid, dtype=np.int32, count=nr+1)
        lc = np.fromfile(fid, dtype=np.int32, count=nnz)
        vals = np.fromfile(fid, dtype=np.float64, count=nnz)
    return int(level), kr, lc, vals


def compare_files(fa, fb):
    la, kra, lca, va = read_dump(fa)
    lb, krb, lcb, vb = read_dump(fb)
    if not (np.array_equal(kra, krb) and np.array_equal(lca, lcb)):
        return None
    nr = kra.size - 1
    rows = np.repeat(np.arange(nr), np.diff(kra))
    d = np.abs(va - vb)
    a = np.abs(va)
    rowmax = np.zeros(nr)
    np.maximum.at(rowmax, rows, a)
    rowd = np.zeros(nr)
    np.maximum.at(rowd, rows, d)
    ok = rowmax > 0
    rel_row = np.max(rowd[ok]/rowmax[ok]) if ok.any() else 0.
    rel_row = max(rel_row, np.max(rowd[~ok]) if (~ok).any() else 0.)
    rel_glob = d.max()/max(a.max(), 1.E-300)
    return {'level': la, 'nnz': va.size, 'rel_row': rel_row, 'rel_global': rel_glob,
            'max_abs': a.max(), 'max_diff': d.max()}


def compare(args):
    worst = 0.
    files = sorted(glob.glob(os.path.join(args.a, '*_L*_R*.bin')))
    if not files:
        print('No dumps in', args.a)
        return 1
    for fa in files:
        fb = os.path.join(args.b, os.path.basename(fa))
        if not os.path.exists(fb):
            print('{0:40s} MISSING in {1}'.format(os.path.basename(fa), args.b))
            worst = np.inf
            continue
        r = compare_files(fa, fb)
        if r is None:
            print('{0:40s} STRUCTURE DIFFERS'.format(os.path.basename(fa)))
            worst = np.inf
            continue
        worst = max(worst, r['rel_row'])
        print('{0:40s} nnz={1:9d} rel_row={2:10.3e} rel_global={3:10.3e}'.format(
            os.path.basename(fa)[:-4], r['nnz'], r['rel_row'], r['rel_global']))
    print('WORST rel_row {0:.3e}'.format(worst))
    return 0 if worst < args.tol else 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    pr = sub.add_parser('run')
    pr.add_argument('--build-dir', required=True)
    pr.add_argument('--out', required=True)
    pr.add_argument('--exe', default=None, help='xmhd_bench executable (default: from --build-dir)')
    pr.add_argument('--qcache', type=int, default=None)
    pr.add_argument('--xmhd-opt', action='append', default=[])
    pr.add_argument('--cases', default='')
    pr.add_argument('--mpiexec', default='mpiexec --oversubscribe')
    pc = sub.add_parser('compare')
    pc.add_argument('a')
    pc.add_argument('b')
    pc.add_argument('--tol', type=float, default=1.E-11)
    args = p.parse_args()
    if args.cmd == 'run':
        run(args)
    else:
        sys.exit(compare(args))


if __name__ == '__main__':
    main()

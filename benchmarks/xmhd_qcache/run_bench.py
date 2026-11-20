#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Benchmark harness for the 3D xMHD (MUG) module: Jacobian build, NL residual and time-step timings.

Runs the `xmhd_bench` driver (src/tests/physics/MUG/xmhd_bench.F90, built with OFT_BUILD_TESTS)
over a case matrix and writes machine-readable results (JSON + CSV). See README.md.

Examples:
  # Default case matrix with the default options (quadrature cache, CSR position map)
  python run_bench.py --build-dir <build_dir> --tag cache
  # Same matrix with direct basis evaluation and column-search assembly
  python run_bench.py --build-dir <build_dir> --tag direct --xmhd-opt xmhd_qcache=0 --xmhd-opt xmhd_jac_posmap=0
  # Compare two result files
  python run_bench.py --compare results/direct.json results/cache.json
'''
import argparse
import csv
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

#------------------------------------------------------------------------------
# Input templates
#------------------------------------------------------------------------------
OFT_IN = """&runtime_options
 ppn=1
 debug={debug}
 use_petsc={petsc}
/

&mesh_options
{mesh_options}/

{extra_mesh_groups}
&xmhd_bench_options
 order={order}
 minlev={minlev}
 mode='{mode}'
 prof_reps={prof_reps}
 lam={lam}
{bench_extra}/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=T
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_ohmic=T
 xmhd_visc_heat=T
 xmhd_upwind=T
 bbc='bc'
 vbc='all'
 nbc='d'
 tbc='d'
 dt={dt}
 eta=1.0
 nu_visc=10.0
 nu_par=100.0
 d_dens=1.0
 kappa_par=1000.0,1000.0
 kappa_perp=1.0,1.0
 nsteps={nsteps}
 rst_freq=100000
 lin_tol=1.E-8
 nl_tol=1.E-5
 nu_xmhd={nu_xmhd}
 nclean=100000
 rst_ind=0
 maxextrap=2
 ittarget=40
 xmhd_mfnk={mfnk}
 xmhd_prefreq=100
{xmhd_extra}/
"""

# Single-level preconditioner: GMRES(its) with block-Jacobi/{solver} (LU as the MUG tests, ILU as the Spheromak examples)
XML_SINGLE = """<oft>
  <xmhd>
    <pre type="gmres">
      <its>{its}</its>
      <nrits>{its}</nrits>
      <pre type="block_jacobi">
        <nlocal>1</nlocal>
        <solver type="{solver}"></solver>
      </pre>
    </pre>
  </xmhd>
</oft>
"""

# p-multigrid preconditioner (as src/tests/physics/MUG/oft_in_xmhd_mg.xml, with a selectable local solver)
XML_MG = """<oft>
  <xmhd>
    <pre type="mg">
      <smoother direction="up">
        <solver type="gmres">
          <its>2</its>
          <nrits>2</nrits>
          <pre type="block_jacobi">
            <nlocal>1</nlocal>
            <solver type="{solver}"></solver>
          </pre>
        </solver>
      </smoother>
    </pre>
  </xmhd>
</oft>
"""


#------------------------------------------------------------------------------
# Case matrix
#------------------------------------------------------------------------------
def cube_mesh(ni, hex_mesh=False, nlevels=1, nbase=1):
    '''Triply periodic box [-0.5,0.5]^2 x [-1,1] from the built-in cube generator (cad_type=92)'''
    mesh = " meshname='cube'\n cad_type=92\n nlevels={0}\n nbase={1}\n grid_order=1\n".format(nlevels, nbase)
    cube = "&cube_options\n mesh_type={0}\n ni={1},{1},{2}\n rscale=1.,1.,2.\n shift=-0.5,-0.5,-1.\n ref_per=T,T,T\n/\n".format(
        2 if hex_mesh else 1, ni, 2*ni)
    return {'mesh_options': mesh, 'extra_mesh_groups': cube, 'lam': '1.,1.,2.', 'mesh_files': [], 'mgdim': nlevels}


def cyl_mesh(nlevels=1, nbase=1):
    '''Curved (quadratic) tetrahedral cylinder mesh from the Spheromak_tilt example (R=1, L=2)'''
    mesh = (" meshname='cyl'\n cad_type=0\n nlevels={0}\n nbase={1}\n grid_order=2\n fix_boundary=T\n".format(nlevels, nbase))
    native = "&native_mesh_options\n filename='cyl_tilt.h5'\n/\n"
    return {'mesh_options': mesh, 'extra_mesh_groups': native, 'lam': '2.,2.,2.',
            'mesh_files': ['examples/MUG/Spheromak_tilt/cyl_tilt.h5'], 'mgdim': nlevels}


def make_case(name, mesh, order, mg=False, mfnk=False, mode='profile', nsteps=0, threads=None, nprocs=1,
              prof_reps='3,5,10', dt='2.E-7'):
    return {'name': name, 'mesh': mesh, 'order': order, 'mg': mg, 'mfnk': mfnk, 'mode': mode, 'nsteps': nsteps,
            'threads': threads, 'nprocs': nprocs, 'prof_reps': prof_reps, 'dt': dt}


def case_matrix(ncores, size=1.0):
    """Default case matrix (sized for ~1-3 minutes per case on a 4-core Xeon); `size` scales cube resolution"""
    def ni(n):
        # Tet cube meshes need a power of 2
        return max(1, 2**int(round(math.log2(max(1., n*size)))))
    c = []
    # --- Profile mode, tets: p=2 (single level), p=3/4 with and without p-multigrid (coarsest p=2)
    c.append(make_case('tet_p2', cube_mesh(ni(CUBE_NI[2])), 2))
    for order in (3, 4):
        mesh = cube_mesh(ni(CUBE_NI[order]))
        c.append(make_case('tet_p{0}'.format(order), mesh, order))
        c.append(make_case('tet_p{0}_mg'.format(order), mesh, order, mg=True))
    # --- Hex (p-MG down to p=2 is disabled for hex in the MUG tests) and curved tets
    c.append(make_case('hex_p3', cube_mesh(ni(CUBE_NI['hex3']), hex_mesh=True), 3))
    c.append(make_case('curved_p3_mg', cyl_mesh(), 3, mg=True))
    # --- Thread scaling (profile mode, p=3 MG)
    c.append(make_case('tet_p3_mg_omp1', cube_mesh(ni(CUBE_NI[3])), 3, mg=True, threads=1))
    # --- MPI (2 ranks x ncores/2 threads; level 1 is the serial base mesh, level 2 its distributed copy)
    c.append(make_case('tet_p3_mg_mpiref', cube_mesh(ni(CUBE_NI['mpi'])), 3, mg=True))  # Serial reference for mpi2
    c.append(make_case('tet_p3_mg_mpi2', cube_mesh(ni(CUBE_NI['mpi']), nlevels=2, nbase=1), 3, mg=True, nprocs=2,
                       threads=max(1, ncores//2)))
    # --- Time stepping (smaller mesh): single level p=2, and p=3 MG with assembled Newton-Krylov vs matrix-free
    c.append(make_case('tet_p2_run', cube_mesh(ni(CUBE_NI['run'])), 2, mode='run', nsteps=NSTEPS))
    c.append(make_case('tet_p3_mg_run', cube_mesh(ni(CUBE_NI['run'])), 3, mg=True, mode='run', nsteps=NSTEPS))
    c.append(make_case('tet_p3_mg_mfnk_run', cube_mesh(ni(CUBE_NI['run'])), 3, mg=True, mfnk=True, mode='run',
                       nsteps=NSTEPS))
    return c


def quick_matrix(ncores, size=1.0):
    """Tiny smoke-test matrix"""
    return [make_case('quick_tet_p2', cube_mesh(2), 2, prof_reps='1,2,2'),
            make_case('quick_tet_p3_mg_run', cube_mesh(2), 3, mg=True, mode='run', nsteps=2)]


# Cube resolution (boxes in x,y; 2x in z) per case type, and time steps for run cases
CUBE_NI = {2: 8, 3: 8, 4: 4, 'hex3': 4, 'mpi': 4, 'run': 4}
NSTEPS = 5


#------------------------------------------------------------------------------
# Running and parsing
#------------------------------------------------------------------------------
def fe_minlev(case):
    """FE minimum level: -1 (finest only) or the p=2 level on the finest mesh (mgdim+1) for p-multigrid.

    A p=1 coarse level (minlev=mgdim) gives NaNs in the MG preconditioner for xMHD, and the MUG regression tests
    never use it, so p-multigrid here coarsens to p=2 (as in test_alfven.py/test_sound.py)."""
    if not case['mg']:
        return -1
    return case['mesh']['mgdim'] + case.get('mg_coarse_p', 2) - 1


# Smoother iterations for each p-level above the coarsest (coarsest first), following the MUG tests
NU_PLEVELS = {2: ['2'], 3: ['2', '2'], 4: ['10', '4', '2'], 5: ['10', '4', '2', '2']}


def nu_xmhd(case):
    """Smoother iterations per FE level (index 1..mgdim+order-1); zero below minlev"""
    if not case['mg']:
        return '-1'
    minlev = fe_minlev(case)
    nlev = case['mesh']['mgdim'] + case['order'] - 1
    vals = ['0']*(minlev-1) + NU_PLEVELS[case['order']][-(nlev-minlev+1):]
    return ','.join(vals)


def write_inputs(case, run_dir, src_dir, args):
    mesh = case['mesh']
    xmhd_extra = ''.join(' {0}\n'.format(o) for o in args.xmhd_opt)
    bench_extra = ''.join(' {0}\n'.format(o) for o in args.bench_opt)
    with open(os.path.join(run_dir, 'oft.in'), 'w') as fid:
        fid.write(OFT_IN.format(debug=args.debug, petsc='T' if args.petsc else 'F', mesh_options=mesh['mesh_options'],
                                extra_mesh_groups=mesh['extra_mesh_groups'], order=case['order'],
                                minlev=fe_minlev(case), mode=case['mode'], prof_reps=case['prof_reps'],
                                lam=mesh['lam'], bench_extra=bench_extra, dt=case['dt'], nsteps=max(1, case['nsteps']),
                                nu_xmhd=nu_xmhd(case), mfnk='T' if case['mfnk'] else 'F', xmhd_extra=xmhd_extra))
    with open(os.path.join(run_dir, 'oft_in.xml'), 'w') as fid:
        if case['mg']:
            fid.write(XML_MG.format(solver=args.pre))
        else:
            fid.write(XML_SINGLE.format(solver=args.pre, its=(1 if args.pre == 'lu' else 8)))
    for f in mesh['mesh_files']:
        shutil.copy(os.path.join(src_dir, f), run_dir)


FLOAT = r'([-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?)'


def parse_output(text):
    r = {}
    for key, pat in (('dot_s', r'^\s*Dot\s*=\s*' + FLOAT), ('matvec_s', r'^\s*MatVec\s*=\s*' + FLOAT),
                     ('nlop_s', r'^\s*NL Op\s*=\s*' + FLOAT), ('build_s', r'^\s*Build\s*=\s*' + FLOAT),
                     ('b_dofs', r'# B-DOFs\s*=\s*' + FLOAT), ('v_dofs', r'# V-DOFs\s*=\s*' + FLOAT),
                     ('n_dofs', r'# N-DOFs\s*=\s*' + FLOAT)):
        m = re.search(pat, text, re.M)
        if m:
            r[key] = float(m.group(1))
    for m in re.finditer(r'^XMHD_BENCH (\w+)=\s*(\S+)', text, re.M):
        try:
            r[m.group(1)] = float(m.group(2)) if '.' in m.group(2) or 'E' in m.group(2) else int(m.group(2))
        except ValueError:
            r[m.group(1)] = m.group(2)
    for m in re.finditer(r'^XMHD_PROF (setops|nlop|preupdate) ncalls,total_s\s+(\d+)\s+' + FLOAT, text, re.M):
        r['prof_{0}_n'.format(m.group(1))] = int(m.group(2))
        r['prof_{0}_total_s'.format(m.group(1))] = float(m.group(3))
    m = re.search(r'^XMHD_PROF qcache mode,setup_s,mem_mb\s+(\d+)\s+' + FLOAT + r'\s+' + FLOAT, text, re.M)
    if m:
        r['qcache_mode'], r['qcache_setup_s'], r['qcache_mem_mb'] = int(m.group(1)), float(m.group(2)), float(m.group(3))
    for m in re.finditer(r'^XMHD_PROF (nlop_norm|matvec_norm)\s+' + FLOAT, text, re.M):
        r[m.group(1)] = float(m.group(2))
    levels = []
    for m in re.finditer(r'^XMHD_PROF build level,order,finemesh,ncalls,total_s,avg_s,nrows,nnz\s+(\d+)\s+(\d+)\s+([TF])\s+(\d+)\s+'
                         + FLOAT + r'\s+' + FLOAT + r'\s+(-?\d+)\s+(-?\d+)', text, re.M):
        levels.append({'level': int(m.group(1)), 'order': int(m.group(2)), 'finemesh': m.group(3) == 'T',
                       'ncalls': int(m.group(4)), 'total_s': float(m.group(5)), 'avg_s': float(m.group(6)),
                       'nrows': int(m.group(7)), 'nnz': int(m.group(8))})
    for m in re.finditer(r'^XMHD_PROF jac level,qc_builds,posmap_mb,posmap_setup_s\s+(\d+)\s+(\d+)\s+' + FLOAT + r'\s+' + FLOAT,
                         text, re.M):
        for lv in levels:
            if lv['level'] == int(m.group(1)):
                lv['qc_builds'] = int(m.group(2))
                lv['posmap_mb'] = float(m.group(3))
                lv['posmap_setup_s'] = float(m.group(4))
    m = re.search(r'^XMHD_PROF jcache coarse_mem_mb,posmap\s+' + FLOAT + r'\s+(\d+)', text, re.M)
    if m:
        r['jcache_coarse_mem_mb'] = float(m.group(1))
        r['jac_posmap'] = int(m.group(2))
    if levels:
        r['posmap_mem_mb'] = sum(lv.get('posmap_mb', 0.) for lv in levels)
        r['posmap_setup_s'] = sum(lv.get('posmap_setup_s', 0.) for lv in levels)
        r['levels'] = levels
        top = max(lv['level'] for lv in levels)
        r['build_fine_avg_s'] = sum(lv['avg_s'] for lv in levels if lv['level'] == top)
        r['build_coarse_avg_s'] = sum(lv['avg_s'] for lv in levels if lv['level'] != top)
        r['jac_nnz_total'] = sum(max(lv['nnz'], 0) for lv in levels)
    steps = []
    for m in re.finditer(r'^\s*Timestep\s+(\d+)\s+' + FLOAT + r'\s+(\d+)\s+(\d+)\s+' + FLOAT + r'\s+' + FLOAT, text, re.M):
        steps.append({'step': int(m.group(1)), 't': float(m.group(2)), 'lits': int(m.group(3)),
                      'nlits': int(m.group(4)), 'time_s': float(m.group(5)), 'dt': float(m.group(6))})
    if steps:
        r['steps'] = steps
        r['step_time_avg_s'] = sum(s['time_s'] for s in steps)/len(steps)
        r['step_time_avg_excl1_s'] = (sum(s['time_s'] for s in steps[1:])/(len(steps)-1)) if len(steps) > 1 else None
        r['lits_total'] = sum(s['lits'] for s in steps)
        r['nlits_total'] = sum(s['nlits'] for s in steps)
    m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
    if m:
        r['max_rss_mb'] = int(m.group(1))/1024.
    else:
        # BSD/macOS time -l reports bytes
        m = re.search(r'^\s*(\d+)\s+maximum resident set size', text, re.M)
        if m:
            r['max_rss_mb'] = int(m.group(1))/1024.**2
    m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)', text)
    if m:
        parts = [float(p) for p in m.group(1).split(':')]
        secs = 0.
        for p in parts:
            secs = secs*60. + p
        r['wall_s'] = secs
    return r


def mpiexec_cmd(mpiexec):
    """MPI launcher command; OpenMPI binds each rank to one core by default for np<=2, which would leave the
    OpenMP threads of a rank on a single core, so binding is disabled for OpenMPI unless set explicitly"""
    cmd = mpiexec.split()
    if '--bind-to' in cmd:
        return cmd
    try:
        ver = subprocess.run([cmd[0], '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode()
    except OSError:
        ver = ''
    if 'Open MPI' in ver or 'OpenRTE' in ver:
        cmd += ['--bind-to', 'none']
    return cmd


def time_cmd():
    """Prefix that reports the peak memory of a run: GNU time (`-v`) or BSD/macOS time (`-l`), or none"""
    if not os.path.isfile('/usr/bin/time'):
        return []
    for flag in ('-v', '-l'):
        p = subprocess.run(['/usr/bin/time', flag, 'true'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if p.returncode == 0:
            return ['/usr/bin/time', flag]
    return []


def run_case(case, args, src_dir, exe):
    threads = case['threads'] if case['threads'] is not None else args.threads
    run_dir = os.path.join(args.work_dir, args.tag, case['name'])
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)
    write_inputs(case, run_dir, src_dir, args)
    cmd = [exe, 'oft.in', 'oft_in.xml']
    if case['nprocs'] > 1:
        cmd = mpiexec_cmd(args.mpiexec) + ['-np', str(case['nprocs'])] + cmd
    if args.perf and case['nprocs'] == 1:
        cmd = [args.perf_exe, 'record', '-F', str(args.perf_freq), '-g', '-o', 'perf.data', '--'] + cmd
    cmd = time_cmd() + cmd
    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = str(threads)
    env.setdefault('OMPI_ALLOW_RUN_AS_ROOT', '1')
    env.setdefault('OMPI_ALLOW_RUN_AS_ROOT_CONFIRM', '1')
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=run_dir, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=args.timeout)
        out = p.stdout.decode(errors='replace')
        rc = p.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b'').decode(errors='replace') + '\nTIMEOUT\n'
        rc = -999
    with open(os.path.join(run_dir, 'output.log'), 'w') as fid:
        fid.write(out)
    rec = {'name': case['name'], 'tag': args.tag, 'pre': args.pre, 'order': case['order'], 'mg': case['mg'], 'mfnk': case['mfnk'],
           'mode': case['mode'], 'nsteps': case['nsteps'], 'threads': threads, 'nprocs': case['nprocs'],
           'xmhd_opt': list(args.xmhd_opt), 'returncode': rc, 'harness_wall_s': time.time()-t0}
    rec.update(parse_output(out))
    if args.perf and case['nprocs'] == 1 and rc == 0:
        rec['perf_report'] = os.path.join(run_dir, 'perf_report.txt')
        with open(rec['perf_report'], 'w') as fid:
            subprocess.run([args.perf_exe, 'report', '-i', 'perf.data', '--no-children', '--stdio',
                            '--percent-limit', '0.3', '-g', 'none'], cwd=run_dir, stdout=fid, stderr=subprocess.DEVNULL)
        rec['perf_children'] = os.path.join(run_dir, 'perf_children.txt')
        with open(rec['perf_children'], 'w') as fid:
            subprocess.run([args.perf_exe, 'report', '-i', 'perf.data', '--children', '--stdio',
                            '--percent-limit', '0.3', '-g', 'none'], cwd=run_dir, stdout=fid, stderr=subprocess.DEVNULL)
    return rec


CSV_KEYS = ['tag', 'name', 'order', 'mg', 'mfnk', 'mode', 'threads', 'nprocs', 'returncode', 'ncells_global', 'cell_np',
            'mesh_order', 'b_dofs', 'v_dofs', 'n_dofs', 'build_s', 'build_fine_avg_s', 'build_coarse_avg_s', 'nlop_s',
            'matvec_s', 'dot_s', 'jac_nnz_total', 'prof_preupdate_n', 'prof_preupdate_total_s', 'prof_setops_n',
            'prof_setops_total_s', 'prof_nlop_n', 'prof_nlop_total_s', 'qcache_mode', 'qcache_setup_s', 'qcache_mem_mb',
            'jac_posmap', 'posmap_mem_mb', 'posmap_setup_s', 'jcache_coarse_mem_mb', 'step_time_avg_s', 'step_time_avg_excl1_s',
            'lits_total', 'nlits_total', 'run_time_s', 'max_rss_mb', 'wall_s', 'nlop_norm', 'matvec_norm',
            'final_norm_B', 'final_norm_V', 'final_norm_Ti']


def write_results(records, meta, path):
    with open(path, 'w') as fid:
        json.dump({'meta': meta, 'results': records}, fid, indent=1)
    with open(os.path.splitext(path)[0] + '.csv', 'w', newline='') as fid:
        w = csv.DictWriter(fid, fieldnames=CSV_KEYS, extrasaction='ignore')
        w.writeheader()
        for r in records:
            w.writerow(r)


def machine_info(args):
    info = {'hostname_hash': hash(platform.node()) & 0xffff, 'platform': platform.platform(),
            'ncores': os.cpu_count(), 'python': sys.version.split()[0], 'date': time.strftime('%Y-%m-%d %H:%M:%S')}
    try:
        with open('/proc/cpuinfo') as fid:
            for line in fid:
                if line.startswith('model name'):
                    info['cpu'] = line.split(':', 1)[1].strip()
                    break
    except OSError:
        pass
    try:
        info['git_sha'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=HERE).decode().strip()
        info['git_dirty'] = bool(subprocess.check_output(['git', 'status', '--porcelain', '--', '../../src'], cwd=HERE).strip())
    except (OSError, subprocess.CalledProcessError):
        pass
    try:
        with open(os.path.join(args.build_dir, 'CMakeCache.txt')) as fid:
            for line in fid:
                if line.startswith(('CMAKE_Fortran_COMPILER:', 'CMAKE_BUILD_TYPE:', 'OFT_USE_MPI:', 'CMAKE_HOME_DIRECTORY:')):
                    k, v = line.strip().split('=', 1)
                    info[k.split(':')[0]] = v
    except OSError:
        pass
    return info

#------------------------------------------------------------------------------
# Comparison
#------------------------------------------------------------------------------
def compare(path_a, path_b):
    with open(path_a) as fid:
        a = {r['name']: r for r in json.load(fid)['results']}
    with open(path_b) as fid:
        b = {r['name']: r for r in json.load(fid)['results']}
    keys = ['build_s', 'build_fine_avg_s', 'build_coarse_avg_s', 'nlop_s', 'prof_nlop_total_s',
            'matvec_s', 'step_time_avg_excl1_s', 'run_time_s', 'max_rss_mb']
    print('{0:24s} {1:>22s} {2:>12s} {3:>12s} {4:>8s}'.format('case', 'metric', 'A', 'B', 'A/B'))
    for name in a:
        if name not in b:
            continue
        for k in keys:
            va, vb = a[name].get(k), b[name].get(k)
            if va is None or vb is None or vb == 0:
                continue
            print('{0:24s} {1:>22s} {2:12.4e} {3:12.4e} {4:8.2f}'.format(name, k, va, vb, va/vb))
        for k in ('nlop_norm', 'matvec_norm', 'final_norm_B', 'final_norm_V', 'final_norm_Ti'):
            va, vb = a[name].get(k), b[name].get(k)
            if va is not None and vb is not None:
                rel = abs(va-vb)/max(abs(va), 1.e-300)
                print('{0:24s} {1:>22s} {2:12.4e} {3:12.4e} relerr={4:.2e}'.format(name, k, va, vb, rel))
        for k in ('lits_total', 'nlits_total'):
            va, vb = a[name].get(k), b[name].get(k)
            if va is not None and vb is not None:
                print('{0:24s} {1:>22s} {2:12d} {3:12d}'.format(name, k, va, vb))
        if ('steps' in a[name]) and ('steps' in b[name]):
            its_a = [(st['lits'], st['nlits']) for st in a[name]['steps']]
            its_b = [(st['lits'], st['nlits']) for st in b[name]['steps']]
            print('{0:24s} {1:>22s} {2}'.format(name, 'per-step its match', its_a == its_b))

#------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--build-dir', help='OFT build tree (contains tests/physics/MUG/xmhd_bench)')
    parser.add_argument('--exe', help='Explicit path to xmhd_bench executable (overrides --build-dir)')
    parser.add_argument('--tag', default='baseline', help='Label for this run (results/<tag>.json)')
    parser.add_argument('--size', type=float, default=1.0, help='Scale factor on cube mesh resolution')
    parser.add_argument('--matrix', default='default', choices=('default', 'quick'), help='Case matrix')
    parser.add_argument('--cases', default=None, help='Comma separated subset of case names (or regex)')
    parser.add_argument('--threads', type=int, default=os.cpu_count(), help='OMP threads for cases without explicit setting')
    parser.add_argument('--xmhd-opt', action='append', default=[],
                        help='Extra entry for &xmhd_options, eg. "xmhd_qcache=0" (repeatable; applies to all cases)')
    parser.add_argument('--bench-opt', action='append', default=[],
                        help='Extra entry for &xmhd_bench_options (repeatable)')
    parser.add_argument('--pre', default='ilu', choices=('ilu', 'lu'),
                        help='Local solver in block-Jacobi preconditioner/smoother (lu: SuperLU, slow at p>=3)')
    parser.add_argument('--petsc', action='store_true',
                        help='Use the PETSc linear algebra backend (use_petsc=T; OFT must be built with PETSc)')
    parser.add_argument('--prof-reps', default=None,
                        help='Override profile-mode repetitions "builds,nlops_per_build,matvecs_per_nlop" '
                             '(eg. 5,0,0 for builds only; 1,100,0 for NL residual only)')
    parser.add_argument('--debug', type=int, default=0, help='OFT debug level (1 prints per-level build times)')
    parser.add_argument('--perf', action='store_true', help='Record perf profiles (serial cases only)')
    parser.add_argument('--perf-exe', default=shutil.which('perf') or 'perf')
    parser.add_argument('--perf-freq', type=int, default=499)
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--timeout', type=float, default=1800.)
    parser.add_argument('--work-dir', default=os.path.join(HERE, 'runs'))
    parser.add_argument('--out', default=None, help='Results file (default results/<tag>.json)')
    parser.add_argument('--compare', nargs=2, metavar=('A.json', 'B.json'), help='Compare two result files and exit')
    parser.add_argument('--list', action='store_true', help='List cases and exit')
    args = parser.parse_args()

    if args.compare:
        compare(*args.compare)
        return
    ncores = os.cpu_count()
    cases = (case_matrix if args.matrix == 'default' else quick_matrix)(ncores, args.size)
    if args.cases:
        sel = args.cases.split(',')
        cases = [c for c in cases if any(re.fullmatch(s, c['name']) for s in sel)]
    if args.prof_reps:
        for c in cases:
            c['prof_reps'] = args.prof_reps
    if args.list:
        for c in cases:
            print(c['name'], c['order'], 'mg' if c['mg'] else '', 'mfnk' if c['mfnk'] else '', c['mode'], c['nprocs'])
        return
    exe = args.exe or os.path.join(args.build_dir, 'tests', 'physics', 'MUG', 'xmhd_bench')
    exe = os.path.abspath(exe)
    if not os.path.isfile(exe):
        sys.exit('xmhd_bench executable not found: {0} (build OFT with OFT_BUILD_TESTS=ON)'.format(exe))
    src_dir = os.path.abspath(os.path.join(HERE, '..', '..', 'src'))
    out = args.out or os.path.join(HERE, 'results', args.tag + '.json')
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    meta = machine_info(args)
    meta.update({'prof_reps': args.prof_reps, 'pre': args.pre, 'petsc': args.petsc, 'size': args.size, 'tag': args.tag, 'exe': exe, 'xmhd_opt': args.xmhd_opt, 'bench_opt': args.bench_opt})
    records = []
    for case in cases:
        print('Running {0} ...'.format(case['name']), flush=True)
        rec = run_case(case, args, src_dir, exe)
        records.append(rec)
        print('  rc={0} build={1} nlop={2} matvec={3} step={4} rss={5}MB wall={6}'.format(
            rec['returncode'], rec.get('build_s'), rec.get('nlop_s'), rec.get('matvec_s'),
            rec.get('step_time_avg_s'), rec.get('max_rss_mb'), rec.get('wall_s')), flush=True)
        write_results(records, meta, out)
    print('Results written to', out)


if __name__ == '__main__':
    main()

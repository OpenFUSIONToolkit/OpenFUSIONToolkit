#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: equivalence of the xMHD nonlinear residual implementations (xmhd_qcache=0: direct basis
# evaluation, 1: quadrature-point cache in FULL mode, 2: cache in REFERENCE mode). The xmhd_bench driver
# (mode='qcheck') evaluates the residual of the same nonlinear state with each implementation for dt>0,
# dt>0 with fields at the start of the step (back_step), and the dt<0 diagnostics call, and reports the
# differences of every field block, the diagnostic values and the floor flags.
import os
import re
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT

oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
{mesh}/

{mesh_groups}
&xmhd_bench_options
 order={order}
 minlev={minlev}
 mode='qcheck'
 lam={lam}
{bench}/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=T
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_ohmic=T
 xmhd_visc_heat=T
 xmhd_upwind=T
 bbc='bc'
 vbc='all'
 nbc='d'
 tbc='d'
 dt=2.E-7
 eta=1.0
 nu_visc=10.0
 nu_par=100.0
 d_dens=1.0
 kappa_par=1000.0,1000.0
 kappa_perp=1.0,1.0
 nsteps=0
 rst_freq=100000
 nclean=100000
{xmhd}/
"""

# Single-level preconditioner definition (not used by the check, but read by xmhd_run)
oft_in_xml = """<oft>
  <xmhd>
    <pre type="gmres">
      <its>4</its>
      <nrits>4</nrits>
      <pre type="block_jacobi">
        <nlocal>1</nlocal>
        <solver type="ilu"></solver>
      </pre>
    </pre>
  </xmhd>
</oft>
"""

# Differences of field blocks relative to the largest entry of the same equation (B: curl and gradient
# parts, V: all components). Summation order differs between the implementations, so blocks agree to
# ~1E-15. The Grad(H^1) block of B alone is not a good scale: it is a near-cancellation (div(B)~0) whose
# entries are 1E-4 to 1E-7 of the individual quadrature contributions. Implementation errors give
# differences of order one.
block_tol = 1.E-12
equation_groups = {'Bc': 'B', 'Bg': 'B', 'Vx': 'V', 'Vy': 'V', 'Vz': 'V'}
# Differences relative to the largest entry of the whole residual (the norm the nonlinear solver uses)
global_tol = 1.E-13
# Diagnostic integrals and floor flags (no cancellation)
diag_tol = 1.E-12


def cube_mesh(ni, hex_mesh=False, periodic=True, nlevels=1, nbase=1):
    '''Box [-0.5,0.5]^2 x [-1,1] from the built-in cube generator'''
    mesh = " meshname='cube'\n cad_type=92\n nlevels={0}\n nbase={1}\n grid_order=1\n".format(nlevels, nbase)
    per = 'T' if periodic else 'F'
    groups = ("&cube_options\n mesh_type={0}\n ni={1},{1},{2}\n rscale=1.,1.,2.\n shift=-0.5,-0.5,-1.\n"
              " ref_per={3},{3},{3}\n/\n").format(2 if hex_mesh else 1, ni, 2*ni, per)
    return {'mesh': mesh, 'mesh_groups': groups, 'lam': '1.,1.,2.'}


def curved_mesh():
    '''Quadratic (curved) tetrahedral cylinder mesh from the Marklin tests'''
    mesh = " meshname='cyl'\n cad_type=0\n nlevels=1\n nbase=1\n grid_order=2\n fix_boundary=T\n"
    groups = "&native_mesh_options\n filename='../Marklin/cyl_Marklin.h5'\n/\n"
    return {'mesh': mesh, 'mesh_groups': groups, 'lam': '2.,2.,2.'}


def revolved_mesh(cell_type):
    '''Wedge about the x-axis, reflected and made periodic (revolved periodicity), from test_vlag_periodic'''
    sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'fem')))
    from test_vlag_periodic import write_mesh, x_wedge
    mesh_file = 'qcache_revolved_{0}.h5'.format(cell_type)
    write_mesh(mesh_file, x_wedge, cell_type)
    mesh = " meshname='wedge'\n cad_type=0\n nlevels=1\n nbase=1\n grid_order=1\n"
    groups = "&native_mesh_options\n filename='{0}'\n reflect=T\n ref_periodic=T\n/\n".format(mesh_file)
    return {'mesh': mesh, 'mesh_groups': groups, 'lam': '2.,2.,2.'}


def parse_qcheck(output):
    '''Parse XMHD_QCHECK lines into (block results, diagnostic results, flag results, max_rel)'''
    blocks, diags, flags = [], [], []
    max_rel = None
    for line in output.splitlines():
        if not line.startswith('XMHD_QCHECK'):
            continue
        vals = line.split()
        if vals[1] == 'field':
            case = re.search(r'case=(\S+)', line).group(1)
            modes = re.search(r'modes=(\S+)', line).group(1)
            block = re.search(r'block=(\S+)', line).group(1)
            ref_max, max_diff, rel = (float(v) for v in vals[-3:])
            blocks.append((case, modes, block, ref_max, max_diff, rel))
        elif vals[1] == 'diag':
            diags.append((vals[2], int(vals[4]), float(vals[-1])))
        elif vals[1] == 'flags':
            flags.append((vals[2], vals[3], float(vals[-1])))
        elif vals[1] == 'max_rel':
            max_rel = float(vals[2])
    return blocks, diags, flags, max_rel


def qcheck_case(mesh, order, xmhd='', bench='', minlev=-1, nproc=1):
    '''Write the input, run the equivalence check, and validate the reported differences'''
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(order=order, minlev=minlev, xmhd=xmhd, bench=bench, **mesh))
    with open('oft_in_qcache.xml', 'w+') as fid:
        fid.write(oft_in_xml)
    status, output = run_OFT("./xmhd_bench oft.in oft_in_qcache.xml", nproc, 300, return_stdout=True)
    assert status
    blocks, diags, flags, max_rel = parse_qcheck(output)
    # 3 cases x 3 mode pairs for every field block, 7 diagnostics x 3 mode pairs, 3 x 3 flag comparisons
    assert max_rel is not None
    assert len(blocks) > 0 and len(blocks) % 9 == 0
    assert len(diags) == 21
    assert len(flags) == 9
    passed = True
    for case_modes in set((b[0], b[1]) for b in blocks):
        group = [b for b in blocks if (b[0], b[1]) == case_modes]
        scale = max(b[3] for b in group)
        eq_scale = {}
        for b in group:
            eq = equation_groups.get(b[2], b[2])
            eq_scale[eq] = max(eq_scale.get(eq, 0.), b[3])
        for case, modes, block, ref_max, max_diff, rel in group:
            rel_eq = max_diff/max(eq_scale[equation_groups.get(block, block)], 1.E-300)
            if (rel_eq > block_tol) or (max_diff > global_tol*scale):
                print('FAILED: case={0} modes={1} block={2} rel={3:.3E} rel_eq={4:.3E} rel_global={5:.3E}'.format(
                    case, modes, block, rel, rel_eq, max_diff/scale))
                passed = False
    for modes, index, rel in diags:
        if rel > diag_tol:
            print('FAILED: diagnostic {0} {1} rel={2:.3E}'.format(modes, index, rel))
            passed = False
    for case, modes, max_diff in flags:
        if max_diff > diag_tol:
            print('FAILED: floor flags {0} {1} max_diff={2:.3E}'.format(case, modes, max_diff))
            passed = False
    return passed


# Physics option sets (added to the default two-temperature, upwinded, viscous, conducting model above)
physics_cases = {
    'base': ('', ''),
    'single_temp_hall': (" xmhd_hall=T\n te_factor=0.8\n", " two_temp=F\n"),
    'two_temp_hall_equil': (" xmhd_hall=T\n xmhd_therm_equil=T\n", ''),
    'iso_gyro': (" visc_type='iso'\n nu_gyro=1.0\n", ''),
    'kin_gyro_nopar': (" nu_par=-1.0\n nu_gyro=0.5\n", ''),
    'brag_two_temp': (" xmhd_brag=T\n kappa_par=0.5,0.5\n kappa_perp=2.0,2.0\n", ''),
    'brag_single_temp': (" xmhd_brag=T\n kappa_par=0.5,0.5\n kappa_perp=2.0,2.0\n", " two_temp=F\n"),
    'hyper': (" eta_hyper=1.E-3\n d2_dens=1.E-3\n", ''),
    'hyper_centered': (" eta_hyper=1.E-3\n d2_dens=1.E-3\n xmhd_diss_centered=T\n", ''),
    'shaping_eta_temp': (" eta_temp=20.0\n", " shaping=T\n shaping_mix=T\n"),
    'eta_profile': (" eta_temp=20.0\n", " res_prof=T\n"),
    'no_upwind_iso': (" xmhd_upwind=F\n visc_type='iso'\n", ''),
    'no_advec_den': (" xmhd_advec=F\n xmhd_adv_den=F\n xmhd_upwind=F\n", ''),
    'no_temp': (" xmhd_adv_temp=F\n xmhd_jcb=F\n", ''),
    'floors_solid': ('', " den_floor=1.02E19\n temp_floor=20.4\n g_accel=9.8\n qcheck_solid=T\n qcheck_negsrc=T\n"),
}


@pytest.mark.parametrize("physics", list(physics_cases.keys()))
def test_qcache_physics_tet_p2(physics):
    '''Residual implementations agree on straight tets (p=2) for each physics option set'''
    xmhd, bench = physics_cases[physics]
    assert qcheck_case(cube_mesh(2), 2, xmhd=xmhd, bench=bench)


@pytest.mark.parametrize("physics", ['base', 'single_temp_hall', 'hyper', 'shaping_eta_temp', 'floors_solid'])
def test_qcache_physics_tet_p3(physics):
    '''Residual implementations agree on straight tets (p=3) for a subset of the physics option sets'''
    xmhd, bench = physics_cases[physics]
    assert qcheck_case(cube_mesh(2), 3, xmhd=xmhd, bench=bench)


@pytest.mark.parametrize("order", [2, 3, 4])
def test_qcache_tet(order):
    '''Residual implementations agree on straight tets with Hall, hyper-resistivity/diffusivity and shaping'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n nu_gyro=1.0\n"
    assert qcheck_case(cube_mesh(1 if order == 4 else 2), order, xmhd=xmhd, bench=" shaping=T\n")


@pytest.mark.parametrize("order", [2, 3])
def test_qcache_hex(order):
    '''Residual implementations agree on hexes with Hall, hyper-resistivity/diffusivity and shaping'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n nu_gyro=1.0\n"
    assert qcheck_case(cube_mesh(2, hex_mesh=True), order, xmhd=xmhd, bench=" shaping=T\n")


def test_qcache_curved_p3():
    '''Residual implementations agree on a curved (quadratic) tet mesh with shaping, floors and solid cells'''
    bench = " shaping=T\n den_floor=1.02E19\n temp_floor=20.4\n qcheck_solid=T\n qcheck_negsrc=T\n"
    assert qcheck_case(curved_mesh(), 3, xmhd=" xmhd_hall=T\n", bench=bench)


@pytest.mark.parametrize("vbc", ['norm', 'tang'])
def test_qcache_vbc(vbc):
    '''Residual implementations agree on a non-periodic box with normal/tangential velocity and Neumann BCs'''
    xmhd = " vbc='{0}'\n nbc='n'\n tbc='n'\n".format(vbc)
    assert qcheck_case(cube_mesh(2, periodic=False), 2, xmhd=xmhd)


@pytest.mark.parametrize("cell_type", ['tet', 'hex'])
def test_qcache_revolved(cell_type):
    '''Residual implementations agree on a revolved periodic wedge (rotated velocity DOFs)'''
    assert qcheck_case(revolved_mesh(cell_type), 2, xmhd=" vbc='norm'\n")


def test_qcache_pmg():
    '''Residual implementations agree on the finest level of a p-multigrid hierarchy (p=3 down to p=2)'''
    assert qcheck_case(cube_mesh(2), 3, minlev=2)


@pytest.mark.mpi
def test_qcache_mpi2():
    '''Residual implementations agree on a mesh distributed over 2 MPI ranks'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n"
    assert qcheck_case(cube_mesh(2, nlevels=2), 3, xmhd=xmhd, bench=" shaping=T\n qcheck_solid=T\n", nproc=2)

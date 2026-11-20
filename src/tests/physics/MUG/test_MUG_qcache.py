#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: equivalence of the xMHD nonlinear residual and Jacobian implementations (xmhd_qcache=0: direct
# basis evaluation, 1: quadrature-point cache). The xmhd_bench driver (mode='qcheck') evaluates the residual of
# the same nonlinear state with both implementations for dt>0, dt>0 with fields at the start of the step
# (back_step), and the dt<0 diagnostics call, and reports the differences of every field block, the diagnostic
# values and the floor flags. It then builds the Jacobian of every level with xmhd_qcache=0 (searching the CSR
# column lists, the reference), xmhd_qcache=1 with the precomputed CSR position map (the default) and
# xmhd_qcache=0 with the position map, and reports entrywise differences relative to the row maximum and the
# differences of products with random vectors. The *_petsc variants run the same checks with the PETSc linear
# algebra backend (use_petsc=T): the residual blocks on PETSc vectors, and the Jacobians through products with
# random vectors only (the CSR values of PETSc matrices are not accessible).
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
 use_petsc={petsc}
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
# Jacobian entries relative to the largest entry of the same row (summation order differs; curved cells give
# ~3E-12 for H(Curl) curl couplings, whose physical values are rebuilt from cancelling logical components) and
# products with random vectors relative to their norm
jac_row_tol = 1.E-11
jac_matvec_tol = 1.E-12


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


def parse_jcheck(output):
    '''Parse XMHD_JCHECK lines into a list of (level, cfgs, matrix, rel_row, rel_global, rel_matvec)

    matrix is 'native' (entrywise and product differences) or 'generic' (products only, rel_row=rel_global=-1)'''
    results = []
    for line in output.splitlines():
        if line.startswith('XMHD_JCHECK level='):
            level = int(re.search(r'level=\s*(\d+)', line).group(1))
            cfgs = re.search(r'cfgs=(\S+)', line).group(1)
            matrix = re.search(r'matrix=(\S+)', line).group(1)
            rel_row, rel_glob, rel_mv = (float(v) for v in line.split()[-3:])
            results.append((level, cfgs, matrix, rel_row, rel_glob, rel_mv))
    return results


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


def qcheck_case(mesh, order, xmhd='', bench='', minlev=-1, nproc=1, nlevels=1, jac=True, petsc=False):
    '''Write the input, run the equivalence check, and validate the reported differences

    nlevels is the expected number of Jacobian levels (polynomial levels on the finest mesh); jac=False
    skips the Jacobian comparison (memory-heavy configurations); petsc=True uses the PETSc linear algebra
    backend (skipped if OFT is not compiled with PETSc), where the Jacobians are compared through products
    with random vectors only'''
    if not jac:
        bench += " jac_check=F\n"
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(order=order, minlev=minlev, xmhd=xmhd, bench=bench,
                                         petsc='T' if petsc else 'F', **mesh))
    with open('oft_in_qcache.xml', 'w+') as fid:
        fid.write(oft_in_xml)
    status, output = run_OFT("./xmhd_bench oft.in oft_in_qcache.xml", nproc, 300, return_stdout=True)
    assert status
    blocks, diags, flags, max_rel = parse_qcheck(output)
    # 3 cases for every field block, 7 diagnostics (dt<0 case), 3 flag comparisons
    assert max_rel is not None
    assert len(blocks) > 0 and len(blocks) % 3 == 0
    assert len(diags) == 7
    assert len(flags) == 3
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
    # Jacobian: two configurations (cache + position map, direct + position map) against the reference per level
    jres = parse_jcheck(output)
    if jac:
        assert len(jres) == 2*nlevels
        assert len(set(r[0] for r in jres)) == nlevels
    else:
        assert len(jres) == 0
    # The backend must be the requested one (PETSc matrices have no entrywise comparison)
    matrix_exp = 'generic' if petsc else 'native'
    for level, cfgs, matrix, rel_row, rel_glob, rel_mv in jres:
        if matrix != matrix_exp:
            print('FAILED: Jacobian level={0} matrix type {1}, expected {2}'.format(level, matrix, matrix_exp))
            passed = False
        if petsc:
            rel_row = 0.
        if (rel_row > jac_row_tol) or (rel_mv > jac_matvec_tol):
            print('FAILED: Jacobian level={0} cfgs={1} rel_row={2:.3E} rel_matvec={3:.3E}'.format(
                level, cfgs, rel_row, rel_mv))
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
    '''Residual implementations agree on a curved (quadratic) tet mesh with shaping, floors and solid cells

    The Jacobian comparison is skipped here (the p=3 Jacobian of this mesh needs ~8 GB for the check), see
    test_qcache_curved_p2'''
    bench = " shaping=T\n den_floor=1.02E19\n temp_floor=20.4\n qcheck_solid=T\n qcheck_negsrc=T\n"
    assert qcheck_case(curved_mesh(), 3, xmhd=" xmhd_hall=T\n", bench=bench, jac=False)


def test_qcache_curved_p2():
    '''Residual and Jacobian implementations agree on a curved (quadratic) tet mesh with shaping, floors and solid cells'''
    bench = " shaping=T\n den_floor=1.02E19\n temp_floor=20.4\n qcheck_solid=T\n qcheck_negsrc=T\n"
    assert qcheck_case(curved_mesh(), 2, xmhd=" xmhd_hall=T\n", bench=bench)


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
    '''Residual implementations agree on the finest level, Jacobians on all levels of a p-multigrid hierarchy (p=3, 2)'''
    assert qcheck_case(cube_mesh(2), 3, minlev=2, nlevels=2)


def test_qcache_pmg_p4():
    '''Jacobians agree on all levels of a p-multigrid hierarchy (p=4, 3, 2) with Hall, hyper terms and shaping'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n"
    assert qcheck_case(cube_mesh(1), 4, xmhd=xmhd, bench=" shaping=T\n qcheck_solid=T\n", minlev=2, nlevels=3)


def test_qcache_pmg_coarse_quad():
    '''Jacobians agree on all levels of a p-multigrid hierarchy with the experimental coarse-level quadrature'''
    assert qcheck_case(cube_mesh(2), 3, xmhd=" xmhd_coarse_quad=T\n xmhd_hall=T\n", bench=" shaping=T\n",
                       minlev=2, nlevels=2)


@pytest.mark.mpi
def test_qcache_mpi2():
    '''Residual implementations agree on a mesh distributed over 2 MPI ranks'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n"
    assert qcheck_case(cube_mesh(2, nlevels=2), 3, xmhd=xmhd, bench=" shaping=T\n qcheck_solid=T\n", nproc=2)


#============================================================================
# PETSc linear algebra backend: the residual blocks on PETSc vectors and the Jacobians (all levels) through products
# with random vectors. Two-temperature (base) and single-temperature models, Hall, viscosity with parallel and
# gyroviscous terms, floors and solid cells, hexes, revolved periodicity, p-multigrid and 2 MPI ranks.
@pytest.mark.parametrize("physics", ['base', 'single_temp_hall', 'two_temp_hall_equil', 'iso_gyro', 'kin_gyro_nopar',
                                     'floors_solid'])
def test_qcache_petsc_tet_p2(physics):
    '''Residual and Jacobian implementations agree on straight tets (p=2) with PETSc'''
    xmhd, bench = physics_cases[physics]
    assert qcheck_case(cube_mesh(2), 2, xmhd=xmhd, bench=bench, petsc=True)


def test_qcache_petsc_hex_p2():
    '''Residual and Jacobian implementations agree on hexes (p=2) with Hall, hyper terms and shaping with PETSc'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n nu_gyro=1.0\n"
    assert qcheck_case(cube_mesh(2, hex_mesh=True), 2, xmhd=xmhd, bench=" shaping=T\n", petsc=True)


@pytest.mark.parametrize("cell_type", ['tet', 'hex'])
def test_qcache_petsc_revolved(cell_type):
    '''Residual and Jacobian implementations agree on a revolved periodic wedge with PETSc'''
    assert qcheck_case(revolved_mesh(cell_type), 2, xmhd=" vbc='norm'\n", petsc=True)


def test_qcache_petsc_pmg():
    '''Jacobians agree on all levels of a p-multigrid hierarchy (p=3, 2) with PETSc'''
    assert qcheck_case(cube_mesh(2), 3, minlev=2, nlevels=2, petsc=True)


def test_qcache_petsc_pmg_p4():
    '''Jacobians agree on all levels of a p-multigrid hierarchy (p=4, 3, 2) with Hall, hyper terms, shaping and
    solid cells with PETSc'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n"
    assert qcheck_case(cube_mesh(1), 4, xmhd=xmhd, bench=" shaping=T\n qcheck_solid=T\n", minlev=2, nlevels=3,
                       petsc=True)


@pytest.mark.mpi
@pytest.mark.parametrize("pmg", [False, True])
def test_qcache_petsc_mpi2(pmg):
    '''Residual and Jacobian implementations agree on a mesh distributed over 2 MPI ranks with PETSc, on the finest
    level (p=3) or on all levels of a p-multigrid hierarchy (p=3, 2)'''
    xmhd = " xmhd_hall=T\n eta_hyper=1.E-3\n d2_dens=1.E-3\n"
    minlev, nlevels = (3, 2) if pmg else (-1, 1)
    assert qcheck_case(cube_mesh(2, nlevels=2), 3, xmhd=xmhd, bench=" shaping=T\n qcheck_solid=T\n", nproc=2,
                       minlev=minlev, nlevels=nlevels, petsc=True)

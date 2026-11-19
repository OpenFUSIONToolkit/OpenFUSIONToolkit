import os
import sys
import shlex
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
python_dir = os.path.abspath(os.path.join(test_dir, '..', '..', 'python'))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT, run_command

# Input file for the built-in cube mesh
oft_in_cube = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=1
/

&cube_options
 mesh_type={3}
 ni={4},{4},{5}
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per={6},{6},{6}
/

&test_qcache_options
 order={2}
/
"""

# Input file for curved (quadratic) Cubit meshes converted to native format
oft_in_native = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='sphere'
 cad_type=0
 nlevels={1}
 nbase={0}
 grid_order=2
/

&native_mesh_options
 filename='{3}.h5'
/

&test_qcache_options
 order={2}
/
"""

# Maximum allowed normalized error (cache vs direct evaluation)
TOL = 1.E-12


def check_results(order):
    '''Check all errors in qcache.results and that all expected entries are present'''
    results = {}
    with open('qcache.results', 'r') as fid:
        for line in fid:
            name, val = line.split()
            results[name] = float(val)
    expected = ['geom_detw', 'geom_pt']
    spaces = [('lag', ('val', 'grad')), ('h1', ('val', 'grad')), ('hcurl', ('val', 'curl'))]
    suffixes = ['']
    if order > 1:
        suffixes.append('_pm1')
    for suffix in suffixes:
        for space, quants in spaces:
            for mode in ('ref', 'refc'):
                expected.append('{0}{1}_{2}_mat'.format(space, suffix, mode))
                checks = ['eval', 'interp', 'accum']
                if mode != 'refc':
                    checks += ['icell', 'acell']
                for quant in quants:
                    for check in checks:
                        expected.append('{0}{1}_{2}_{3}_{4}'.format(space, suffix, mode, quant, check))
    retval = True
    # Accepted orientation transforms must be exact; spaces with untransformable entity states
    # (eg. H(Curl) faces on tets, order > 1) use the gather path instead of the canonical frame
    for suffix in suffixes:
        for space, _ in spaces:
            canon = results.get('{0}{1}_canonical'.format(space, suffix))
            resid = results.get('{0}{1}_xform_resid'.format(space, suffix))
            nnone = results.get('{0}{1}_xform_none'.format(space, suffix))
            if canon is None or resid is None or nnone is None:
                print('FAILED: missing transform results for "{0}{1}"'.format(space, suffix))
                retval = False
            elif (canon > 0.5) != (nnone < 0.5) or not (resid < TOL):
                print('FAILED: inconsistent transform state for "{0}{1}"'.format(space, suffix))
                retval = False
    for name in expected:
        if name not in results:
            print('FAILED: missing result "{0}"'.format(name))
            retval = False
        elif not (results[name] < TOL):
            print('FAILED: {0} = {1} (tolerance {2})'.format(name, results[name], TOL))
            retval = False
    return retval


def run_cube(order, hex_mesh=False, ni=2, periodic=False, mpi=False):
    nlevels = 2
    nbase = 2
    nproc = 1
    if mpi:
        nbase = 1
        nproc = 2
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_cube.format(nbase, nlevels, order, 2 if hex_mesh else 1, ni, 2*ni,
                                     'T' if periodic else 'F'))
    assert run_OFT("./test_fem_qcache", nproc, 300)
    assert check_results(order)


def run_curved(order, prefix, mpi=False):
    nlevels = 1
    nbase = 1
    nproc = 1
    if mpi:
        nlevels = 2
        nproc = 2
    os.chdir(test_dir)
    convert_cmd = [sys.executable, os.path.join(python_dir, 'OFT_convert_cubit.py'),
                   '--in_file={0}.g'.format(prefix)]
    outs, errs, errcode = run_command(" ".join(shlex.quote(arg) for arg in convert_cmd))
    if errcode != 0:
        print(outs.decode())
        print(errs.decode())
    assert errcode == 0
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_native.format(nbase, nlevels, order, prefix))
    assert run_OFT("./test_fem_qcache", nproc, 300)
    assert check_results(order)


#============================================================================
# Straight tetrahedra (cube)
@pytest.mark.parametrize("order", (1, 2, 3, 4))
def test_tet(order):
    run_cube(order)


@pytest.mark.parametrize("order", (2, 3))
def test_tet_periodic(order):
    run_cube(order, periodic=True)


@pytest.mark.mpi
@pytest.mark.parametrize("order", (2, 3))
def test_tet_mpi(order):
    run_cube(order, mpi=True)


#============================================================================
# Hexahedra (cube)
@pytest.mark.parametrize("order", (1, 2, 3))
def test_hex(order):
    run_cube(order, hex_mesh=True)


@pytest.mark.mpi
def test_hex_mpi():
    run_cube(2, hex_mesh=True, mpi=True)


#============================================================================
# Curved (quadratic) tetrahedra and hexahedra (sphere)
@pytest.mark.parametrize("order", (1, 2, 3, 4))
def test_curved_tet(order):
    run_curved(order, 'sphere_tet10_test')


@pytest.mark.parametrize("order", (1, 2, 3))
def test_curved_hex(order):
    run_curved(order, 'sphere_hex27_test')


@pytest.mark.mpi
def test_curved_tet_mpi():
    run_curved(3, 'sphere_tet10_test', mpi=True)

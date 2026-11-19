#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Shared helpers for ThickCurr tests: structured box mesh generation, input templates, runner, and analytic
# decay times for the test geometries
import os
import sys
import json
import math
import numpy
import h5py
from scipy.optimize import brentq
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing
from oft_testing import run_OFT

MU0 = 4.E-7*math.pi


def write_box_mesh(filename, n, slab=None):
    '''Write a native tetrahedral mesh of the unit cube with n=(nx,ny,nz) divisions

    Each hexahedral block is split into 6 tetrahedra using the Kuhn triangulation, which is conforming
    between neighboring blocks. If `slab=(z1,z2)` is given, cells with centroids in z1 < z < z2 are placed
    in region 2 and all others in region 1; otherwise all cells are in region 1.
    '''
    nx, ny, nz = n
    x = numpy.linspace(0.0, 1.0, nx+1)
    y = numpy.linspace(0.0, 1.0, ny+1)
    z = numpy.linspace(0.0, 1.0, nz+1)
    r = numpy.array([[x[i], y[j], z[k]] for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)])

    def pind(i, j, k):
        return i + j*(nx+1) + k*(nx+1)*(ny+1)
    perms = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
    lc = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                for perm in perms:
                    ind = [i, j, k]
                    cell = [pind(*ind)]
                    for axis in perm:
                        ind[axis] += 1
                        cell.append(pind(*ind))
                    lc.append(cell)
    lc = numpy.array(lc, dtype=numpy.int32)
    # Ensure positive orientation
    for cell in lc:
        dr = r[cell[1:], :] - r[cell[0], :]
        if numpy.linalg.det(dr) < 0.0:
            cell[[2, 3]] = cell[[3, 2]]
    reg = numpy.ones((lc.shape[0],), dtype=numpy.int32)
    if slab is not None:
        zc = r[lc, 2].mean(axis=1)
        reg[(zc > slab[0]) & (zc < slab[1])] = 2
    with h5py.File(os.path.join(test_dir, filename), 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        dset = h5_file.create_dataset('mesh/LC', data=lc+1, dtype='i4')
        dset.attrs["TYPE"] = "tet".encode('ascii')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')


PRE_XML = {
    # Direct solve (LU) of the full system
    'lu': """<oft>
  <thickcurr>
    <pre type="block_jacobi">
      <nlocal>1</nlocal>
      <solver type="lu"></solver>
    </pre>
  </thickcurr>
</oft>
""",
    # p-multigrid with GMRES + block Jacobi smoothing and a direct solve on the coarse (linear) level
    'mg': """<oft>
  <thickcurr>
    <pre type="mg">
      <smoother direction="both">
        <solver type="gmres">
          <its>2</its>
          <nrits>2</nrits>
          <pre type="block_jacobi">
            <nlocal>4</nlocal>
            <solver type="lu"></solver>
          </pre>
        </solver>
      </smoother>
      <coarse>
        <solver type="gmres">
          <its>1</its>
          <nrits>1</nrits>
          <pre type="block_jacobi">
            <nlocal>1</nlocal>
            <solver type="lu"></solver>
          </pre>
        </solver>
      </coarse>
    </pre>
  </thickcurr>
</oft>
"""
}


def run_driver(config, nproc=1, timeout=180):
    '''Run the Python ThickCurr driver (`thickcurr_driver.py`) with the given configuration, using MPI if
    `nproc > 1`, and return the results written by the driver (None on failure)'''
    # Coverage configuration path is relative to the root test directory (not the current directory, which may be
    # changed by other tests), make absolute for the driver process
    cov_config = os.environ.get('COVERAGE_PROCESS_START', None)
    if (cov_config is not None) and (not os.path.isabs(cov_config)):
        root_test_dir = os.path.dirname(os.path.abspath(oft_testing.__file__))
        os.environ['COVERAGE_PROCESS_START'] = os.path.join(root_test_dir, cov_config)
    os.chdir(test_dir)
    for filename in ('thickcurr_results.json', 'thickcurr.hist', 'floops.hist'):
        if os.path.exists(filename):
            os.remove(filename)
    pre = config.get('pre', 'lu')
    if pre is not None:
        with open('oft_in.xml', 'w+') as fid:
            fid.write(PRE_XML[pre])
        config['xml_filename'] = 'oft_in.xml'
    with open('thickcurr_config.json', 'w+') as fid:
        json.dump(config, fid)
    command = "{0} {1} thickcurr_config.json".format(sys.executable, os.path.join(test_dir, 'thickcurr_driver.py'))
    if not run_OFT(command, nproc, timeout):
        return None
    if not os.path.exists('thickcurr_results.json'):
        print("FAILED: No results file written by driver")
        return None
    with open('thickcurr_results.json', 'r') as fid:
        return json.load(fid)


def check_value(name, value, expected, rtol):
    '''Compare a value to its expected value with a relative tolerance'''
    rel_err = abs(value-expected)/abs(expected)
    if (not math.isfinite(value)) or (rel_err > rtol):
        print("FAILED: {0} outside tolerance".format(name))
        print("  Expected = {0:.8E}".format(expected))
        print("  Actual =   {0:.8E}".format(value))
        print("  Rel err =  {0:.3E} (tol = {1:.1E})".format(rel_err, rtol))
        return False
    return True


def box_tau(eta, mode_numbers):
    '''Decay time of the diffusion eigenmode with indices (l,m,n) in a perfectly conducting unit cube'''
    k2 = math.pi**2*sum(m**2 for m in mode_numbers)
    return MU0/(eta*k2)


def slab_tau(eta, z1, z2):
    '''Decay time of the slowest mode of a conducting slab z1 < z < z2 between vacuum layers in a
    perfectly conducting unit cube

    The slowest mode varies as cos(pi*x) (or cos(pi*y)), so the horizontal wavenumber is kappa = pi. Matching
    a potential field in each vacuum layer (with B_z = 0 on the walls) to the even poloidal mode in the slab,
    h(z) = cos(s*(z-zc)), gives s*tan(s*d/2) = kappa*coth(kappa*z1) for a slab centered in the box. The decay
    rate is then eta*(s^2 + kappa^2)/mu0.
    '''
    d = z2 - z1
    if abs(z1 - (1.0-z2)) > 1.E-12:
        raise ValueError("Slab must be centered in box")
    kappa = math.pi
    rhs = kappa/math.tanh(kappa*z1)
    s = brentq(lambda s: s*math.tan(s*d/2.0) - rhs, 1.E-6, math.pi/d*(1.0-1.E-10))
    return MU0/(eta*(s**2 + kappa**2))

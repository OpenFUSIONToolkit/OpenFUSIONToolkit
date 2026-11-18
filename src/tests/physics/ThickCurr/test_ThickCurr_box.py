#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThickCurr tests: resistive diffusion in a fully conducting unit cube with perfectly conducting walls,
# compared with the analytic eigenmodes of the curl-curl operator (B.n = 0 on the walls). All cases are run
# through the Python interface (see thickcurr_driver.py), including the MPI cases.
import math
import pytest
from thickcurr_testing import run_driver, write_box_mesh, check_value, box_tau

ETA = 1.256637E-6  # Resistivity [Ohm-m], gives a magnetic diffusivity of ~1 m^2/s
MESH_FILE = 'thickcurr_box.h5'
# Range of maximum linear iterations per step with the XML p-multigrid definition (reference: 14 in serial, 20 on
# two processes). The default (non-XML) smoothers require only ~3 iterations on this small case, so the lower bound
# detects the XML definition being ignored.
MG_LITS_RANGE = (8, 20)
MG_LITS_RANGE_MPI = (8, 30)


def box_signals():
    '''Analytic sensor signals for A = sin(pi*x)*sin(pi*y) z_hat (see thickcurr_driver.py)'''
    pi = math.pi

    def bfield(x, y):
        return (pi*math.sin(pi*x)*math.cos(pi*y), -pi*math.cos(pi*x)*math.sin(pi*y))
    b1 = bfield(0.3, 0.2)
    b2 = bfield(0.7, 0.35)
    b3 = bfield(0.15, 0.6)
    floop = math.sin(0.3*pi)*(math.sin(0.5*pi)-math.sin(0.1*pi))*0.6
    return [b1[0], b2[1], (b3[0]+b3[1])/math.sqrt(2.0), floop]


def run_decay(order, nproc=1, pre='lu', eig_shift=None, in_memory=False, pmg=False, lits_range=None):
    write_box_mesh(MESH_FILE, (4, 4, 4))
    # Crank-Nicolson for second-order time accuracy (no stiff pseudo-vacuum modes here)
    config = {'case': 'decay', 'mesh_file': MESH_FILE, 'eta': [ETA], 'order': order, 'dt': 1.E-3, 'theta': 0.5,
              'nsteps': 20, 'plot_freq': 10, 'pre': pre, 'eig_shift': eig_shift, 'in_memory': in_memory,
              'pmg': pmg}
    results = run_driver(config, nproc=nproc)
    if results is None:
        return False
    tau = box_tau(ETA, (1, 1, 0))
    t_final = results['t_final']
    retval = check_value("t_final", t_final, 2.E-2, 1.E-10)
    # Decay rate from magnetic energy
    e0, e1 = results['energy']
    retval &= check_value("Decay rate", -math.log(e1/e0)/(2.0*t_final), 1.0/tau, 2.E-2)
    # Initial sensor signals compared with analytic values. Point (Mirnov) probes sample B = curl(A) directly
    # and include pointwise discretization error on this coarse mesh, while the flux loop (integral of A.dl) is
    # much more accurate.
    rtols = [8.E-2]*3 + [1.E-2]
    for i, (val, val_exp, rtol) in enumerate(zip(results['signals0'], box_signals(), rtols)):
        retval &= check_value("Initial signal {0}".format(i+1), val, val_exp, rtol)
    # Mode shape is conserved, so all signals decay at the same rate (point probes also see small local
    # contributions from faster decaying modes introduced by projection error)
    rtols = [5.E-2]*3 + [1.E-2]
    for i, (val0, val1, rtol) in enumerate(zip(results['signals0'], results['signals1'], rtols)):
        retval &= check_value("Signal decay {0}".format(i+1), val1/val0, math.exp(-t_final/tau), rtol)
    # Final sensor history entry matches final signals
    retval &= check_value("History time", results['hist_t1'], t_final, 1.E-10)
    for i, (val, val_exp) in enumerate(zip(results['hist_signals1'], results['signals1'])):
        retval &= check_value("History signal {0}".format(i+1), val, val_exp, 1.E-8)
    # Check that the expected preconditioner is used (maximum number of linear iterations per step)
    if lits_range is not None:
        if (results['max_lits'] < lits_range[0]) or (results['max_lits'] > lits_range[1]):
            print("FAILED: Linear iterations outside expected range ({0} not in [{1}, {2}])".format(
                results['max_lits'], *lits_range))
            retval = False
    return retval


def run_eigs(order, nproc=1, dt=5.E-2, eig_shift=None):
    write_box_mesh(MESH_FILE, (4, 4, 4))
    config = {'case': 'eigs', 'mesh_file': MESH_FILE, 'eta': [ETA], 'order': order, 'dt': dt, 'nev': 5,
              'eig_shift': eig_shift}
    results = run_driver(config, nproc=nproc)
    if results is None:
        return False
    # Lowest modes: (1,1,0) family (3-fold degenerate), then (1,1,1) (2-fold, TE and TM)
    tau_exp = [box_tau(ETA, (1, 1, 0))]*3 + [box_tau(ETA, (1, 1, 1))]*2
    retval = True
    for i, (val, val_exp) in enumerate(zip(results['tau'], tau_exp)):
        retval &= check_value("tau_{0}".format(i+1), val, val_exp, 2.E-2)
    return retval


@pytest.mark.parametrize("order", (2, 3))
@pytest.mark.coverage
def test_decay(order):
    '''Time advance of the slowest diffusion eigenmode of a conducting cube, checking the decay rate and sensor
    signals (point Mirnov probes and a flux loop) against analytic values.'''
    assert run_decay(order)


@pytest.mark.coverage
@pytest.mark.mpi
def test_decay_mpi():
    '''Time advance of a conducting cube eigenmode on two MPI processes, driven from Python.'''
    assert run_decay(2, nproc=2)


@pytest.mark.coverage
def test_decay_iterative():
    '''Time advance of a conducting cube eigenmode using the default (diagonal) preconditioner.'''
    assert run_decay(2, pre=None)


@pytest.mark.coverage
def test_decay_in_memory():
    '''Time advance of a conducting cube eigenmode with the mesh passed from Python as arrays.'''
    assert run_decay(2, in_memory=True)


@pytest.mark.coverage
def test_decay_pmg():
    '''Time advance of a conducting cube eigenmode using a p-multigrid preconditioner (default smoothers).'''
    assert run_decay(2, pre=None, pmg=True)


@pytest.mark.mpi
def test_decay_pmg_mpi():
    '''Time advance of a conducting cube eigenmode using a p-multigrid preconditioner on two MPI processes.'''
    assert run_decay(2, nproc=2, pre=None, pmg=True)


@pytest.mark.coverage
def test_decay_pmg_xml():
    '''Time advance of a conducting cube eigenmode using a p-multigrid preconditioner defined in XML, checking
    that the XML definition is used (iteration count). The XML is passed directly to ThickCurr from Python, so
    no global XML is loaded.'''
    assert run_decay(2, pre='mg', pmg=True, lits_range=MG_LITS_RANGE)


@pytest.mark.mpi
def test_decay_pmg_xml_mpi():
    '''Time advance of a conducting cube eigenmode using a p-multigrid preconditioner defined in XML on two MPI
    processes.'''
    assert run_decay(2, nproc=2, pre='mg', pmg=True, lits_range=MG_LITS_RANGE_MPI)


@pytest.mark.coverage
def test_decay_after_eigs():
    '''Time advance of a conducting cube eigenmode after an eigenvalue solve with a different shift, checking
    that the time advance operators are restored.'''
    assert run_decay(2, eig_shift=5.E-2)


@pytest.mark.parametrize("order", (2, 3))
@pytest.mark.coverage
def test_eigs(order):
    '''Slowest decay times of a conducting cube from the eigenvalue solver, compared with the analytic diffusion
    eigenvalues including their degeneracy.'''
    assert run_eigs(order)


@pytest.mark.coverage
@pytest.mark.mpi
def test_eigs_mpi():
    '''Slowest decay times of a conducting cube on two MPI processes, driven from Python.'''
    assert run_eigs(2, nproc=2)


@pytest.mark.coverage
def test_eigs_shift():
    '''Slowest decay times of a conducting cube using a shift for the eigenvalue solve that differs from the
    time step (requires rebuilding the LHS matrices and updating the preconditioner).'''
    assert run_eigs(2, dt=1.E-3, eig_shift=5.E-2)


@pytest.mark.coverage
def test_recreate():
    '''Create, destroy, and recreate ThickCurr models (different meshes and solver settings) in the same process,
    checking that state from the first model (XML preconditioner, plotting operators, and restart files) does not
    leak into the second model.'''
    write_box_mesh(MESH_FILE, (4, 4, 4))
    write_box_mesh('thickcurr_box3.h5', (3, 3, 3))
    config = {'case': 'recreate', 'mesh_file': MESH_FILE, 'mesh_file2': 'thickcurr_box3.h5', 'eta': [ETA],
              'order': 2, 'dt': 5.E-2}
    results = run_driver(config)
    assert results is not None
    tau_exp = box_tau(ETA, (1, 1, 0))
    retval = True
    for name in ('tau1', 'tau2'):
        for i, val in enumerate(results[name]):
            retval &= check_value("{0}_{1}".format(name, i+1), val, tau_exp, 5.E-2)
    # Second model has no XML preconditioner, so the default (diagonal) preconditioner should be used
    if results['max_lits'] <= 1:
        print("FAILED: Preconditioner from previous model used ({0} iterations)".format(results['max_lits']))
        retval = False
    if results['stale_rst']:
        print("FAILED: Restart file from previous run not removed")
        retval = False
    assert retval

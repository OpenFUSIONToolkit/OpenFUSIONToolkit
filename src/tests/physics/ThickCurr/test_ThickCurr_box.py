#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThickCurr tests: resistive diffusion in a fully conducting unit cube with perfectly conducting walls,
# compared with the analytic eigenmodes of the curl-curl operator (B.n = 0 on the walls)
import pytest
from thickcurr_testing import thickcurr_setup, write_box_mesh, read_results, check_value, box_tau

ETA = 1.256637E-6  # Resistivity [Ohm-m], gives a magnetic diffusivity of ~1 m^2/s
MESH_FILE = 'thickcurr_box.h5'


def write_mesh(nproc):
    # MPI runs add a level that distributes the same mesh over the ranks
    write_box_mesh(MESH_FILE, (4, 4, 4))
    return 2 if nproc > 1 else 1


def run_decay(order, nproc=1, direct=True):
    nlevels = write_mesh(nproc)
    # Crank-Nicolson for second-order time accuracy (no stiff pseudo-vacuum modes here)
    assert thickcurr_setup('decay', MESH_FILE, 1, ETA, order=order, nlevels=nlevels, dt=1.E-3, theta=0.5,
                           nsteps=20, rst_freq=10, direct=direct, nproc=nproc)
    rate, shape_err = read_results()
    rate_exp = 1.0/box_tau(ETA, (1, 1, 0))
    retval = check_value("Decay rate", rate, rate_exp, 2.E-2)
    if shape_err > 1.E-2:
        print("FAILED: Mode shape changed during decay, err = {0:.3E}".format(shape_err))
        retval = False
    return retval


def run_eigs(order, nproc=1):
    nlevels = write_mesh(nproc)
    assert thickcurr_setup('eigs', MESH_FILE, 1, ETA, order=order, nlevels=nlevels, dt=5.E-2, nev=5,
                           nproc=nproc)
    tau = read_results()
    # Lowest modes: (1,1,0) family (3-fold degenerate), then (1,1,1) (2-fold, TE and TM)
    tau_exp = [box_tau(ETA, (1, 1, 0))]*3 + [box_tau(ETA, (1, 1, 1))]*2
    retval = True
    for i, (val, val_exp) in enumerate(zip(tau, tau_exp)):
        retval &= check_value("tau_{0}".format(i+1), val, val_exp, 2.E-2)
    return retval


@pytest.mark.parametrize("order", (2, 3))
def test_decay(order):
    '''Time advance (thickcurr_run) of the slowest diffusion eigenmode of a conducting cube, checking the decay
    rate against the analytic value and the conservation of the mode shape.'''
    assert run_decay(order)


@pytest.mark.mpi
def test_decay_mpi():
    '''Time advance of a conducting cube eigenmode on two MPI processes.'''
    assert run_decay(2, nproc=2)


@pytest.mark.mpi
def test_eigs_mpi():
    '''Slowest decay times of a conducting cube on two MPI processes.'''
    assert run_eigs(2, nproc=2)


def test_decay_iterative():
    '''Time advance of a conducting cube eigenmode using the default (diagonal) preconditioner.'''
    assert run_decay(2, direct=False)


@pytest.mark.parametrize("order", (2, 3))
def test_eigs(order):
    '''Slowest decay times of a conducting cube from the eigenvalue solver (thickcurr_eigs), compared with the
    analytic diffusion eigenvalues including their degeneracy.'''
    assert run_eigs(order)

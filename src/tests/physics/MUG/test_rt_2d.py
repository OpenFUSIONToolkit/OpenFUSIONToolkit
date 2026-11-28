#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: gravity and the Rayleigh-Taylor instability with the 2D xMHD solver. An exponentially stratified
# plasma between two conducting, free-slip walls (z = +-1/2) in a box periodic in x, with gravity along -z and a
# constant-Alfven-speed field along x, compared with the incompressible growth rate of the m=1 mode
# (S. Chandrasekhar, Hydrodynamic and Hydromagnetic Stability, 1961, sections 92 and 97)
import os
import sys
import math
import numpy as np
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT
from OpenFUSIONToolkit.io import histfile

# Input file template, filled in by the setup function
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='test'
 cad_type=92
 nlevels={nlevels}
 nbase=1
 grid_order=2
 fix_boundary=T
/

&cube_options
 mesh_type=1
 ni=4,4,0
 rscale=1.,1.,0.
 shift=-0.5,-0.5,0.
 ref_per=T,F,F
/

&rt_2d_options
 order=3
 linear={linear}
 nsteps={nsteps}
 n_0=1.E19
 T_0=100.
 L_n={L_n}
 g={g}
 b_frac={b_frac}
 V1={V1}
 L_x=1.
 L_z=1.
 dt={dt}
 nu=1.
 eta=1.
 lin_tol=1.E-10
 nl_tol=1.E-8
/
"""

G_ACCEL = 1.E8 # Gravitational acceleration [m/s^2]
L_N = 0.5 # Density scale length [m]
DT = 2.E-6 # Time step [s], gamma*dt = 0.025 for B_0 = 0
N_0 = 1.E19 # Density at z=0 [m^-3]
M_I = 1.672621777E-27 # Ion (proton) mass [kg]


def gamma_theory(b_frac, L_n=L_N, g=G_ACCEL):
    '''Incompressible, ideal growth rate of the m=1 mode with k = 2*pi in the unit box, for B_0 = b_frac*B_c'''
    k = 2.0*math.pi
    K2 = k**2 + math.pi**2 + 1.0/(4.0*L_n**2)
    return k*math.sqrt((g/L_n)/K2*(1.0-b_frac**2))


def rt_2d_setup(linear, nsteps, b_frac=0.0, L_n=L_N, V1=1.0, nproc=1):
    # The level above nbase is the base mesh distributed over the processes (a transfer level, not a refinement),
    # so MPI runs reproduce the serial result
    nlevels = 2 if nproc > 1 else 1
    os.chdir(test_dir)
    if os.path.exists('xmhd_2d.hist'):
        os.remove('xmhd_2d.hist')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels=nlevels, linear=('T' if linear else 'F'), nsteps=nsteps,
                                         L_n=L_n, g=G_ACCEL, b_frac=b_frac, V1=V1, dt=DT))
    # Block-Jacobi/LU preconditioner (as in the Taylor-Green tests); the default diagonal one does not converge
    return run_OFT("./test_rt_2d oft.in oft_in_taylor_green.xml", nproc, 300)


def kinetic_energy():
    hist = histfile('xmhd_2d.hist')
    return np.asarray(hist['time'], dtype=float), np.asarray(hist['ven'], dtype=float)


def growth_rate(time, ven, istart, iend):
    '''Growth rate of the velocity, (1/2) d ln(ven)/dt, from a least-squares fit over steps [istart, iend)'''
    return np.polyfit(time[istart:iend], 0.5*np.log(ven[istart:iend]), 1)[0]


def check_rate(name, gamma, gamma_ref, tol):
    print('{0}: gamma = {1:.5E}, reference = {2:.5E}, ratio = {3:.5F}'.format(name, gamma, gamma_ref,
                                                                           gamma/gamma_ref))
    return math.isfinite(gamma) and abs(gamma/gamma_ref-1.0) < tol

#============================================================================
def test_stable_at_rest():
    '''Non-linear run of a stable stratification (density decreasing upward) with no seed: the plasma must stay
    at rest. The kinetic energy after 40 steps is compared with that of free fall over the same time, so a wrong
    sign or magnitude of gravity (relative to the hydrostatic pressure gradient) fails the test.'''
    nsteps = 40
    assert rt_2d_setup(False, nsteps, L_n=-L_N, V1=0.0)
    time, ven = kinetic_energy()
    ke_free_fall = 0.5*M_I*N_0*(G_ACCEL*DT*nsteps)**2
    # The remaining motion is the adjustment of the projected equilibrium to discrete balance
    print('Max kinetic energy / free fall = {0:.3E}'.format(ven.max()/ke_free_fall))
    assert np.all(np.isfinite(ven))
    assert ven.max() < 1.E-5*ke_free_fall


@pytest.mark.linear
@pytest.mark.parametrize("b_frac", (0.0, 0.5))
def test_lin_growth(b_frac):
    '''Linearized Rayleigh-Taylor growth of the m=1 mode with no field and with B_0 = B_c/2 along x, compared
    with the incompressible ideal growth rate. The velocity-only seed grows as cosh(gamma*t), so the rate is fit
    once gamma*t > 3; the remaining differences are compressibility, viscosity and resistivity (<1%).'''
    nsteps = 200
    assert rt_2d_setup(True, nsteps, b_frac=b_frac)
    time, ven = kinetic_energy()
    assert check_rate('Linear', growth_rate(time, ven, 150, nsteps), gamma_theory(b_frac), 0.03)


@pytest.mark.mpi
@pytest.mark.linear
def test_lin_growth_mpi():
    '''Linearized Rayleigh-Taylor growth with no field on two MPI processes, compared with the incompressible
    ideal growth rate.'''
    nsteps = 200
    assert rt_2d_setup(True, nsteps, nproc=2)
    time, ven = kinetic_energy()
    assert check_rate('Linear (MPI)', growth_rate(time, ven, 150, nsteps), gamma_theory(0.0), 0.03)


def test_nl_growth():
    '''Non-linear Rayleigh-Taylor growth with no field, while the perturbation is small, compared with a
    linearized run from the same seed over the same window. Both start from the velocity-only seed, so the
    cosh(gamma*t) transient cancels in the comparison. The 10 m/s seed is far above the kinetic energy from the
    discrete imbalance of the projected equilibrium and keeps k*xi < 0.03.'''
    nsteps = 100
    assert rt_2d_setup(True, nsteps, V1=10.0)
    time, ven = kinetic_energy()
    gamma_lin = growth_rate(time, ven, 60, nsteps)
    assert rt_2d_setup(False, nsteps, V1=10.0)
    time, ven = kinetic_energy()
    gamma_nl = growth_rate(time, ven, 60, nsteps)
    assert check_rate('Non-linear', gamma_nl, gamma_lin, 0.02)

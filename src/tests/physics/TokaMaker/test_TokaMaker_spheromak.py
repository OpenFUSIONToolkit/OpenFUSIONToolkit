#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: fixed-boundary spheromak equilibrium (analytic Bessel-function eigenmode)
import pytest
import numpy as np
from scipy.special import jv, jn_zeros
from tokamaker_testing import mp_run
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain


def run_sph_case(mesh_resolution,fe_order,mp_q):
    def spheromak_psi(r_grid,z_grid,a,h):
        gamma_11 = jn_zeros(1,1)[0]*r_grid/a
        x_01 = jn_zeros(0,1)[0]
        norm = x_01*jv(1,x_01)
        return gamma_11*jv(1,gamma_11)*np.sin(np.pi*z_grid/h)/norm
    gs_mesh = gs_Domain()
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.add_rectangle(0.5,0.5,1.0,1.0,'plasma')
    mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc)
    mygs.settings.free_boundary = False
    mygs.setup(order=fe_order)
    mygs.p_scale=0.0
    ffp_prof={
        'type': 'linterp',
        'x': [0.0,1.0],
        'y': [1.0,0.0],
    }
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof={'type': 'flat'})
    mygs.settings.nl_tol = 1.E-12
    mygs.settings.maxits = 100
    mygs.urf = 0.0
    mygs.update_settings()
    mygs.init_psi()
    try:
        mygs.solve()
    except ValueError:
        mp_q.put(None)
        return
    psi_TM = mygs.get_psi(False)
    psi_eig_TM = spheromak_psi(mygs.r[:,0], mygs.r[:,1],1.0,1.0)
    # The eigenmode amplitude is arbitrary, so project onto the analytic mode before comparing
    psi_TM /= psi_TM.dot(psi_eig_TM)/psi_eig_TM.dot(psi_eig_TM)
    psi_err = np.linalg.norm(psi_TM-psi_eig_TM)/np.linalg.norm(psi_eig_TM)
    mp_q.put([psi_err])
    oftpy_dump_cov()


def validate_sph(results,psi_err_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs(results[0]) > abs(psi_err_exp)*1.1:
        print("FAILED: psi error too high!")
        print("  Expected = {0}".format(psi_err_exp))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    return test_result


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h1(order):
    '''Fixed-boundary spheromak equilibrium (force-free Taylor state, J. B. Taylor, Rev. Mod. Phys. 58, 741 (1986)) in
    a unit cylinder, whose flux function is the analytic Bessel-function eigenmode r J1(x11 r) sin(pi z): relative
    psi error at the base mesh resolution for element orders 2-4.'''
    errs = [2.039674417912789e-05, 5.103597862537552e-07, 8.088772274705608e-09]
    results = mp_run(run_sph_case,(0.05,order))
    assert validate_sph(results,errs[order-2])
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h2(order):
    '''Fixed-boundary spheromak equilibrium (force-free Taylor state, J. B. Taylor, Rev. Mod. Phys. 58, 741 (1986)) in
    a unit cylinder, whose flux function is the analytic Bessel-function eigenmode r J1(x11 r) sin(pi z): relative
    psi error at half the base mesh resolution for element orders 2-4.'''
    errs = [2.5203856661960034e-06, 3.279268054674832e-08, 2.5185712724779513e-10]
    results = mp_run(run_sph_case,(0.05/2.0,order))
    assert validate_sph(results,errs[order-2])
@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h3(order):
    '''Fixed-boundary spheromak equilibrium (force-free Taylor state, J. B. Taylor, Rev. Mod. Phys. 58, 741 (1986)) in
    a unit cylinder, whose flux function is the analytic Bessel-function eigenmode r J1(x11 r) sin(pi z): relative
    psi error at a quarter of the base mesh resolution for element orders 2-4.'''
    errs = [3.257155111957006e-07, 2.090369020180253e-09, 8.601148342547016e-12]
    results = mp_run(run_sph_case,(0.05/4.0,order))
    assert validate_sph(results,errs[order-2])

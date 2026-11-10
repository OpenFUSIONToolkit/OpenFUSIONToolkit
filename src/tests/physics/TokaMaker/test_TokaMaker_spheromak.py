#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: fixed-boundary spheromak equilibrium (analytic Bessel-function eigenmode)
import pytest
import numpy as np
from scipy.special import jv, jn_zeros
from tokamaker_testing import mp_run, validate_errors
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain


def spheromak_psi(r_grid,z_grid,a,h):
    '''Lowest Grad-Shafranov eigenmode in a cylinder of radius `a` and height `h`, normalized to unit peak'''
    gamma_11 = jn_zeros(1,1)[0]*r_grid/a
    x_01 = jn_zeros(0,1)[0]
    norm = x_01*jv(1,x_01)
    return gamma_11*jv(1,gamma_11)*np.sin(np.pi*z_grid/h)/norm


def run_sph_case(mesh_resolution,fe_order,mp_q):
    gs_mesh = gs_Domain()
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.add_rectangle(0.5,0.5,1.0,1.0,'plasma')
    mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.settings.free_boundary = False
    mygs.setup_model(mesh_pts, mesh_lc, order=fe_order)
    # Force-free with FF' linear in psi, so the equilibrium is an eigenmode of Delta*
    mygs.p_scale=0.0
    ffp_prof={
        'type': 'linterp',
        'x': [0.0,1.0],
        'y': [1.0,0.0],
    }
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof={'type': 'flat'})
    mygs.settings.nl_tol = 1.E-12
    mygs.settings.maxits = 100
    # No under-relaxation of the nonlinear iteration
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
    return validate_errors(results,(psi_err_exp,),('psi',))


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h1(order):
    '''Fixed-boundary force-free spheromak in a square cylindrical domain on the coarsest mesh: the relative flux
    error against the analytic Bessel-function eigenmode must not exceed its recorded value.

    Reference: C. Hansen et al., Comput. Phys. Commun. 298, 109111 (2024), Sec. 4.1.1.
    '''
    errs = [2.039674417912789e-05, 5.103597862537552e-07, 8.088772274705608e-09]
    results = mp_run(run_sph_case,(0.05,order))
    assert validate_sph(results,errs[order-2])


@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h2(order):
    '''Spheromak eigenmode with the mesh size halved, checking convergence with resolution and element order'''
    errs = [2.5203856661960034e-06, 3.279268054674832e-08, 2.5185712724779513e-10]
    results = mp_run(run_sph_case,(0.05/2.0,order))
    assert validate_sph(results,errs[order-2])


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h3(order):
    '''Spheromak eigenmode with the mesh size quartered, checking convergence with resolution and element order'''
    errs = [3.257155111957006e-07, 2.090369020180253e-09, 8.601148342547016e-12]
    results = mp_run(run_sph_case,(0.05/4.0,order))
    assert validate_sph(results,errs[order-2])

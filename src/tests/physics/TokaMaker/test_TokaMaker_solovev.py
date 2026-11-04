#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: fixed-boundary Solov'ev equilibrium (analytic flux and X-points)
import pytest
import numpy as np
from tokamaker_testing import mp_run
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain


def run_solo_case(mesh_resolution,fe_order,mp_q):
    def solovev_psi(r_grid, z_grid,R,a,b,c0):
        # psi = np.zeros_like(r_grid)
        zeta = (np.power(r_grid,2)-np.power(R,2))/(2.0*R)
        psi_x = (a-c0)*np.power(b+c0,2)*np.power(R,4)/(8.0*np.power(c0,2))
        zeta_x = -(b+c0)*R/(2*c0)
        Z_x = np.sqrt((b+c0)*(a-c0)/(2*c0*c0))*R
        psi_grid = (b+c0)*np.power(R,2)*np.power(z_grid,2)/2.0 + c0*R*zeta*np.power(z_grid,2) + (a-c0)*np.power(R,2)*np.power(zeta,2)/2.0
        return psi_grid, psi_x, [np.sqrt(zeta_x*2*R+R*R), Z_x]
    # Fixed parameters
    R=1.0
    a=1.2
    b=-1.0
    c0=1.1
    # Build mesh
    gs_mesh = gs_Domain()
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.add_rectangle(R,0.0,0.12,0.15,'plasma')
    mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
    # Run EQ
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc)
    mygs.settings.free_boundary = False
    mygs.setup(order=fe_order,F0=1.0,full_domain=True)
    mygs.p_scale=a
    mygs.ffp_scale=b*R*R*2.0
    mygs.set_profiles(ffp_prof={'type': 'flat'},pp_prof={'type': 'flat'})
    mygs.init_psi()
    psi_solovev_TM, _, rz_x = solovev_psi(mygs.r[:,0], mygs.r[:,1],R,a,b,c0)
    mygs.set_psi(-psi_solovev_TM)
    mygs.settings.nl_tol = 1.E-14
    mygs.update_settings()
    try:
        mygs.solve()
    except ValueError:
        mp_q.put(None)
        return
    psi_TM = mygs.get_psi(False)
    # Compute error in psi
    psi_err = np.linalg.norm(psi_TM+psi_solovev_TM)
    # Compute error in X-points
    x_points, _ = mygs.get_xpoints()
    X_err = 0.0
    for i in range(2):
        diff = x_points[i,:]-rz_x
        if x_points[i,1] < 0.0:
            diff[1] = x_points[i,1]+rz_x[1]
        X_err += np.linalg.norm(diff)
    mp_q.put([psi_err, X_err])
    oftpy_dump_cov()


def validate_solo(results,psi_err_exp,X_err_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs(results[0]) > abs(psi_err_exp)*1.1:
        print("FAILED: psi error too high!")
        print("  Expected = {0}".format(psi_err_exp))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    if abs(results[1]) > abs(X_err_exp)*1.1:
        print("FAILED: X-point error too high!")
        print("  Expected = {0}".format(X_err_exp))
        print("  Actual =   {0}".format(results[1]))
        test_result = False
    return test_result


# Test runners for Solov'ev cases
@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
def test_solo_h1(order):
    errs = [
        [3.2048631614233643e-07,0.00014929412629149645],
        [8.919954733021135e-10,4.659825491095631e-07],
        [5.084454462338564e-15,4.224329766330554e-12]
    ]
    results = mp_run(run_solo_case,(0.015,order))
    assert validate_solo(results,errs[order-2][0],errs[order-2][1])
@pytest.mark.parametrize("order", (2,3,4))
def test_solo_h2(order):
    errs = [
        [7.725262474858205e-08,4.9243688140144384e-05],
        [1.1190059530634016e-10,2.919838380657025e-08],
        [1.0424769098635496e-14,1.3388421173180407e-12]
    ]
    results = mp_run(run_solo_case,(0.015/2.0,order))
    assert validate_solo(results,errs[order-2][0],errs[order-2][1])
@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
def test_solo_h3(order):
    errs = [
        [2.0607919004158514e-08,5.955338556344096e-06],
        [1.3950375633902016e-11,1.154542061756696e-09],
        [2.0552832098467707e-14,1.1263537155091759e-12]
    ]
    results = mp_run(run_solo_case,(0.015/4.0,order))
    assert validate_solo(results,errs[order-2][0],errs[order-2][1])


#============================================================================
def run_nonfinite_case(mode,mp_q):
    gs_mesh = gs_Domain()
    gs_mesh.define_region('plasma',0.015,'plasma')
    gs_mesh.add_rectangle(1.0,0.0,0.12,0.15,'plasma')
    mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc)
    mygs.settings.free_boundary = False
    mygs.settings.maxits = 200
    mygs.setup(order=2,F0=1.0,full_domain=True)
    mygs.p_scale=np.nan if mode == 'source' else 1.2
    mygs.ffp_scale=-2.0
    mygs.set_profiles(ffp_prof={'type': 'flat'},pp_prof={'type': 'flat'})
    mygs.init_psi()
    if mode == 'psi':
        psi = mygs.get_psi(False)
        psi[0] = np.nan
        mygs.set_psi(psi)
    try:
        mygs.solve()
    except ValueError as err:
        mp_q.put(str(err))
        return
    mp_q.put('converged')


# A non-finite initial flux or source must fail at once (error -10), not run to "maxits"
@pytest.mark.parametrize("mode", ('psi','source'))
def test_nonfinite_abort(mode):
    result = mp_run(run_nonfinite_case,(mode,))
    assert result is not None
    assert 'Non-finite state' in result

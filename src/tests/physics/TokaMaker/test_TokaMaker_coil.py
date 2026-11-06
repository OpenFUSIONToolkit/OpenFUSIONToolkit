#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: vacuum field and mutual inductance of rectangular coils (analytic Green's function)
import pytest
import numpy as np
from scipy.integrate import dblquad
from tokamaker_testing import mp_run, validate_errors
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.TokaMaker.util import eval_green


def run_coil_case(mesh_resolution,fe_order,dist,mp_q):
    # Plasma region (unused here), driven coil COIL1, and passive coil COIL2 used for the mutual inductance
    px1,py1,pdx,pdy = 0.4,0.4,0.2,0.2
    cx1,cy1,cdx,cdy = 0.8,0.8,0.1,0.1
    cx2,cy2 = 0.8,0.4
    def coil_green(zc,rc,r,z):
        '''Green's function from a filament at (rc,zc), weighted by the coil current distribution if set'''
        if dist is None:
            return eval_green(np.array([[r,z]]),np.array([rc,zc]))[0]
        else:
            return eval_green(np.array([[r,z]]),np.array([rc,zc]))[0]*dist(rc,zc)
    def masked_err(point_mask,gs_obj,psi,sort_ind):
        '''Analytic and computed flux at the boundary nodes selected by `point_mask`, sorted along one coordinate'''
        bdry_points = gs_obj.r[point_mask,:]
        sort_ind = bdry_points[:,sort_ind].argsort()
        psi_bdry = psi[point_mask]
        psi_bdry = psi_bdry[sort_ind]
        bdry_points = bdry_points[sort_ind]
        green = np.zeros((bdry_points.shape[0],))
        for i in range(bdry_points.shape[0]):
            green[i], _ = dblquad(coil_green,cx1-cdx/2,cx1+cdx/2,cy1-cdy/2,cy1+cdy/2,args=(bdry_points[i,0],bdry_points[i,1]))
        return green, psi_bdry
    def analytic_mutual():
        '''Mutual inductance of the two coils by integrating the Green's function over both cross-sections'''
        def mutual_integrand(z,r):
            integrand, _ = dblquad(coil_green,cx1-cdx/2,cx1+cdx/2,cy1-cdy/2,cy1+cdy/2,args=(r,z))
            return integrand
        mutual, _ = dblquad(mutual_integrand,cx2-cdx/2,cx2+cdx/2,cy2-cdy/2,cy2+cdy/2)
        return mutual*2.0*np.pi/(cdx*cdy)/(cdx*cdy)
    gs_mesh = gs_Domain(rextent=1.0,zextents=[0.0,1.0])
    gs_mesh.define_region('air',mesh_resolution,'boundary')
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    # Fixed resolution in the driven coil so convergence reflects the vacuum region only
    gs_mesh.define_region('coil1',0.01,'coil')
    gs_mesh.define_region('coil2',mesh_resolution,'coil')
    gs_mesh.add_rectangle(px1,py1,pdx,pdy,'plasma')
    gs_mesh.add_rectangle(cx1,cy1,cdx,cdy,'coil1')
    gs_mesh.add_rectangle(cx2,cy2,cdx,cdy,'coil2')
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    coil_dict = gs_mesh.get_coils()
    cond_dict = gs_mesh.get_conductors()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=fe_order)
    # Total current equal to the coil area, so the current density is unity
    mygs.set_coil_currents({'COIL1': cdx*cdy})
    if dist is not None:
        mygs.set_coil_current_dist('COIL1',dist(mygs.r[:,0],mygs.r[:,1]))
    try:
        vac_Eq = mygs.vac_solve()
        psi0 = vac_Eq.get_psi(False)
    except ValueError:
        mp_q.put(None)
        return
    # TokaMaker's psi has the opposite sign to the Green's function convention
    Lmat = mygs.get_coil_Lmat()
    Mcc = analytic_mutual()
    mutual_err = abs((Lmat[0,1]+Mcc)/Mcc)
    # Compare on the top, outer, and bottom domain boundaries
    green1, psi1 = masked_err(mygs.r[:,1]==1.0,mygs,psi0,0)
    green2, psi2 = masked_err(mygs.r[:,0]==1.0,mygs,psi0,1)
    green3, psi3 = masked_err(mygs.r[:,1]==0.0,mygs,psi0,0)
    # Drop the nodes on the axis (R = 0), where both fluxes vanish
    green_full = np.hstack((green1[1:], green2, green3[1:]))
    psi_full = np.hstack((psi1[1:], psi2, psi3[1:]))
    psi_err = np.linalg.norm(green_full+psi_full)/np.linalg.norm(green_full)
    mp_q.put([psi_err,mutual_err])
    oftpy_dump_cov()


def validate_coil(results,psi_err_exp,mutual_err_exp):
    return validate_errors(results,(psi_err_exp,mutual_err_exp),('psi','coil mutual'))


def coil_dist(r,z):
    '''Linearly varying current distribution for the non-uniform coil case'''
    return r-z


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h1(order,dist_coil):
    '''Free-boundary vacuum field of a rectangular coil with uniform or linearly varying current density on the
    coarsest mesh: the boundary flux and the coil mutual inductance must match integrals of the analytic Green's
    function to within their recorded errors.

    Reference: C. Hansen et al., Comput. Phys. Commun. 298, 109111 (2024), Sec. 4.1.3.
    '''
    if dist_coil:
        errs = np.r_[2.34993E-02, 7.68143E-03, 7.54414E-04]
        mutual_errs = np.r_[9.88100E-03, 5.16977E-03, 3.54599E-04]
        results = mp_run(run_coil_case,(0.1,order,coil_dist))
    else:
        errs = np.r_[9.62847E-03, 7.45252E-04, 4.29408E-05]
        mutual_errs = np.r_[1.01747E-03, 3.92000E-04, 5.27038E-05]
        results = mp_run(run_coil_case,(0.1,order,None))
    assert validate_coil(results,errs[order-2],mutual_errs[order-2])


@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h2(order,dist_coil):
    '''Vacuum coil field and mutual inductance with the mesh size halved, checking convergence'''
    if dist_coil:
        errs = np.r_[4.32354E-03, 2.15974E-04, 1.97917E-05]
        mutual_errs = np.r_[4.54609E-04, 1.25920E-04, 4.84935E-06]
        results = mp_run(run_coil_case,(0.1/2.0,order,coil_dist))
    else:
        errs = np.r_[3.70634E-03, 5.16587E-05, 3.55922E-06]
        mutual_errs = np.r_[1.14499E-03, 2.03449E-05, 1.49148E-06]
        results = mp_run(run_coil_case,(0.1/2.0,order,None))
    assert validate_coil(results,errs[order-2],mutual_errs[order-2])


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h3(order,dist_coil):
    '''Vacuum coil field and mutual inductance with the mesh size quartered, checking convergence'''
    if dist_coil:
        errs = np.r_[1.30686E-03, 1.00862E-05, 4.33982E-07]
        mutual_errs = np.r_[4.30803E-04, 3.04001E-06, 1.69941E-07]
        results = mp_run(run_coil_case,(0.1/4.0,order,coil_dist))
    else:
        errs = np.r_[9.12771E-04, 1.58232E-06, 3.04067E-07]
        mutual_errs = np.r_[3.31022E-04, 2.88318E-07, 2.19675E-07]
        results = mp_run(run_coil_case,(0.1/4.0,order,None))
    assert validate_coil(results,errs[order-2],mutual_errs[order-2])

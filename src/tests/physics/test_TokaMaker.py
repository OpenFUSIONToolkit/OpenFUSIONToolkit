import os
import sys
import time
import multiprocessing
import json
import pytest
import numpy as np
from scipy.special import jv, jn_zeros
from scipy.integrate import dblquad
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..','..','python')))
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.util import mu0, eC
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh, load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_isoflux, eval_green, create_power_flux_fun, create_isoflux_xpts, xpoints_from_moments


def mp_run(target,args,timeout=30):
    if os.environ.get('OFT_DEBUG_TEST', 0):
        timeout *= 4
    os.chdir(test_dir)
    mp_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=target, args=args + (mp_q,))
    p.start()
    start = time.time()
    while time.time() - start <= timeout:
        if not p.is_alive():
            break
        time.sleep(.5)
    else: # Reached timeout
        print("Timeout reached")
        p.terminate()
        p.join()
        return None
    # Completed successfully
    try:
        test_result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return None
    p.join()
    return test_result


def validate_dict(results,dict_exp,tol_dict=None):
    default_tol = {
        'delta': 5.E-2,
        'deltaU': 5.E-2,
        'deltaL': 5.E-2,
        # R at a Z extremum is taken from the nearest traced point, so it carries
        # O(dtheta) sampling error and moves with element order and mesh resolution
        'fsa_R_at_Zmin': 2.E-2,
        'fsa_R_at_Zmax': 2.E-2
    }
    if tol_dict is None:
        tol_dict = default_tol
    else:
        # Merge user-provided tolerances on top of the defaults
        merged = dict(default_tol)
        merged.update(tol_dict)
        tol_dict = merged
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    for key, exp_val in dict_exp.items():
        result_val = results[0].get(key,None)
        if result_val is None:
            print('FAILED: key "{0}" not present!'.format(key))
            result_val = False
        else:
            if type(exp_val) is list:
                for i in range(len(exp_val)):
                    if exp_val[i] is None:
                        continue
                    if abs((result_val[i]-exp_val[i])/exp_val[i]) > tol_dict.get(key,1.E-2):
                        print("FAILED: {0} ({1}) error too high!".format(key,i))
                        print("  Expected = {0:E}".format(exp_val[i]))
                        print("  Actual =   {0:E}".format(result_val[i]))
                        test_result = False
            else:
                if abs((result_val-exp_val)/exp_val) > tol_dict.get(key,1.E-2):
                    print("FAILED: {0} error too high!".format(key))
                    print("  Expected = {0:E}".format(exp_val))
                    print("  Actual =   {0:E}".format(result_val))
                    test_result = False
    return test_result


def validate_eqdsk(file_test,file_ref,helicity=1.0):
    from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk
    try:
        test_data = read_eqdsk(file_test)
    except:
        print("FAILED: Could not read result EQDSK")
        return False
    try:
        ref_data = read_eqdsk(file_ref)
    except:
        print("FAILED: Could not read reference EQDSK")
        return False
    ref_data['bcentr'] *= helicity
    ref_data['fpol'] *= helicity
    ref_data['qpsi'] *= helicity
    test_result = True
    for key, exp_val in ref_data.items():
        result_val = test_data.get(key,None)
        if result_val is None:
            print('FAILED: key "{0}" not present in result!'.format(key))
            result_val = False
        else:
            if key == 'case':
                continue
            if isinstance(exp_val,np.ndarray):
                if np.linalg.norm(exp_val-result_val)/np.linalg.norm(exp_val) > 1.E-2:
                    print("FAILED: {0} error too high!".format(key))
                    print("  Actual =   {0}".format(np.linalg.norm(exp_val-result_val)/np.linalg.norm(exp_val)))
                    test_result = False
            else:
                if abs((result_val-exp_val)/exp_val) > 1.E-2:
                    print("FAILED: {0} error too high!".format(key))
                    print("  Expected = {0}".format(exp_val))
                    print("  Actual =   {0}".format(result_val))
                    test_result = False
    return test_result


def validate_ifile(ifile_test,ifile_ref,helicity=1.0):
    from OpenFUSIONToolkit.TokaMaker.util import read_ifile
    try:
        test_data = read_ifile(ifile_test)
    except:
        print("FAILED: Could not read result i-file")
        return False
    try:
        ref_data = read_ifile(ifile_ref)
    except:
        print("FAILED: Could not read reference i-file")
        return False
    ref_data['f'] *= helicity
    ref_data['q'] *= helicity
    test_result = True
    for key, exp_val in ref_data.items():
        result_val = test_data.get(key,None)
        if result_val is None:
            print('FAILED: key "{0}" not present in result!'.format(key))
            result_val = False
        else:
            if isinstance(exp_val,np.ndarray):
                if np.linalg.norm(exp_val-result_val)/np.linalg.norm(exp_val) > 1.E-2:
                    print("FAILED: {0} error too high!".format(key))
                    print("  Actual =   {0}".format(np.linalg.norm(exp_val-result_val)/np.linalg.norm(exp_val)))
                    test_result = False
            else:
                if result_val != exp_val:
                    print("FAILED: {0} error too high!".format(key))
                    print("  Expected = {0}".format(exp_val))
                    print("  Actual =   {0}".format(result_val))
                    test_result = False
    return test_result


#============================================================================
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
def run_sph_case(mesh_resolution,fe_order,mp_q):
    def spheromak_psi(r_grid,z_grid,a,h):
        gamma_11 = jn_zeros(1,1)[0]*r_grid/a
        x_01 = jn_zeros(0,1)[0]
        norm = x_01*jv(1,x_01)
        return gamma_11*jv(1,gamma_11)*np.sin(np.pi*z_grid/h)/norm
    # Build mesh
    gs_mesh = gs_Domain()
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.add_rectangle(0.5,0.5,1.0,1.0,'plasma')
    mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
    # Run EQ
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
    psi_TM /= psi_TM.dot(psi_eig_TM)/psi_eig_TM.dot(psi_eig_TM)
    # Compute error in psi
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


# Test runners for Spheromak cases
@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h1(order):
    errs = [2.039674417912789e-05, 5.103597862537552e-07, 8.088772274705608e-09]
    results = mp_run(run_sph_case,(0.05,order))
    assert validate_sph(results,errs[order-2])
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h2(order):
    errs = [2.5203856661960034e-06, 3.279268054674832e-08, 2.5185712724779513e-10]
    results = mp_run(run_sph_case,(0.05/2.0,order))
    assert validate_sph(results,errs[order-2])
@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
def test_spheromak_h3(order):
    errs = [3.257155111957006e-07, 2.090369020180253e-09, 8.601148342547016e-12]
    results = mp_run(run_sph_case,(0.05/4.0,order))
    assert validate_sph(results,errs[order-2])


#============================================================================
def run_coil_case(mesh_resolution,fe_order,dist,mp_q):
    px1,py1,pdx,pdy = 0.4,0.4,0.2,0.2
    cx1,cy1,cdx,cdy = 0.8,0.8,0.1,0.1
    cx2,cy2 = 0.8,0.4
    def coil_green(zc,rc,r,z):
        if dist is None:
            return eval_green(np.array([[r,z]]),np.array([rc,zc]))[0]
        else:
            return eval_green(np.array([[r,z]]),np.array([rc,zc]))[0]*dist(rc,zc)
    def masked_err(point_mask,gs_obj,psi,sort_ind):
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
        def mutual_integrand(z,r):
            integrand, _ = dblquad(coil_green,cx1-cdx/2,cx1+cdx/2,cy1-cdy/2,cy1+cdy/2,args=(r,z))
            return integrand
        mutual, _ = dblquad(mutual_integrand,cx2-cdx/2,cx2+cdx/2,cy2-cdy/2,cy2+cdy/2)
        return mutual*2.0*np.pi/(cdx*cdy)/(cdx*cdy)
    # Build mesh
    gs_mesh = gs_Domain(rextent=1.0,zextents=[0.0,1.0])
    gs_mesh.define_region('air',mesh_resolution,'boundary')
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.define_region('coil1',0.01,'coil')
    gs_mesh.define_region('coil2',mesh_resolution,'coil')
    gs_mesh.add_rectangle(px1,py1,pdx,pdy,'plasma')
    gs_mesh.add_rectangle(cx1,cy1,cdx,cdy,'coil1')
    gs_mesh.add_rectangle(cx2,cy2,cdx,cdy,'coil2')
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    coil_dict = gs_mesh.get_coils()
    cond_dict = gs_mesh.get_conductors()
    # Run EQ
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=fe_order)
    mygs.set_coil_currents({'COIL1': cdx*cdy})
    if dist is not None:
        mygs.set_coil_current_dist('COIL1',dist(mygs.r[:,0],mygs.r[:,1]))
    try:
        vac_Eq = mygs.vac_solve()
        psi0 = vac_Eq.get_psi(False)
    except ValueError:
        mp_q.put(None)
        return
    # Get coil mutual matrix
    Lmat = mygs.get_coil_Lmat()
    Mcc = analytic_mutual()
    mutual_err = abs((Lmat[0,1]+Mcc)/Mcc)

    # Get analytic result
    green1, psi1 = masked_err(mygs.r[:,1]==1.0,mygs,psi0,0)
    green2, psi2 = masked_err(mygs.r[:,0]==1.0,mygs,psi0,1)
    green3, psi3 = masked_err(mygs.r[:,1]==0.0,mygs,psi0,0)
    # Compute error in psi
    green_full = np.hstack((green1[1:], green2, green3[1:]))
    psi_full = np.hstack((psi1[1:], psi2, psi3[1:]))
    psi_err = np.linalg.norm(green_full+psi_full)/np.linalg.norm(green_full)
    mp_q.put([psi_err,mutual_err])
    oftpy_dump_cov()


def validate_coil(results,psi_err_exp,mutual_err_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs(results[0]) > 1.1*abs(psi_err_exp):
        print("FAILED: psi error too high!")
        print("  Expected = {0:.5E}".format(psi_err_exp))
        print("  Actual =   {0:.5E}".format(results[0]))
        test_result = False
    if abs(results[1]) > 1.1*abs(mutual_err_exp):
        print("FAILED: coil mutual error too high!")
        print("  Expected = {0:.5E}".format(mutual_err_exp))
        print("  Actual =   {0:.5E}".format(results[1]))
        test_result = False
    return test_result


# Test runners for vacuum coil cases
def coil_dist(r,z):
    return r-z

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h1(order,dist_coil):
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
    if dist_coil:
        errs = np.r_[1.30686E-03, 1.00862E-05, 4.33982E-07]
        mutual_errs = np.r_[4.30803E-04, 3.04001E-06, 1.69941E-07]
        results = mp_run(run_coil_case,(0.1/4.0,order,coil_dist))
    else:
        errs = np.r_[9.12771E-04, 1.58232E-06, 3.04067E-07]
        mutual_errs = np.r_[3.31022E-04, 2.88318E-07, 2.19675E-07]
        results = mp_run(run_coil_case,(0.1/4.0,order,None))
    assert validate_coil(results,errs[order-2],mutual_errs[order-2])


#============================================================================
def run_ITER_case(mesh_resolution,fe_orders,test_type,helicity,mp_q):
    def create_mesh():
        with open('ITER_geom.json','r') as fid:
            ITER_geom = json.load(fid)
        plasma_dx = 0.15/mesh_resolution
        coil_dx = 0.2/mesh_resolution
        vv_dx = 0.3/mesh_resolution
        vac_dx = 0.6/mesh_resolution
        gs_mesh = gs_Domain()
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('vacuum1',vv_dx,'vacuum')
        gs_mesh.define_region('vacuum2',vv_dx,'vacuum')
        gs_mesh.define_region('vv1',vv_dx,'conductor',eta=6.9E-7)
        gs_mesh.define_region('vv2',vv_dx,'conductor',eta=6.9E-7)
        for key, coil in ITER_geom['coils'].items():
            if not key.startswith('VS'):
                gs_mesh.define_region(key,coil_dx,'coil')
        gs_mesh.define_region('VSU',coil_dx,'coil',coil_set='VS',nTurns=1.0)
        gs_mesh.define_region('VSL',coil_dx,'coil',coil_set='VS',nTurns=-1.0)
        gs_mesh.add_polygon(ITER_geom['limiter'],'plasma',parent_name='vacuum1')             # Define the shape of the limiter
        gs_mesh.add_annulus(ITER_geom['inner_vv'][0],'vacuum1',ITER_geom['inner_vv'][1],'vv1',parent_name='vacuum2') # Define the shape of the VV
        gs_mesh.add_annulus(ITER_geom['outer_vv'][0],'vacuum2',ITER_geom['outer_vv'][1],'vv2',parent_name='air') # Define the shape of the VV
        for key, coil in ITER_geom['coils'].items():
            if key.startswith('VS'):
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='vacuum1')
            else:
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='air')
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        coil_dict = gs_mesh.get_coils()
        cond_dict = gs_mesh.get_conductors()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict,'ITER_mesh.h5')
    if not os.path.exists('ITER_mesh.h5'):
        try:
            create_mesh()
        except Exception as e:
            print(e)
            mp_q.put(None)
            return
    # Run EQ
    mygs = None
    myOFT = OFT_env(nthreads=-1)
    for fe_order in fe_orders:
        mygs_last = mygs
        mygs = TokaMaker(myOFT)
        mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('ITER_mesh.h5')
        mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
        mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
        mygs.setup(order=fe_order,F0=helicity*5.3*6.2)
        #
        if test_type.startswith('eig'):
            if test_type == 'eig_dep':
                eig_vals, _ = mygs.eig_wall(10)
                mp_q.put([{'Tau_w': 1.0/eig_vals[:5,0]}])
            else:
                eig_vals, _ = mygs.compute_wall_modes(10)
                mp_q.put([{'Tau_w': eig_vals[:5]}])
            oftpy_dump_cov()
            return
        #
        mygs.set_coil_vsc({'VS': 1.0})
        #
        coil_bounds = {key: [-50.E6, 50.E6] for key in mygs.coil_sets}
        mygs.set_coil_bounds(coil_bounds)
        #
        Ip_target=15.6E6
        P0_target=6.2E5
        mygs.set_targets(Ip=Ip_target, pax=P0_target)
        isoflux_pts = np.array([
            [ 8.20,  0.41],
            [ 8.06,  1.46],
            [ 7.51,  2.62],
            [ 6.14,  3.78],
            [ 4.51,  3.02],
            [ 4.26,  1.33],
            [ 4.28,  0.08],
            [ 4.49, -1.34],
            [ 7.28, -1.89],
            [ 8.00, -0.68]
        ])
        x_point = np.array([[5.125, -3.4],])
        mygs.set_isoflux(np.vstack((isoflux_pts,x_point)))
        mygs.set_saddles(x_point)
        # Set regularization weights
        regularization_terms = []
        for name in mygs.coil_sets:
            if name.startswith('CS'):
                if name.startswith('CS1'):
                    regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=2.E-2))
                else:
                    regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-2))
            elif name.startswith('PF'):
                regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-2))
            elif name.startswith('VS'):
                regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-2))
        regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E2))
        mygs.set_coil_reg(reg_terms=regularization_terms)
        #
        ffp_prof = create_power_flux_fun(40,1.5,2.0)
        pp_prof = create_power_flux_fun(40,4.0,1.0)
        mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
        #
        R0 = 6.3
        Z0 = 0.5
        a = 2.0
        kappa = 1.4
        delta = 0.0
        try:
            mygs.init_psi(R0, Z0, a, kappa, delta)
            EQ_obj, nl_its = mygs.solve(return_its=True)
        except ValueError:
            mp_q.put(None)
            return
        if test_type.startswith('stab'):
            if test_type == 'stab_dep':
                eig_vals, eig_modes = mygs.eig_td(-1.E2,10,False)
                growth_rates = -eig_vals[:,0]
            else:
                growth_rates, eig_modes = mygs.compute_linear_stability(1.E2,10,False)
            # Run brief nonlinear evolution
            psi0 = mygs.get_psi(False)
            eig_sign = eig_modes[0,(mygs.r[:,1]-R0)>0.0][abs(eig_modes[0,(mygs.r[:,1]-R0)>0.0]).argmax()]
            psi_ic = psi0-0.01*eig_modes[0,:]*(mygs.psi_bounds[1]-mygs.psi_bounds[0])/eig_sign
            mygs.set_psi(psi_ic)
            mygs.set_saddles(None)
            mygs.set_isoflux(None)
            dt = 0.1/abs(growth_rates[0])
            mygs.setup_td(dt,1.E-13,1.E-11)
            sim_time = 0.0
            for i in range(5):
                sim_time, _, _, _, _ = mygs.step_td(sim_time,dt)
            psi1 = mygs.get_psi(False)
            mp_q.put([{'gamma': growth_rates[:5], 'nl_change': np.linalg.norm(psi1-psi0)}])
            oftpy_dump_cov()
            return
        mygs.save_eqdsk('test.eqdsk',lcfs_pressure=6.E4)
        if test_type == 'io':
            EQ_obj.save_TokaMaker('test_eq.h5')
            mygs.init_psi(R0, Z0, a, kappa, delta)
            mygs.replace_eq(source_file='test_eq.h5')
            EQ_obj = mygs.copy_eq()
            _, nl_its = mygs.solve(return_its=True)
        eq_info = EQ_obj.get_stats(li_normalization='ITER')
        Lmat = mygs.get_coil_Lmat()
        eq_info['LCS1'] = Lmat[mygs.coil_sets['CS1U']['id'],mygs.coil_sets['CS1U']['id']]
        eq_info['MCS1_plasma'] = Lmat[mygs.coil_sets['CS1U']['id'],-1]
        eq_info['Lplasma'] = Lmat[-1,-1]
        eq_info['nl_its'] = nl_its
        # Flux surface averages and per-surface shape parameters (see `get_fsa`)
        fsa = EQ_obj.get_fsa(psi=np.r_[0.25, 0.5, 0.9])
        for key in ('q', 'F', 'dV/dPsi', '<|grad psi|>', '<|grad psi|^2>', '<Bp^2>',
                    '<1/B^2>', 'R_min', 'R_max', 'Z_min', 'Z_max', 'R_at_Zmin', 'R_at_Zmax'):
            eq_info['fsa_' + key] = fsa[key].tolist()
    if test_type.startswith('recon'):
        import random
        from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
        # Sample constraint locations
        with open('ITER_geom.json','r') as fid:
            ITER_geom = json.load(fid)
            probe_contour = np.asarray(ITER_geom['inner_vv'][0])
        dl_contour = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(probe_contour,axis=0),axis=1))]
        probe_contour_new = np.zeros((100,2))
        probe_contour_new[:,0] = np.interp(np.linspace(0.0,dl_contour[-1],100),dl_contour,probe_contour[:,0])
        probe_contour_new[:,1] = np.interp(np.linspace(0.0,dl_contour[-1],100),dl_contour,probe_contour[:,1])
        B_locs = []
        for i, pt in enumerate(probe_contour_new):
            if i % 5 == 0:
                B_locs.append(pt)
        B_locs = np.asarray(B_locs)
        # Setup constraints
        random.seed(42)
        myrecon = reconstruction(mygs)
        noise_amp = random.gauss(0.0,1.0)
        Ip_noised = eq_info['Ip']*(1.0+noise_amp*0.05)
        myrecon.set_Ip(Ip_noised, err=0.05*eq_info['Ip'])
        flux_vals = []
        field_eval = mygs.get_field_eval('PSI')
        for i in range(B_locs.shape[0]):
            B_tmp = field_eval.eval(B_locs[i,:])
            noise_amp = random.gauss(0.0,1.0)
            flux_vals.append(B_tmp[0])
            psi_val = B_tmp[0]*2.0*np.pi
            myrecon.add_flux_loop(B_locs[i,:], psi_val*(1.0 + noise_amp*0.05), err=abs(psi_val*0.05))
        field_eval = mygs.get_field_eval('B')
        for i in range(B_locs.shape[0]):
            B_tmp = field_eval.eval(B_locs[i,:])
            noise_amp = random.gauss(0.0,1.0)
            myrecon.add_Mirnov(B_locs[i,:], np.r_[1.0,0.0,0.0], B_tmp[0] + noise_amp*abs(B_tmp[0]*0.05), err=abs(B_tmp[0]*0.05))
            noise_amp = random.gauss(0.0,1.0)
            myrecon.add_Mirnov(B_locs[i,:], np.r_[0.0,0.0,1.0], B_tmp[2] + noise_amp*abs(B_tmp[2]*0.05), err=abs(B_tmp[2]*0.05))
        coil_currents, _ = mygs.get_coil_currents()
        for key in coil_currents:
            noise_amp = random.gauss(0.0,1.0)
            coil_currents[key] *= 1.0+noise_amp*0.05
        # Compute starting equilibrium
        mygs.set_isoflux(None)
        mygs.set_saddles(None)
        mygs.set_targets(Ip=Ip_noised,Ip_ratio=2.0)
        if test_type == 'recon':
            mygs.settings.ffp_target_weight=1.0/((0.05*abs(Ip_noised))*mu0)
            mygs.settings.pp_target_weight=1.0
            mygs.update_settings()
        mygs.set_flux(B_locs,np.array(flux_vals))
        regularization_terms = []
        for name in coil_currents:
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=coil_currents[name],weight=1.E-1))
        regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E2))
        mygs.set_coil_reg(reg_terms=regularization_terms)
        coil_err = {name: 0.05*abs(current) for name, current in coil_currents.items()}
        myrecon.set_coil_currents(coil_currents,coil_err)
        R0 = 6.3
        Z0 = 0.5
        a = 1.0
        kappa = 1.0
        delta = 0.0
        mygs.init_psi(R0, Z0, a, kappa, delta)
        if test_type == 'recon':
            mygs.settings.maxits=300
        elif test_type == 'recon_legacy':
            mygs.settings.maxits=100
        mygs.update_settings()
        mygs.solve()
        # Perform reconstruction
        if test_type == 'recon':
            mygs.settings.pm = False
            mygs.update_settings()
            mygs.set_isoflux(None)
            mygs.set_saddles(None)
            mygs.set_targets(Ip=Ip_noised)
            mygs.settings.ffp_target_weight=1.0/((0.05*abs(Ip_noised))*mu0)
            myrecon.setup_constraints()
            # Configure reconstruction settings
            myrecon.settings.fitF = False
            myrecon.settings.fitP = False
            myrecon.settings.fit_Pscale = True
            myrecon.settings.fit_FFPscale = False
            myrecon.settings.fitR0 = False
            myrecon.settings.fitZ0 = False
            myrecon.settings.fixedCentering = False
            myrecon.settings.dx = 1.E-2
            x0_nl, _ = myrecon.setup_get_opt()
            from scipy.optimize import least_squares
            opt_result = least_squares(myrecon.opt_error, x0_nl, jac=myrecon.opt_error_jacobian, method='lm', ftol=1.E-3, args=(myrecon,False))
            mygs.settings.pm = True
            mygs.update_settings()
        elif test_type == 'recon_legacy':
            mygs.set_isoflux(None)
            mygs.set_flux(None,None)
            mygs.set_saddles(None)
            mygs.set_targets(R0=mygs.o_point[0],V0=mygs.o_point[1])
            myrecon.settings.fit_Pscale = False
            myrecon.settings.fitR0 = True
            myrecon.settings.fitCoils = True
            myrecon.settings.pm = False
            _ = myrecon.reconstruct()
        #
        eq_info = mygs.get_stats(li_normalization='ITER')
        eq_info['LCS1'] = Lmat[mygs.coil_sets['CS1U']['id'],mygs.coil_sets['CS1U']['id']]
        eq_info['MCS1_plasma'] = Lmat[mygs.coil_sets['CS1U']['id'],-1]
        eq_info['Lplasma'] = Lmat[-1,-1]
    # Test deletion if multiple cases
    if mygs_last is not None:
        del mygs_last
    # Save equilibrium to gEQDSK and i-file format
    mygs.save_eqdsk('tokamaker.eqdsk',nr=64,nz=64,lcfs_pad=0.001)
    mygs.save_ifile('tokamaker.ifile',npsi=64,ntheta=64,lcfs_pad=0.001)
    # Save final one
    mp_q.put([eq_info])
    oftpy_dump_cov()


# Test runners for ITER test cases
@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_eig(order):
    exp_dict = {
        'Tau_w': [6.619977E-01, 3.479492E-01, 2.554444E-01, 1.910381E-01, 1.782464E-01]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),'eig',1.0))
    assert validate_dict(results,exp_dict)
    # Test deprecated interface
    results = mp_run(run_ITER_case,(1.0,(order,),'eig_dep',1.0))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_stability(order):
    exp_dict = {
        'gamma': [12.3620, -1.83981, -3.41613, -5.12470, -6.53393],
        'nl_change': [225.4421413167051, 338.0113029638385][order-2]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),'stab',1.0))
    assert validate_dict(results,exp_dict)
    # Test deprecated interface
    results = mp_run(run_ITER_case,(1.0,(order,),'stab_dep',1.0))
    assert validate_dict(results,exp_dict)

ITER_eq_dict = {
    'Ip': 15599996.692463942,
    'Ip_centroid': [6.20273409, 0.52959503],
    'kappa': 1.8728151512244395,
    'kappaU': 1.7634853298971116,
    'kappaL': 1.9821449725517677,
    'delta': 0.4721203463868737,
    # 'deltaU': 0.40521771808760293,
    'deltaL': 0.5390229746861446,
    'R_geo': 6.222328618622752,
    'a_geo': 1.9835670211775267,
    'vol': 820.212921617247,
    'q_0': 0.8232444101221106,
    'q_95': 2.7602989886308738,
    'P_ax': 619225.017325726,
    'W_MHD': 242986393.97329777,
    'beta_pol': 42.427927348488936,
    'dflux': 1.540293464599462,
    'tflux': 121.86081608235014,
    'l_i': 0.9054096856166233,
    'beta_tor': 1.7798144109869558,
    'beta_n': 1.1951205307278518,
    'LCS1': 2.4858609418809336e-06,
    'MCS1_plasma': 8.931779419000401e-07,
    'Lplasma': 1.1900576990802187e-05,
    # Flux surface averages and shape at psi_N = [0.25, 0.5, 0.9] (see `get_fsa`)
    'fsa_q': [0.8957740190095893, 1.0925386896870009, 2.244697561979323],
    'fsa_F': [33.974461335230956, 33.24150425404669, 32.863767384960504],
    'fsa_dV/dPsi': [-40.809941802973256, -48.739421299593786, -88.3635676135186],
    'fsa_<|grad psi|>': [6.631303992418605, 8.12645638205315, 6.722632628048097],
    'fsa_<|grad psi|^2>': [44.962611735546105, 68.65358375003274, 53.596719328235814],
    'fsa_<Bp^2>': [1.1232107307935824, 1.7334030132908238, 1.3698497417461675],
    'fsa_<1/B^2>': [0.03394218780262852, 0.03459374214812227, 0.03358354423001478],
    'fsa_R_min': [5.462226095551001, 5.036393890349691, 4.409277529477673],
    'fsa_R_max': [7.220524711028818, 7.593373222623317, 8.085943200676986],
    'fsa_Z_min': [-0.7383764856169046, -1.3583537534271921, -2.508416991476664],
    'fsa_Z_max': [1.8005752842848168, 2.414360593511412, 3.5092131308630288],
    'fsa_R_at_Zmin': [6.304387105351003, 6.204424153705717, 5.796505368206957],
    'fsa_R_at_Zmax': [6.272540830880864, 6.157453542216596, 5.750257254804224],
}

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
@pytest.mark.parametrize("helicity", (1.0,-1.0))
def test_ITER_eq(order,helicity):
    eq_dict = ITER_eq_dict.copy()
    eq_dict['tflux'] *= helicity
    eq_dict['dflux'] *= helicity
    eq_dict['q_0'] *= helicity
    eq_dict['q_95'] *= helicity
    # `q` and `F` follow the sign of F0; the remaining `fsa_*` entries are invariant.
    # Rebind rather than scale in place: `ITER_eq_dict.copy()` is shallow, so mutating
    # these lists would leak into the other parametrized runs.
    eq_dict['fsa_q'] = [val*helicity for val in eq_dict['fsa_q']]
    eq_dict['fsa_F'] = [val*helicity for val in eq_dict['fsa_F']]
    results = mp_run(run_ITER_case,(1.0,(order,),'',helicity))
    assert validate_dict(results,eq_dict)
    assert validate_eqdsk('tokamaker.eqdsk','ITER_test.eqdsk',helicity)
    assert validate_ifile('tokamaker.ifile','ITER_test.ifile',helicity)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
@pytest.mark.parametrize("helicity", (1.0,-1.0))
def test_ITER_eq_io(order,helicity):
    eq_dict = ITER_eq_dict.copy()
    eq_dict['nl_its'] = 1
    eq_dict['tflux'] *= helicity
    eq_dict['dflux'] *= helicity
    eq_dict['q_0'] *= helicity
    eq_dict['q_95'] *= helicity
    # `q` and `F` follow the sign of F0; the remaining `fsa_*` entries are invariant.
    # Rebind rather than scale in place: `ITER_eq_dict.copy()` is shallow, so mutating
    # these lists would leak into the other parametrized runs.
    eq_dict['fsa_q'] = [val*helicity for val in eq_dict['fsa_q']]
    eq_dict['fsa_F'] = [val*helicity for val in eq_dict['fsa_F']]
    results = mp_run(run_ITER_case,(1.0,(order,),'io',helicity))
    assert validate_dict(results,eq_dict)

@pytest.mark.coverage
def test_ITER_recon():
    ITER_recon_dict = ITER_eq_dict.copy()
    ITER_recon_dict['q_0'] = 0.8417344
    ITER_recon_dict['P_ax'] = 6.575446E5
    ITER_recon_dict['W_MHD'] = 2.591990E8
    ITER_recon_dict['beta_pol'] = 46.02736
    ITER_recon_dict['dflux'] = 1.440565
    ITER_recon_dict['l_i'] = 0.8872565
    ITER_recon_dict['beta_tor'] = 1.906180
    ITER_recon_dict['beta_n'] = 1.284312
    results = mp_run(run_ITER_case,(1.0,(2,),'recon',1.0))
    assert validate_dict(results,ITER_recon_dict)

@pytest.mark.coverage
def test_ITER_recon_legacy():
    ITER_recon_dict = ITER_eq_dict.copy()
    ITER_recon_dict['vol'] = 8.119402E2
    ITER_recon_dict['q_95'] = 2.718779
    ITER_recon_dict['P_ax'] = 6.656184E5
    ITER_recon_dict['W_MHD'] = 2.607023E8
    ITER_recon_dict['beta_pol'] = 45.63260
    ITER_recon_dict['dflux'] = 1.462672
    ITER_recon_dict['tflux'] = 1.201402E2
    ITER_recon_dict['l_i'] = 0.8906732
    ITER_recon_dict['beta_tor'] = 1.934024
    ITER_recon_dict['beta_n'] = 1.294385
    results = mp_run(run_ITER_case,(1.0,(2,),'recon_legacy',1.0))
    assert validate_dict(results,ITER_recon_dict)

def test_ITER_concurrent():
    results = mp_run(run_ITER_case,(1.0,(2,3),'',1.0))
    assert validate_dict(results,ITER_eq_dict)

#============================================================================
def run_LTX_case(fe_order,test_type,mp_q):
    def create_mesh():
        with open('LTX_geom.json','r') as fid:
            LTX_geom = json.load(fid)
        plasma_dx = 0.02
        coil_dx = 0.02
        vv_dx = 0.015
        vac_dx = 0.05
        gs_mesh = gs_Domain()
        #
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('shellU',vv_dx,'conductor',eta=4.E-7,noncontinuous=True)
        gs_mesh.define_region('shellL',vv_dx,'conductor',eta=4.E-7,noncontinuous=True)
        for i, vv_segment in enumerate(LTX_geom['vv']):
            gs_mesh.define_region('vv{0}'.format(i),vv_dx,'conductor',eta=vv_segment[1])
        for key, coil in LTX_geom['coils'].items():
            if key.startswith('OH'):
                gs_mesh.define_region(key,coil_dx,'coil',nTurns=coil['nturns'],coil_set='OH')
            else:
                gs_mesh.define_region(key,coil_dx,'coil',nTurns=coil['nturns'])
        #
        gs_mesh.add_polygon(LTX_geom['limiter'],'plasma',parent_name='air')
        gs_mesh.add_polygon(LTX_geom['shell'],'shellU',parent_name='air')
        shell_lower = np.array(LTX_geom['shell'].copy()); shell_lower[:,1] *= -1.0
        gs_mesh.add_polygon(shell_lower,'shellL',parent_name='air')
        for i, vv_segment in enumerate(LTX_geom['vv']):
            gs_mesh.add_polygon(vv_segment[0],'vv{0}'.format(i),parent_name='air')
        for key, coil in LTX_geom['coils'].items():
            gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='air')
        #
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        coil_dict = gs_mesh.get_coils()
        cond_dict = gs_mesh.get_conductors()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict,'LTX_mesh.h5')
    if not os.path.exists('LTX_mesh.h5'):
        try:
            create_mesh()
        except Exception as e:
            print(e)
            mp_q.put(None)
            return
    # Run EQ
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('LTX_mesh.h5')
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=fe_order,F0=0.10752)
    #
    if test_type == 'eig':
        eig_vals, _ = mygs.compute_wall_modes(10)
        mp_q.put([{'Tau_w': eig_vals[:5]}])
        oftpy_dump_cov()
        return
    #
    mygs.set_coil_vsc({'INTERNALU': 1.0, 'INTERNALL': -1.0})
    #
    Ip_target = 8.0E4
    mygs.set_targets(Ip=Ip_target,Ip_ratio=2.0)
    isoflux_pts = create_isoflux(20,0.40,0.0,0.22,1.5,0.1)
    mygs.set_isoflux(isoflux_pts)
    # Set regularization weights
    disable_list = ('YELLOW',)
    regularization_terms = []
    for name in mygs.coil_sets:
        if name[:-1] in disable_list:
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E4))
            continue
        if name == 'OH': # OH coil has no mirror
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-1))
            continue
        elif name[-1] == 'L':
            continue
        regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-1))
        regularization_terms.append(mygs.coil_reg_term({name: 1.0, name[:-1]+'L': -1.0},target=0.0,weight=1.E2))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-4))
    mygs.set_coil_reg(reg_terms=regularization_terms)
    #
    ffp_prof = create_power_flux_fun(50,1.5,2.0)
    pp_prof = create_power_flux_fun(50,4.0,1.0)
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
    #
    mygs.init_psi(0.42,0.0,0.15,1.5,0.6)
    mygs.settings.pm=True
    mygs.update_settings()
    mygs.solve()
    if test_type == 'stab':
        eig_vals, _ = mygs.compute_linear_stability(1.E3,10,False)
        mp_q.put([{'gamma': eig_vals[:5]}])
        oftpy_dump_cov()
        return
    #
    psi_last = mygs.get_psi(False)
    mygs.set_psi_dt(psi_last,5.E-3)
    Ip_target = 9.0E4
    mygs.set_targets(Ip=Ip_target,Ip_ratio=2.0)
    mygs.solve()
    mygs.save_eqdsk('test.eqdsk')
    #
    mp_q.put([mygs.get_stats()])
    oftpy_dump_cov()

# Test runners for LTX test cases
@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_eig(order):
    exp_dict = {
        'Tau_w': [5.152566E-03, 3.953030E-03, 2.536384E-03, 2.172948E-03, 1.853882E-03]
    }
    results = mp_run(run_LTX_case,(order,'eig'))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_stability(order):
    exp_dict = {
        'gamma': [234.1051, -214.4196, -282.0877, -388.7592, -388.7592]
    }
    results = mp_run(run_LTX_case,(order,'stab'))
    assert validate_dict(results,exp_dict)

LTX_eq_dict = {
    'Ip': 90002.51679781199,
    'Ip_centroid': [ 4.05458767e-01, None],
    'kappa': 1.525595596236063,
    'kappaU': 1.5256161060199729,
    'kappaL': 1.5255750864521527,
    # 'delta': 0.1274386709874723,
    # 'deltaU': 0.13292909306765727,
    # 'deltaL': 0.1219482489072871,
    'R_geo': 0.39198831687969443,
    'a_geo': 0.2378387910900877,
    'vol': 0.6507554668762836,
    'q_0': 1.3280540982807334,
    'q_95': 5.901997513881755,
    'P_ax': 1720.958666106632,
    'W_MHD': 563.2852958452944,
    'beta_pol': 41.4052788464515,
    'dflux': 0.0009601908294886685,
    'tflux': 0.08551496495989133,
    'l_i': 1.0271521431711803,
    'beta_tor': 1.9276444168027145,
    'beta_n': 1.3972402635015146
}

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_eq(order):
    results = mp_run(run_LTX_case,(order,''))
    assert validate_dict(results,LTX_eq_dict)

#============================================================================
# Bootstrap current test (ITER-based)
#============================================================================
def run_ITER_bootstrap_case(mesh_resolution, fe_order, mp_q):
    from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap, Hmode_profiles

    # --- Mesh creation (identical to run_ITER_case) ---
    def create_mesh():
        with open('ITER_geom.json','r') as fid:
            ITER_geom = json.load(fid)
        plasma_dx = 0.15/mesh_resolution
        coil_dx = 0.2/mesh_resolution
        vv_dx = 0.3/mesh_resolution
        vac_dx = 0.6/mesh_resolution
        gs_mesh = gs_Domain()
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('vacuum1',vv_dx,'vacuum')
        gs_mesh.define_region('vacuum2',vv_dx,'vacuum')
        gs_mesh.define_region('vv1',vv_dx,'conductor',eta=6.9E-7)
        gs_mesh.define_region('vv2',vv_dx,'conductor',eta=6.9E-7)
        for key, coil in ITER_geom['coils'].items():
            if not key.startswith('VS'):
                gs_mesh.define_region(key,coil_dx,'coil')
        gs_mesh.define_region('VSU',coil_dx,'coil',coil_set='VS',nTurns=1.0)
        gs_mesh.define_region('VSL',coil_dx,'coil',coil_set='VS',nTurns=-1.0)
        gs_mesh.add_polygon(ITER_geom['limiter'],'plasma',parent_name='vacuum1')
        gs_mesh.add_annulus(ITER_geom['inner_vv'][0],'vacuum1',ITER_geom['inner_vv'][1],'vv1',parent_name='vacuum2')
        gs_mesh.add_annulus(ITER_geom['outer_vv'][0],'vacuum2',ITER_geom['outer_vv'][1],'vv2',parent_name='air')
        for key, coil in ITER_geom['coils'].items():
            if key.startswith('VS'):
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='vacuum1')
            else:
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='air')
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        coil_dict = gs_mesh.get_coils()
        cond_dict = gs_mesh.get_conductors()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict,'ITER_mesh.h5')

    if not os.path.exists('ITER_mesh.h5'):
        try:
            create_mesh()
        except Exception as e:
            print(e)
            mp_q.put(None)
            return

    # --- Set up GS solver (same as run_ITER_case) ---
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=fe_order, F0=5.3*6.2)

    mygs.set_coil_vsc({'VS': 1.0})
    coil_bounds = {key: [-50.E6, 50.E6] for key in mygs.coil_sets}
    mygs.set_coil_bounds(coil_bounds)

    Ip_target = 15.6E6
    P0_target = 6.2E5
    mygs.set_targets(Ip=Ip_target, pax=P0_target)

    isoflux_pts = np.array([
        [ 8.20,  0.41],
        [ 8.06,  1.46],
        [ 7.51,  2.62],
        [ 6.14,  3.78],
        [ 4.51,  3.02],
        [ 4.26,  1.33],
        [ 4.28,  0.08],
        [ 4.49, -1.34],
        [ 7.28, -1.89],
        [ 8.00, -0.68]
    ])
    x_point = np.array([[5.125, -3.4],])
    mygs.set_isoflux(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddles(x_point)

    regularization_terms = []
    for name in mygs.coil_sets:
        if name.startswith('CS'):
            if name.startswith('CS1'):
                regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=2.E-2))
            else:
                regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
        elif name.startswith('PF'):
            regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
        elif name.startswith('VS'):
            regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=regularization_terms)

    ffp_prof = create_power_flux_fun(40, 1.5, 2.0)
    pp_prof = create_power_flux_fun(40, 4.0, 1.0)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    R0 = 6.3
    Z0 = 0.5
    a = 2.0
    kappa = 1.4
    delta = 0.0
    try:
        mygs.init_psi(R0, Z0, a, kappa, delta)
        mygs.solve()
    except ValueError:
        mp_q.put(None)
        return

    # --- Define kinetic and current profiles for bootstrap solve ---
    n_sample = 257
    psi_sample = np.linspace(0.0, 1.0, n_sample)
    psi_pad = 1.E-3

    # Inductive j_phi profile shape
    jphi_prof = create_power_flux_fun(len(psi_sample), 2.25, 2.5)
    inductive_jphi = jphi_prof['y']

    # H-mode kinetic profiles
    xphalf = 0.965
    widthp_Te = 0.1
    widthp_ne = 0.35

    ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n_sample,
                        expin=1.6, expout=1.6, widthp=widthp_ne, xphalf=xphalf) * 1e20
    Te = Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=n_sample,
                        expin=1.3, expout=1.7, widthp=widthp_Te, xphalf=xphalf)
    ni = ne.copy()       # Assuming quasineutrality
    Ti = Te.copy()       # Assuming isothermal
    Zeff = np.full(n_sample, 1.7)

    # --- Solve with bootstrap current ---
    try:
        bs_results = solve_with_bootstrap(
            mygs,
            ne, Te, ni, Ti, Zeff,
            Ip_target,
            inductive_jphi,
            scale_jBS=1.0,
            isolate_edge_jBS=False,
            psi_pad=psi_pad,
            iterations=2,
            diagnostic_plots=False,
            use_python_solve=True,
        )
    except Exception as e:
        print("Bootstrap solve failed: {0}".format(e))
        mp_q.put(None)
        return

    # --- Collect results ---
    eq_info = mygs.get_stats(li_normalization='ITER')

    # Bootstrap-specific diagnostics
    j_BS = bs_results['j_BS']
    j_total = bs_results['total_j_phi']
    j_ind = bs_results['j_inductive']

    eq_info['j_BS_max'] = float(np.max(np.abs(j_BS)))
    eq_info['j_BS_axis'] = float(j_BS[0])
    eq_info['jphi_axis'] = float(j_total[0])
    eq_info['jphi_max'] = float(np.max(np.abs(j_total)))
    eq_info['j_ind_axis'] = float(j_ind[0])

    # Bootstrap fraction (psi-space trapezoid estimate)
    bs_frac = np.trapezoid(j_BS, psi_sample) / np.trapezoid(j_total, psi_sample) \
              if np.trapezoid(j_total, psi_sample) != 0 else 0.0
    eq_info['bs_fraction'] = float(bs_frac)

    mp_q.put([eq_info])
    oftpy_dump_cov()

# -----------------------------------------------------------------------
# Expected values dictionary
# -----------------------------------------------------------------------
ITER_bootstrap_eq_dict = {
    'Ip': 15599997.261988742,
    'kappa': 1.8746806271900014,
    'R_geo': 6.222524490498655,
    'a_geo': 1.9807072056151687,
    'q_0': 1.036565329524374,
    'q_95': 2.863029898626678,
    'P_ax': 740023.5117187898,
    'j_BS_max': 219163.55709184994,
    'j_BS_axis': 4696.213397541225,
    'jphi_axis': 1397130.088093367,
    'jphi_max': 1507460.349786304,
    # alpha * seed, now = jphi_axis - j_BS_axis: since 1c6955d the Python alpha closure
    # integrates the plasma only (was -7 % from the limiter-area over-count)
    'j_ind_axis': 1.392583E+06,
    'bs_fraction': 0.17910707149956598,
}

@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_ITER_bootstrap(order):
    results = mp_run(run_ITER_bootstrap_case, (1.0, order), timeout=300)
    assert validate_dict(results, ITER_bootstrap_eq_dict)


# -----------------------------------------------------------------------
# Test: redl_bootstrap() directly (same equilibrium as test_ITER_bootstrap)
# -----------------------------------------------------------------------
def run_Redl_jBS_case(mesh_resolution, fe_order, mp_q):
    from OpenFUSIONToolkit.TokaMaker.bootstrap import (
        redl_bootstrap, calculate_ln_lambda, Hmode_profiles
    )

    # --- Mesh creation (identical to run_ITER_bootstrap_case) ---
    def create_mesh():
        with open('ITER_geom.json','r') as fid:
            ITER_geom = json.load(fid)
        plasma_dx = 0.15/mesh_resolution
        coil_dx = 0.2/mesh_resolution
        vv_dx = 0.3/mesh_resolution
        vac_dx = 0.6/mesh_resolution
        gs_mesh = gs_Domain()
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('vacuum1',vv_dx,'vacuum')
        gs_mesh.define_region('vacuum2',vv_dx,'vacuum')
        gs_mesh.define_region('vv1',vv_dx,'conductor',eta=6.9E-7)
        gs_mesh.define_region('vv2',vv_dx,'conductor',eta=6.9E-7)
        for key, coil in ITER_geom['coils'].items():
            if not key.startswith('VS'):
                gs_mesh.define_region(key,coil_dx,'coil')
        gs_mesh.define_region('VSU',coil_dx,'coil',coil_set='VS',nTurns=1.0)
        gs_mesh.define_region('VSL',coil_dx,'coil',coil_set='VS',nTurns=-1.0)
        gs_mesh.add_polygon(ITER_geom['limiter'],'plasma',parent_name='vacuum1')
        gs_mesh.add_annulus(ITER_geom['inner_vv'][0],'vacuum1',ITER_geom['inner_vv'][1],'vv1',parent_name='vacuum2')
        gs_mesh.add_annulus(ITER_geom['outer_vv'][0],'vacuum2',ITER_geom['outer_vv'][1],'vv2',parent_name='air')
        for key, coil in ITER_geom['coils'].items():
            if key.startswith('VS'):
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='vacuum1')
            else:
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='air')
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        coil_dict = gs_mesh.get_coils()
        cond_dict = gs_mesh.get_conductors()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict,'ITER_mesh.h5')

    if not os.path.exists('ITER_mesh.h5'):
        try:
            create_mesh()
        except Exception as e:
            print(e)
            mp_q.put(None)
            return

    # --- Set up GS solver (same as run_ITER_bootstrap_case) ---
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.setup(order=fe_order, F0=5.3*6.2)

    mygs.set_coil_vsc({'VS': 1.0})
    coil_bounds = {key: [-50.E6, 50.E6] for key in mygs.coil_sets}
    mygs.set_coil_bounds(coil_bounds)

    Ip_target = 15.6E6
    P0_target = 6.2E5
    mygs.set_targets(Ip=Ip_target, pax=P0_target)

    isoflux_pts = np.array([
        [ 8.20,  0.41],
        [ 8.06,  1.46],
        [ 7.51,  2.62],
        [ 6.14,  3.78],
        [ 4.51,  3.02],
        [ 4.26,  1.33],
        [ 4.28,  0.08],
        [ 4.49, -1.34],
        [ 7.28, -1.89],
        [ 8.00, -0.68]
    ])
    x_point = np.array([[5.125, -3.4],])
    mygs.set_isoflux(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddles(x_point)

    regularization_terms = []
    for name in mygs.coil_sets:
        if name.startswith('CS'):
            if name.startswith('CS1'):
                regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=2.E-2))
            else:
                regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
        elif name.startswith('PF'):
            regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
        elif name.startswith('VS'):
            regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=regularization_terms)

    ffp_prof = create_power_flux_fun(40, 1.5, 2.0)
    pp_prof = create_power_flux_fun(40, 4.0, 1.0)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    R0 = 6.3
    Z0 = 0.5
    a = 2.0
    kappa = 1.4
    delta = 0.0
    try:
        mygs.init_psi(R0, Z0, a, kappa, delta)
        mygs.solve()
    except ValueError:
        mp_q.put(None)
        return

    # --- Define kinetic profiles (same as run_ITER_bootstrap_case) ---
    EC = 1.602176634e-19
    n_psi = 257
    psi_N = np.linspace(0.0, 1.0, n_psi)
    psi_pad = 1.E-3

    xphalf = 0.965
    widthp_Te = 0.1
    widthp_ne = 0.35

    ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n_psi,
                        expin=1.6, expout=1.6, widthp=widthp_ne, xphalf=xphalf) * 1e20
    Te = Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=n_psi,
                        expin=1.3, expout=1.7, widthp=widthp_Te, xphalf=xphalf)
    ni = ne.copy()
    Ti = Te.copy()
    Zeff = np.full(n_psi, 1.7)

    pressure = (EC * ne * Te) + (EC * ni * Ti)

    # --- Extract geometry from equilibrium (same as solve_with_bootstrap) ---
    psi_eval = np.clip(psi_N, psi_pad, 1.0 - psi_pad)
    _, f, _, _, _ = mygs.get_profiles(psi=psi_eval)
    _, fc, r_avgs, _, eps = mygs.sauter_fc(psi=psi_eval, return_eps=True)

    ft = 1 - fc
    _, qvals, ravgs_q, _, _, _ = mygs.get_q(psi=psi_eval)
    R_avg = ravgs_q['<R>']

    # --- Gradients (same as solve_with_bootstrap) ---
    # Second-order one-sided stencils at the axis/edge, matching the derivative
    # path used by solve_with_bootstrap
    psi_range = mygs.psi_bounds[1] - mygs.psi_bounds[0]
    psi_range_safe = psi_range if psi_range != 0 else 1e-9

    dn_e_dpsi = np.gradient(ne, psi_N, edge_order=2) / psi_range_safe
    dT_e_dpsi = np.gradient(Te, psi_N, edge_order=2) / psi_range_safe
    dn_i_dpsi = np.gradient(ni, psi_N, edge_order=2) / psi_range_safe
    dT_i_dpsi = np.gradient(Ti, psi_N, edge_order=2) / psi_range_safe

    # --- Coulomb logarithms (same as solve_with_bootstrap) ---
    ln_le, ln_lii = calculate_ln_lambda(
        Te, Ti, ne, ni, Zeff,
        electron_lnLambda_model='NRL',
        ion_lnLambda_model='Zavg',
    )

    # --- Collisionalities (same as solve_with_bootstrap) ---
    Zdom = 1.0
    Zavg = ne / ni
    Zion = (Zdom**2 * Zavg * Zeff)**0.25
    nu_i_star = (4.90e-18 * np.abs(qvals) * R_avg * ni
                 * Zion**4 * ln_lii / (Ti**2 * eps**1.5))
    nu_e_star = (6.921e-18 * np.abs(qvals) * R_avg * ne
                 * Zeff * ln_le / (Te**2 * eps**1.5))

    # --- Call redl_bootstrap (same as solve_with_bootstrap) ---
    try:
        j_BS_neo, coeffs = redl_bootstrap(
            psi_N=psi_N, Te=Te, Ti=Ti, ne=ne, ni=ni,
            pe=EC*(ne*Te), pi=EC*(ni*Ti),
            Zeff=Zeff, R=R_avg, q=qvals, eps=eps, fT=ft, I_psi=f,
            dT_e_dpsi=dT_e_dpsi, dT_i_dpsi=dT_i_dpsi,
            dn_e_dpsi=dn_e_dpsi, dn_i_dpsi=dn_i_dpsi,
            ln_lambda_e=ln_le, ln_lambda_ii=ln_lii,
            nu_e_star_override=nu_e_star,
            nu_i_star_override=nu_i_star,
            use_legacy_L34=False,
            use_sign_q=True,
            formula_form='jboot1',
        )
    except Exception as e:
        print("redl_bootstrap failed: {0}".format(e))
        mp_q.put(None)
        return

    # Convert to j_phi (A/m^2) same as solve_with_bootstrap
    j_BS = j_BS_neo * (R_avg / f)
    j_BS = np.nan_to_num(j_BS, nan=0.0)

    # --- Collect results ---
    results = {}
    results['j_BS_max'] = float(np.max(np.abs(j_BS)))
    results['j_BS_axis'] = float(j_BS[0])
    results['j_BS_edge'] = float(j_BS[-1])
    results['L31_axis'] = float(coeffs['L31'][0])
    results['L32_axis'] = float(coeffs['L32'][0])
    results['alpha_axis'] = float(coeffs['alpha'][0])
    results['nu_e_star_axis'] = float(coeffs['nu_e_star'][0])
    results['nu_i_star_axis'] = float(coeffs['nu_i_star'][0])

    mp_q.put([results])
    oftpy_dump_cov()


#============================================================================
# Validation of the optional `x` grid argument to `solve_with_bootstrap`.
# These pin the Python solve path (`use_python_solve=True`); they exercise input
# checking only, which happens before the solver object is touched, so they are
# fast and run in the default CI selection.
@pytest.mark.coverage
def test_bootstrap_x_validation():
    from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap
    n = 65
    ne = np.full(n, 1.0e20)
    Te = np.full(n, 2.0e3)
    Zeff = np.full(n, 1.7)
    def call(x):
        return solve_with_bootstrap(None, ne, Te, ne.copy(), Te.copy(), Zeff,
                                    1.0e6, x=x, use_python_solve=True)
    good = np.linspace(0.0, 1.0, n)
    # wrong length
    with pytest.raises(ValueError, match="same length"):
        call(np.linspace(0.0, 1.0, n-1))
    # duplicated flux label -> undefined derivative
    dup = good.copy(); dup[32] = dup[31]
    with pytest.raises(ValueError, match="strictly increasing"):
        call(dup)
    # unsorted
    unsorted_grid = good.copy(); unsorted_grid[10], unsorted_grid[11] = good[11], good[10]
    with pytest.raises(ValueError, match="strictly increasing"):
        call(unsorted_grid)
    # non-finite
    nan_grid = good.copy(); nan_grid[5] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        call(nan_grid)
    # out of range
    with pytest.raises(ValueError, match=r"within \[0,1\]"):
        call(np.linspace(-0.1, 1.0, n))
    # psi_pad coarser than the grid would collapse distinct flux surfaces
    with pytest.raises(ValueError, match="larger than the first/last"):
        solve_with_bootstrap(None, ne, Te, ne.copy(), Te.copy(), Zeff, 1.0e6,
                             x=good, psi_pad=0.5, use_python_solve=True)
    # deprecated alias psi_N: warns, then validated as x
    with pytest.warns(DeprecationWarning, match="psi_N"):
        with pytest.raises(ValueError, match="same length"):
            solve_with_bootstrap(None, ne, Te, ne.copy(), Te.copy(), Zeff, 1.0e6,
                                 psi_N=good[:-1], use_python_solve=True)
    with pytest.raises(ValueError, match="x only"):
        solve_with_bootstrap(None, ne, Te, ne.copy(), Te.copy(), Zeff, 1.0e6,
                             x=good, psi_N=good, use_python_solve=True)


@pytest.mark.coverage
def test_bootstrap_derivative_edge_order():
    """Profile derivatives must use a 2nd-order stencil at the axis and edge.

    The first-order default of `numpy.gradient` is badly inaccurate at the
    magnetic axis, where it propagates directly into on-axis j_BS.
    """
    psi = np.linspace(0.0, 1.0, 257)
    # analytic profile with a known slope
    y = np.tanh(6.0*(0.9-psi)) + 0.3*np.cos(3.0*psi)
    exact = -6.0/np.cosh(6.0*(0.9-psi))**2 - 0.9*np.sin(3.0*psi)
    d2 = np.gradient(y, psi, edge_order=2)
    d1 = np.gradient(y, psi, edge_order=1)
    # 2nd-order endpoint is dramatically better at the axis
    assert abs(d2[0]-exact[0]) < 0.1*abs(d1[0]-exact[0])
    assert abs(d2[-1]-exact[-1]) < 0.5*abs(d1[-1]-exact[-1])
    # and matches the analytic slope closely in the interior
    assert np.linalg.norm(d2[1:-1]-exact[1:-1])/np.linalg.norm(exact[1:-1]) < 5.e-3


Redl_jBS_eq_dict = {
    'j_BS_max': 190036.48438360557,
    'j_BS_axis': 5101.043999559047,
    'j_BS_edge': 96587.97236259702,
    'L31_axis': 0.10849481301194958,
    'L32_axis': -0.014581267675566556,
    'alpha_axis': -0.6150406171285213,
    'nu_e_star_axis': 0.3420647094964808,
    'nu_i_star_axis': 0.2975729435421919,
}


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_Redl_jBS(order):
    results = mp_run(run_Redl_jBS_case, (1.0, order), timeout=300)
    assert validate_dict(results, Redl_jBS_eq_dict)

#============================================================================
# Internal bootstrap test (ITER-based)
#============================================================================
def run_ITER_bootstrap_case_internal(mesh_resolution, fe_order, mp_q):
    from OpenFUSIONToolkit.TokaMaker.bootstrap import Hmode_profiles

    # --- Mesh creation (same as run_ITER_case) ---
    def create_mesh():
        with open('ITER_geom.json','r') as fid:
            ITER_geom = json.load(fid)
        plasma_dx = 0.15/mesh_resolution
        coil_dx = 0.2/mesh_resolution
        vv_dx = 0.3/mesh_resolution
        vac_dx = 0.6/mesh_resolution
        gs_mesh = gs_Domain()
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('vacuum1',vv_dx,'vacuum')
        gs_mesh.define_region('vacuum2',vv_dx,'vacuum')
        gs_mesh.define_region('vv1',vv_dx,'conductor',eta=6.9E-7)
        gs_mesh.define_region('vv2',vv_dx,'conductor',eta=6.9E-7)
        for key, coil in ITER_geom['coils'].items():
            if not key.startswith('VS'):
                gs_mesh.define_region(key,coil_dx,'coil')
        gs_mesh.define_region('VSU',coil_dx,'coil',coil_set='VS',nTurns=1.0)
        gs_mesh.define_region('VSL',coil_dx,'coil',coil_set='VS',nTurns=-1.0)
        gs_mesh.add_polygon(ITER_geom['limiter'],'plasma',parent_name='vacuum1')
        gs_mesh.add_annulus(ITER_geom['inner_vv'][0],'vacuum1',ITER_geom['inner_vv'][1],'vv1',parent_name='vacuum2')
        gs_mesh.add_annulus(ITER_geom['outer_vv'][0],'vacuum2',ITER_geom['outer_vv'][1],'vv2',parent_name='air')
        for key, coil in ITER_geom['coils'].items():
            if key.startswith('VS'):
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='vacuum1')
            else:
                gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name='air')
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        coil_dict = gs_mesh.get_coils()
        cond_dict = gs_mesh.get_conductors()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict,'ITER_mesh.h5')

    if not os.path.exists('ITER_mesh.h5'):
        try:
            create_mesh()
        except Exception as e:
            print(e)
            mp_q.put(None)
            return

    # --- Kinetic and current profiles (match ITER_Hmode_bootstrap_ex.py) ---
    n_sample = 257
    psi_sample = np.linspace(0.0, 1.0, n_sample)
    Ip_target = 13.0e6
    Zeff_val = 1.5

    xphalf = 0.965
    ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n_sample,
                        expin=1.6, expout=1.6, widthp=0.35, xphalf=xphalf) * 1e20
    Te = Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=n_sample,
                        expin=1.3, expout=1.7, widthp=0.1, xphalf=xphalf)
    ni = ne.copy()
    Ti = Te.copy()

    inductive_jphi = create_power_flux_fun(n_sample, 2.25, 2.5)['y']

    # --- Set up GS solver ---
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.settings.maxits = 100
    mygs.setup(order=fe_order, F0=5.3*6.2)

    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})

    isoflux_pts = np.array([
        [ 8.20,  0.41], [ 8.06,  1.46], [ 7.51,  2.62],
        [ 6.14,  3.78], [ 4.51,  3.02], [ 4.26,  1.33],
        [ 4.28,  0.08], [ 4.49, -1.34], [ 7.28, -1.89],
        [ 8.00, -0.68],
    ])
    x_point = np.array([[5.125, -3.4]])
    mygs.set_isoflux(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddles(x_point)

    regularization_terms = []
    for name in mygs.coil_sets:
        if name.startswith('CS'):
            w = 2.E-2 if name.startswith('CS1') else 1.E-2
        else:
            w = 1.E-2
        regularization_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=w))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.E2))
    mygs.set_coil_reg(reg_terms=regularization_terms)

    # Initial no-bootstrap solve
    mygs.set_targets(Ip=Ip_target, pax=6.2E5)
    mygs.settings.pm = False
    mygs.update_settings()
    try:
        mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
        mygs.solve()
    except ValueError:
        mp_q.put(None)
        return

    # --- Bootstrap solve ---
    try:
        mygs.solve_bootstrap(
            ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': inductive_jphi},
            te_prof={'type': 'linterp', 'x': psi_sample, 'y': Te / 1e3},
            ne_prof={'type': 'linterp', 'x': psi_sample, 'y': ne},
            ti_prof={'type': 'linterp', 'x': psi_sample, 'y': Ti / 1e3},
            ni_prof={'type': 'linterp', 'x': psi_sample, 'y': ni},
            Zeff=Zeff_val,
            Ip_target=Ip_target,
            scale_jBS=1.0,
            diagnose_bs=True,
        )
    except Exception:
        mp_q.put(None)
        return

    mu0 = 4.0 * np.pi * 1e-7
    eq_info = mygs.get_stats(li_normalization='ITER')

    # --- Extract 1D profiles ---
    psi_i, F_i, Fp_i, P_i, Pp_i = mygs.get_profiles(npsi=n_sample, psi_pad=1e-3)
    _, q_i, ravgs_i, _, _, _     = mygs.get_q(npsi=n_sample, psi_pad=1e-3)
    jtor_i = F_i * Fp_i * ravgs_i['<1/R>'] / mu0 + Pp_i * ravgs_i['<R>']

    eq_info['jphi_axis'] = float(jtor_i[0])
    eq_info['jphi_max']  = float(np.max(np.abs(jtor_i)))
    eq_info['q_axis']    = float(q_i[0])
    sample_idx = np.round(np.linspace(0, len(jtor_i) - 1, 10)).astype(int)
    eq_info['jphi_prof'] = [float(jtor_i[i]) for i in sample_idx]

    # --- Verify the <j_BS.B> -> jphi conversion (doc_tokamaker_current_conventions eqs. A7, A8)
    #     against geometry traced independently from the converged equilibrium ---
    try:
        bp = mygs.get_boot_profs()
        inner = (bp['psi_n'] > 0.05) & (bp['psi_n'] < 0.95)
        psi_c = bp['psi_n'][inner]
        _, F_c, Fp_c, _, Pp_c = mygs.get_profiles(psi=psi_c)
        _, _, ravgs_c, _, _, _ = mygs.get_q(psi=psi_c)
        _, _, _, modb_c = mygs.sauter_fc(psi=psi_c)
        R_c, invR_c, B2_c = ravgs_c['<R>'], ravgs_c['<1/R>'], modb_c[1]
        jdotb = bp['jdotb_bs_raw'][inner]
        field_aligned = jdotb * F_c * invR_c / B2_c
        p_term = Pp_c * (R_c - F_c**2 * invR_c / B2_c)
        scale = np.max(np.abs(bp['j_bs_raw'][inner]))
        err = np.max(np.abs(bp['j_bs_raw'][inner] - (field_aligned + p_term))) / scale
        err_flip = np.max(np.abs(bp['j_bs_raw'][inner] - (field_aligned - p_term))) / scale
        print(f"j_bs_raw vs A7: rel err {err:.3e} (pressure term flipped: {err_flip:.3e}), "
              f"max|P'G|/max|j_BS| = {np.max(np.abs(p_term))/scale:.3e}")
        if not (err < 1.0e-2 and err_flip > 5.0 * err):
            raise AssertionError(f"j_bs_raw does not match eq. A7 (rel err {err:.3e}, flipped {err_flip:.3e})")
        # A8: the equilibrium's own <J.B> = F P' + F'<B^2>/mu0 equals the sum of its components'
        jdotb_eq = F_c * Pp_c + Fp_c * B2_c / mu0
        jdotb_sum = jdotb + (bp['j_ind_final'][inner] + bp['jphi_fixed'][inner]) * B2_c / (F_c * invR_c)
        err_par = np.max(np.abs(jdotb_eq - jdotb_sum)) / np.max(np.abs(jdotb_eq))
        print(f"<J.B> balance (A8): rel err {err_par:.3e}")
        if err_par > 2.0e-2:
            raise AssertionError(f"equilibrium <J.B> != sum of components (rel err {err_par:.3e})")
        # A9d: I_p from TokaMaker's own jphi equals the FEM I_p; the old <R><1/R> measure does not
        psi_f = np.linspace(1.0e-4, 1.0 - 1.0e-4, 1001)
        _, F_f, Fp_f, _, Pp_f = mygs.get_profiles(psi=psi_f)
        _, _, rv_f, _, _, _ = mygs.get_q(psi=psi_f)
        J_f = F_f * Fp_f * rv_f['<1/R>'] / mu0 + Pp_f * rv_f['<R>']
        psi_phys = mygs.psi_bounds[0] + psi_f * (mygs.psi_bounds[1] - mygs.psi_bounds[0])
        w = rv_f['dV/dPsi'] / (2.0 * np.pi)
        ip_exact = abs(np.trapezoid(w * (J_f * rv_f['<1/R^2>'] / rv_f['<1/R>']
                   + Pp_f * (1.0 - rv_f['<R>'] * rv_f['<1/R^2>'] / rv_f['<1/R>'])), psi_phys))
        ip_qtmp = abs(np.trapezoid(w * J_f / rv_f['<R>'], psi_phys))
        print(f"I_p: FEM {eq_info['Ip']:.6e}, A9d {ip_exact:.6e} ({ip_exact/eq_info['Ip']-1:+.3e}), "
              f"old <R><1/R> measure {ip_qtmp:.6e} ({ip_qtmp/eq_info['Ip']-1:+.3e})")
        if abs(ip_exact / eq_info['Ip'] - 1.0) > 2.0e-3:
            raise AssertionError(f"A9d I_p {ip_exact:.6e} != FEM I_p {eq_info['Ip']:.6e}")
    except Exception:
        import traceback
        traceback.print_exc()
        mp_q.put(None)
        return

    # --- Verify that boot_ops, boot_profs round-trips correctly through save/load, and that
    #     replace_eq(source_file=...) correctly syncs the _boot_ops shadow dict ---
    save_file = 'ITER_boot_ops_test.h5'
    try:
        # Capture expected state before save so we can check all fields
        expected_boot_ops = dict(mygs._tMaker_equil._boot_ops)
        expected_boot_profs = mygs.get_boot_profs()
        mygs._tMaker_equil.save_TokaMaker(save_file)
        # Corrupt the shadow dict so we can confirm replace_eq overwrites it from the file
        mygs._tMaker_equil._boot_ops['scale_jBS'] = -999.0
        mygs.replace_eq(source_file=save_file)
        boot_ops = mygs._tMaker_equil._boot_ops
        if boot_ops is None:
            raise AssertionError("_boot_ops is None after replace_eq(source_file=...)")
        # Verify all fields round-trip correctly through save/load
        for key, expected in expected_boot_ops.items():
            val = boot_ops[key]
            if isinstance(expected, bool):
                if val != expected:
                    raise AssertionError(
                        f"_boot_ops['{key}'] = {val} != {expected} after replace_eq(source_file=...)"
                    )
            elif isinstance(expected, float):
                if abs(val - expected) > 1e-10:
                    raise AssertionError(
                        f"_boot_ops['{key}'] = {val} != {expected} after replace_eq(source_file=...)"
                    )
            elif isinstance(expected, int):
                if val != expected:
                    raise AssertionError(
                        f"_boot_ops['{key}'] = {val} != {expected} after replace_eq(source_file=...)"
                    )
        # Verify BOOT_PROFS arrays round-trip correctly through save/load
        if expected_boot_profs is None:
            raise AssertionError("get_boot_profs() returned None before save")
        boot_profs = mygs.get_boot_profs()
        if boot_profs is None:
            raise AssertionError("get_boot_profs() returned None after replace_eq(source_file=...)")
        for key, expected_arr in expected_boot_profs.items():
            if key not in boot_profs:
                raise AssertionError(
                    f"boot_profs key '{key}' missing after replace_eq(source_file=...)"
                )
            if not np.allclose(boot_profs[key], expected_arr, rtol=1e-3):
                reldiff = np.abs(boot_profs[key] - expected_arr) / (np.abs(expected_arr) + 1e-30)
                idx = int(np.argmax(reldiff))
                raise AssertionError(
                    f"boot_profs['{key}'] does not match after replace_eq(source_file=...) "
                    f"max_reldiff={reldiff.max():.3e} at idx={idx} "
                    f"(expected={expected_arr[idx]:.6e}, got={boot_profs[key][idx]:.6e})"
                )
    except Exception as e:
        import traceback
        traceback.print_exc()
        mp_q.put(None)
        return

    # --- Verify that boot_ops, boot_profs round-trips correctly through copy, and that
    #     replace_eq(source_eq=...) correctly syncs the _boot_ops shadow dict ---
    try:
        # Capture expected state before copy so we can check all fields
        expected_boot_ops = dict(mygs._tMaker_equil._boot_ops)
        expected_boot_profs = mygs.get_boot_profs()
        mygs_copied = mygs.copy_eq()
        # Corrupt the shadow dict so we can confirm replace_eq overwrites it from the file
        mygs._tMaker_equil._boot_ops['scale_jBS'] = -999.0
        mygs.replace_eq(source_eq=mygs_copied)
        boot_ops = mygs._tMaker_equil._boot_ops
        if boot_ops is None:
            raise AssertionError("_boot_ops is None after replace_eq(source_eq=...)")
        # Verify all fields round-trip correctly through save/load
        for key, expected in expected_boot_ops.items():
            val = boot_ops[key]
            if isinstance(expected, bool):
                if val != expected:
                    raise AssertionError(
                        f"_boot_ops['{key}'] = {val} != {expected} after replace_eq(source_eq=...)"
                    )
            elif isinstance(expected, float):
                if abs(val - expected) > 1e-10:
                    raise AssertionError(
                        f"_boot_ops['{key}'] = {val} != {expected} after replace_eq(source_eq=...)"
                    )
            elif isinstance(expected, int):
                if val != expected:
                    raise AssertionError(
                        f"_boot_ops['{key}'] = {val} != {expected} after replace_eq(source_eq=...)"
                    )
            print(key,val,expected)
        # Verify BOOT_PROFS arrays round-trip correctly through save/load
        if expected_boot_profs is None:
            raise AssertionError("get_boot_profs() returned None before save")
        boot_profs = mygs.get_boot_profs()
        if boot_profs is None:
            raise AssertionError("get_boot_profs() returned None after replace_eq(source_eq=...)")
        for key, expected_arr in expected_boot_profs.items():
            if key not in boot_profs:
                raise AssertionError(
                    f"boot_profs key '{key}' missing after replace_eq(source_eq=...)"
                )
            if not np.allclose(boot_profs[key], expected_arr, rtol=1e-3):
                reldiff = np.abs(boot_profs[key] - expected_arr) / (np.abs(expected_arr) + 1e-30)
                idx = int(np.argmax(reldiff))
                raise AssertionError(
                    f"boot_profs['{key}'] does not match after replace_eq(source_eq=...) "
                    f"max_reldiff={reldiff.max():.3e} at idx={idx} "
                    f"(expected={expected_arr[idx]:.6e}, got={boot_profs[key][idx]:.6e})"
                )
    except Exception as e:
        import traceback
        traceback.print_exc()
        mp_q.put(None)
        return

    # --- verify that scalar Zeff and a linearly-increasing Zeff profile
    #     produce different bootstrap current profiles ---
    try:
        zeff_common_kwargs = dict(
            ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': inductive_jphi},
            te_prof={'type': 'linterp', 'x': psi_sample, 'y': Te / 1e3},
            ne_prof={'type': 'linterp', 'x': psi_sample, 'y': ne},
            ti_prof={'type': 'linterp', 'x': psi_sample, 'y': Ti / 1e3},
            ni_prof={'type': 'linterp', 'x': psi_sample, 'y': ni},
            Ip_target=Ip_target,
            scale_jBS=1.0,
        )
        profs_scalar = mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        profs_linear = mygs.solve_bootstrap(
            Zeff={'x': psi_sample, 'y': np.linspace(1.0, 2.5, n_sample)},
            **zeff_common_kwargs,
        )
        j_scalar = profs_scalar['j_bs_raw']
        j_linear = profs_linear['j_bs_raw']
        magnitude = 0.5 * (np.abs(j_scalar) + np.abs(j_linear))
        rel_diff = np.where(magnitude > 0, (j_linear - j_scalar) / magnitude, 0.0)
        print("\nZeff scalar vs linear j_bs_raw relative difference (j_linear-j_scalar)/|mean|:")
        print(rel_diff)
        if np.allclose(j_scalar, j_linear):
            raise AssertionError(
                "Bootstrap profiles with scalar Zeff and linearly-increasing Zeff profile "
                "are identical; expected them to differ."
            )
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    # --- verify jphi_fixed_prof is added unscaled to the total and reduces the
    #     inductive share, and that omitting it afterwards resets it to zero ---
    try:
        jfix_y = 2.0e5 * np.exp(-((psi_sample - 0.5) / 0.1)**2)
        profs_fixed = mygs.solve_bootstrap(
            Zeff=Zeff_val,
            jphi_fixed_prof={'type': 'linterp', 'x': psi_sample, 'y': jfix_y},
            **zeff_common_kwargs,
        )
        if not np.allclose(profs_fixed['jphi_fixed'], np.interp(profs_fixed['psi_n'], psi_sample, jfix_y),
                           rtol=1e-6, atol=1e-6*jfix_y.max()):
            raise AssertionError("boot_profs['jphi_fixed'] does not match the input jphi_fixed_prof")
        j_sum = profs_fixed['j_ind_final'] + profs_fixed['j_bs_final'] + profs_fixed['jphi_fixed']
        if not np.allclose(profs_fixed['total_j_phi'], j_sum, rtol=1e-10, atol=1e-6):
            raise AssertionError("total_j_phi != j_ind_final + j_bs_final + jphi_fixed")
        if not np.sum(profs_fixed['j_ind_final']) < np.sum(profs_scalar['j_ind_final']):
            raise AssertionError("jphi_fixed did not reduce the inductive current share")
        profs_reset = mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        if np.any(profs_reset['jphi_fixed'] != 0.0):
            raise AssertionError("jphi_fixed persisted after solve_bootstrap without jphi_fixed_prof")
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    # --- verify p_fixed_prof [Pa] is added to the kinetic pressure: same result as
    #     pres_prof = kinetic + p_fixed, and mutually exclusive with pres_prof ---
    try:
        p_kin = eC * (ne * Te + ni * Ti)
        # Offset keeps pres_prof safely above the kinetic pressure (roundoff) at the edge
        pf = p_kin[0] * (0.2 * np.exp(-((psi_sample - 0.3) / 0.15)**2) + 1.0e-3)
        profs_pfix = mygs.solve_bootstrap(Zeff=Zeff_val, p_fixed_prof={'x': psi_sample, 'y': pf}, **zeff_common_kwargs)
        _, _, _, P_pfix, _ = mygs.get_profiles(npsi=n_sample)
        P_ax_pfix = np.max(mygs.get_profiles(psi=np.array([0.0, 1.0]))[3])
        profs_pres = mygs.solve_bootstrap(Zeff=Zeff_val, pres_prof={'x': psi_sample, 'y': p_kin + pf}, **zeff_common_kwargs)
        _, _, _, P_pres, _ = mygs.get_profiles(npsi=n_sample)
        # Both solves converge j_BS to djBS_tol (1e-4) and may stop an iteration apart
        for key in profs_pres:
            if not np.allclose(profs_pfix[key], profs_pres[key], rtol=2e-4, atol=2e-4*np.max(np.abs(profs_pres[key]))):
                raise AssertionError(f"p_fixed_prof vs pres_prof=kinetic+p_fixed: boot_profs['{key}'] differ")
        if not np.allclose(P_pfix, P_pres, rtol=1e-6, atol=1e-6*np.max(P_pres)):
            raise AssertionError("p_fixed_prof vs pres_prof=kinetic+p_fixed: pressure profiles differ")
        if not np.isclose(P_ax_pfix, p_kin[0] + pf[0], rtol=1e-4):
            raise AssertionError(f"axis pressure {P_ax_pfix:.6e} != kinetic + p_fixed {p_kin[0] + pf[0]:.6e}")
        mygs.solve_bootstrap(Zeff=Zeff_val, **zeff_common_kwargs)
        if np.isclose(np.max(mygs.get_profiles(psi=np.array([0.0, 1.0]))[3]), P_ax_pfix, rtol=1e-4):
            raise AssertionError("p_fixed_prof did not change the axis pressure")
        for bad_kwargs, label in (({'pres_prof': {'x': psi_sample, 'y': p_kin + pf}}, "with pres_prof"),
                                  ({}, "negative")):
            try:
                y_bad = -pf if label == "negative" else pf
                mygs.solve_bootstrap(Zeff=Zeff_val, p_fixed_prof={'x': psi_sample, 'y': y_bad}, **bad_kwargs, **zeff_common_kwargs)
            except ValueError:
                continue
            raise AssertionError(f"p_fixed_prof {label} did not raise ValueError")
    except Exception as e:
        print(e)
        mp_q.put(None)
        return

    mp_q.put([eq_info])
    oftpy_dump_cov()


ITER_bootstrap_internal_eq_dict = {
    'Ip':        13000003.24961143,
    'kappa':     1.8710874971500062,
    'R_geo':     6.2229574417901645,
    'a_geo':     1.9785979198829686,
    'q_0':       1.3748401369017682,
    'q_95':      3.504083741872661,
    'P_ax':      740021.5037035183,
    'beta_pol':  85.19729428055419,
    'beta_tor':  2.492084824025585,
    'jphi_axis': 1010416.1199979713,
    'jphi_max':  1173298.5135544762,
    'q_axis':    1.410565702638605,
    'jphi_prof': [1010416.1199979713, 1155117.92412601,   1162817.6369466858,
                  1072123.3558973635,  898492.2764938722,  683025.6545245483,
                   449301.00087938644, 249118.67207980127, 173771.68344307062,
                   136306.8890837856],
}


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_ITER_bootstrap_internal(order):
    results = mp_run(run_ITER_bootstrap_case_internal, (1.0, order), timeout=300)
    assert validate_dict(results, ITER_bootstrap_internal_eq_dict)

# -----------------------------------------------------------------------
# Test: profiles on normalized toroidal flux (coord='phi_n'/'phi_n_relabel')
#
# Each case solves ITER with psi_N profiles, maps them onto Phi_N with the
# converged q, re-solves with the toroidal-flux profiles and compares.
# Differences are set by profile interpolation in the two coordinates.
# -----------------------------------------------------------------------
def _torflux_ITER_gs(fe_order, Ip_target=15.6E6):
    if not os.path.exists('ITER_mesh.h5'):
        with open('ITER_geom.json','r') as fid:
            ITER_geom = json.load(fid)
        gs_mesh = gs_Domain()
        gs_mesh.define_region('air',0.6,'boundary')
        gs_mesh.define_region('plasma',0.15,'plasma')
        gs_mesh.define_region('vacuum1',0.3,'vacuum')
        gs_mesh.define_region('vacuum2',0.3,'vacuum')
        gs_mesh.define_region('vv1',0.3,'conductor',eta=6.9E-7)
        gs_mesh.define_region('vv2',0.3,'conductor',eta=6.9E-7)
        for key, coil in ITER_geom['coils'].items():
            if not key.startswith('VS'):
                gs_mesh.define_region(key,0.2,'coil')
        gs_mesh.define_region('VSU',0.2,'coil',coil_set='VS',nTurns=1.0)
        gs_mesh.define_region('VSL',0.2,'coil',coil_set='VS',nTurns=-1.0)
        gs_mesh.add_polygon(ITER_geom['limiter'],'plasma',parent_name='vacuum1')
        gs_mesh.add_annulus(ITER_geom['inner_vv'][0],'vacuum1',ITER_geom['inner_vv'][1],'vv1',parent_name='vacuum2')
        gs_mesh.add_annulus(ITER_geom['outer_vv'][0],'vacuum2',ITER_geom['outer_vv'][1],'vv2',parent_name='air')
        for key, coil in ITER_geom['coils'].items():
            parent = 'vacuum1' if key.startswith('VS') else 'air'
            gs_mesh.add_rectangle(coil['rc'],coil['zc'],coil['w'],coil['h'],key,parent_name=parent)
        mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
        save_gs_mesh(mesh_pts,mesh_lc,mesh_reg,gs_mesh.get_coils(),gs_mesh.get_conductors(),'ITER_mesh.h5')
    mygs = TokaMaker(OFT_env(nthreads=-1))
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.settings.maxits = 100
    mygs.setup(order=fe_order,F0=5.3*6.2)
    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})
    mygs.set_targets(Ip=Ip_target, pax=6.2E5)
    isoflux_pts = np.array([
        [ 8.20,  0.41], [ 8.06,  1.46], [ 7.51,  2.62], [ 6.14,  3.78], [ 4.51,  3.02],
        [ 4.26,  1.33], [ 4.28,  0.08], [ 4.49, -1.34], [ 7.28, -1.89], [ 8.00, -0.68]
    ])
    x_point = np.array([[5.125, -3.4],])
    mygs.set_isoflux(np.vstack((isoflux_pts,x_point)))
    mygs.set_saddles(x_point)
    reg_terms = [mygs.coil_reg_term({name: 1.0},target=0.0,weight=(2.E-2 if name.startswith('CS1') else 1.E-2))
                 for name in mygs.coil_sets]
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E2))
    mygs.set_coil_reg(reg_terms=reg_terms)
    return mygs


def _torflux_phi_map(mygs):
    '''Phi_N(psi_N) and dPhi_N/dpsi_N from the converged q profile'''
    from scipy.integrate import cumulative_trapezoid
    psi_n, q, _, _, _, _ = mygs.get_q(npsi=400,psi_pad=1.E-4)
    phi = cumulative_trapezoid(q,psi_n,initial=0.0)
    return psi_n, phi/phi[-1], q/phi[-1]


def _torflux_to_phi(x, pmap):
    x_phi = np.interp(x,pmap[0],pmap[1])
    x_phi[0] = 0.0; x_phi[-1] = 1.0
    return x_phi


def _torflux_compare(mygs, ref, its=None):
    psi = mygs.get_psi(False)
    _, q, _, _, _, _ = mygs.get_q(npsi=100,psi_pad=1.E-3)
    stats = mygs.get_stats()
    return {
        'psi_err': np.linalg.norm(psi-ref['psi'])/np.linalg.norm(ref['psi']),
        'q_err': np.max(np.abs(q-ref['q'])/np.abs(ref['q'])),
        'Ip_err': abs(stats['Ip']-ref['Ip'])/ref['Ip'],
        'its': its
    }


def _torflux_ref(mygs, its=None):
    _, q, _, _, _, _ = mygs.get_q(npsi=100,psi_pad=1.E-3)
    return {'psi': mygs.get_psi(False), 'q': q, 'Ip': mygs.get_stats()['Ip'], 'its': its}


def run_ITER_torflux_case(fe_order, test_type, mp_q):
    try:
        results = {}
        if test_type == 'bootstrap':
            from OpenFUSIONToolkit.TokaMaker.bootstrap import Hmode_profiles, solve_with_bootstrap
            mygs = _torflux_ITER_gs(fe_order,Ip_target=13.0E6)
        else:
            mygs = _torflux_ITER_gs(fe_order)
        init = lambda: mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
        pp = create_power_flux_fun(40,4.0,1.0)
        if test_type == 'ffp':
            # FF' and p' as linterp, both toroidal modes, then copy and save/load
            ffp = create_power_flux_fun(40,1.5,2.0)
            mygs.set_profiles(ffp_prof=ffp,pp_prof=pp)
            init()
            _, its = mygs.solve(return_its=True)
            ref = _torflux_ref(mygs,its)
            pmap = _torflux_phi_map(mygs)
            for coord in ('phi_n', 'phi_n_relabel'):
                profs = {}
                for key, prof in (('ffp', ffp), ('pp', pp)):
                    y = prof['y']/np.interp(prof['x'],pmap[0],pmap[2]) if coord == 'phi_n' else prof['y']
                    profs[key] = {'type': 'linterp', 'x': _torflux_to_phi(prof['x'],pmap), 'y': y, 'coord': coord}
                mygs.set_profiles(ffp_prof=profs['ffp'],pp_prof=profs['pp'])
                init()
                _, its = mygs.solve(return_its=True)
                results[coord] = _torflux_compare(mygs,ref,its)
            results['psi_its'] = ref['its']
            # Solver map against q: dPhi_N/dpsi_N = q/Q and Phi_N increments = int q dpsi_N/Q
            from scipy.integrate import cumulative_trapezoid
            x = np.linspace(0.05,0.95,19)
            phi_x, jac_x = mygs.get_torflux_map(x)
            psi_x, _ = mygs.get_torflux_map(phi_x,inverse=True)
            psi_q, q, _, _, _, _ = mygs.get_q(psi=np.linspace(0.05,0.95,400))
            Qinv = jac_x/np.interp(x,psi_q,q)
            dphi = cumulative_trapezoid(q,psi_q,initial=0.0)*np.mean(Qinv)
            results['jac_err'] = np.max(np.abs(Qinv/np.mean(Qinv)-1.0))
            results['map_err'] = np.max(np.abs(phi_x-phi_x[0]-np.interp(x,psi_q,dphi)))
            results['inv_err'] = np.max(np.abs(psi_x-x))
            # Copy and save/load keep (or rebuild) the map
            prof0 = np.array(mygs.get_profiles(npsi=60)[1:])
            eq_copy = mygs.copy_eq()
            prof_c = np.array(eq_copy.get_profiles(npsi=60)[1:])
            results['copy_err'] = np.max(np.abs(prof_c-prof0))/np.max(np.abs(prof0))
            eq_copy.save_TokaMaker('torflux_eq.h5')
            init()
            mygs.replace_eq(source_file='torflux_eq.h5')
            prof_l = np.array(mygs.get_profiles(npsi=60)[1:])
            results['load_err'] = np.max(np.abs(prof_l-prof0))/np.max(np.abs(prof0))
        elif test_type == 'jphi':
            # jphi-linterp current profile (relabel only)
            x = np.linspace(0.0,1.0,65)
            jphi = {'type': 'jphi-linterp', 'x': x, 'y': create_power_flux_fun(65,2.25,2.5)['y']}
            mygs.set_profiles(ffp_prof=jphi,pp_prof=pp)
            init()
            _, its = mygs.solve(return_its=True)
            ref = _torflux_ref(mygs,its)
            pmap = _torflux_phi_map(mygs)
            mygs.set_profiles(ffp_prof=dict(jphi,x=_torflux_to_phi(x,pmap),coord='phi_n_relabel'),pp_prof=pp)
            init()
            _, its = mygs.solve(return_its=True)
            results['phi_n_relabel'] = _torflux_compare(mygs,ref,its)
            results['psi_its'] = ref['its']
            try:
                mygs.set_profiles(ffp_prof=dict(jphi,x=_torflux_to_phi(x,pmap),coord='phi_n'),pp_prof=pp)
                results['phi_n_rejected'] = False
            except ValueError:
                results['phi_n_rejected'] = True
        elif test_type == 'bootstrap':
            # solve_bootstrap with every profile (jphi, kinetics, pressure, jphi_fixed, p_fixed) on Phi_N
            n = 257
            x = np.linspace(0.0,1.0,n)
            ne = Hmode_profiles(edge=0.35,ped=0.6,core=1.1,rgrid=n,expin=1.6,expout=1.6,widthp=0.35,xphalf=0.965)*1.E20
            Te = Hmode_profiles(edge=1500.,ped=5000.,core=21000.,rgrid=n,expin=1.3,expout=1.7,widthp=0.1,xphalf=0.965)
            jind = create_power_flux_fun(n,2.25,2.5)['y']
            jfix = 1.E5*np.exp(-((x-0.5)/0.15)**2)
            pfix = 0.2*2.0*1.602176634E-19*ne[0]*Te[0]*np.exp(-((x-0.3)/0.15)**2)  # [Pa]
            mygs.set_profiles(ffp_prof=create_power_flux_fun(40,1.5,2.0),pp_prof=pp)
            init()
            mygs.solve()
            def run(x_in, coord):
                return mygs.solve_bootstrap(
                    ffp_prof={'type': 'jphi-split-bootstrap', 'x': x_in, 'y': jind},
                    te_prof={'type': 'linterp', 'x': x_in, 'y': Te/1.E3},
                    ne_prof={'type': 'linterp', 'x': x_in, 'y': ne},
                    ti_prof={'type': 'linterp', 'x': x_in, 'y': Te/1.E3},
                    ni_prof={'type': 'linterp', 'x': x_in, 'y': ne},
                    jphi_fixed_prof={'type': 'linterp', 'x': x_in, 'y': jfix},
                    p_fixed_prof={'x': x_in, 'y': pfix},
                    Zeff=1.5, Ip_target=13.0E6, coord=coord)
            res_psi = run(x,'psi_n')
            ref = _torflux_ref(mygs)
            x_phi = _torflux_to_phi(x,_torflux_phi_map(mygs))
            res_phi = run(x_phi,'phi_n')
            results['phi_n'] = _torflux_compare(mygs,ref)
            # Bootstrap current at matching surfaces (interior, away from edge spike)
            j_ref = np.interp(res_phi['psi_n'],res_psi['psi_n'],res_psi['j_bs_final'])
            mask = (res_phi['psi_n'] > 0.05) & (res_phi['psi_n'] < 0.9)
            results['jbs_err'] = np.max(np.abs(res_phi['j_bs_final'][mask]-j_ref[mask]))/np.max(np.abs(j_ref))
            # jphi_fixed relabelled onto Phi_N nodes returns the input values
            results['jfix_err'] = np.max(np.abs(res_phi['jphi_fixed']-jfix))/np.max(jfix)
            # Node psi_N positions agree with the solver map
            psi_nodes, _ = mygs.get_torflux_map(x_phi,inverse=True)
            results['node_err'] = np.max(np.abs(res_phi['psi_n']-psi_nodes))
            # solve_with_bootstrap: 'psi_n' output; Python solver rejects toroidal coordinates
            res_swb = solve_with_bootstrap(mygs,ne,Te,ne,Te,1.5,13.0E6,inductive_jphi=jind,x=x_phi,coord='phi_n',
                                           jphi_fixed=jfix,p_fixed=pfix)
            psi_nodes, _ = mygs.get_torflux_map(x_phi,inverse=True)
            results['swb_node_err'] = np.max(np.abs(res_swb['psi_n']-psi_nodes))
            results['swb_jfix_err'] = np.max(np.abs(res_swb['j_fixed']-jfix))/np.max(jfix)
            try:
                solve_with_bootstrap(mygs,ne,Te,ne,Te,1.5,13.0E6,inductive_jphi=jind,x=x_phi,
                                     use_python_solve=True,coord='phi_n')
                results['python_rejected'] = False
            except ValueError:
                results['python_rejected'] = True
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    mp_q.put(results)
    oftpy_dump_cov()


def _torflux_check(res, psi_tol, q_tol, its_ref=None):
    assert res['psi_err'] < psi_tol
    assert res['q_err'] < q_tol
    # I_p is set by the exact 1-D profile quadrature (no FEM rescale in jphi_bs_update): residual up to ~4e-4
    assert res['Ip_err'] < 5.E-4
    if its_ref is not None:
        assert res['its'] <= its_ref + 5


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,))
def test_ITER_torflux_ffp(order):
    results = mp_run(run_ITER_torflux_case,(order,'ffp'),timeout=120)
    assert results is not None
    _torflux_check(results['phi_n'],4.E-3,4.E-3,results['psi_its'])
    _torflux_check(results['phi_n_relabel'],2.E-3,3.E-3,results['psi_its'])
    assert results['map_err'] < 1.E-4
    assert results['jac_err'] < 1.E-3
    assert results['inv_err'] < 1.E-10
    assert results['copy_err'] < 1.E-12
    assert results['load_err'] < 5.E-5


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,))
def test_ITER_torflux_jphi(order):
    results = mp_run(run_ITER_torflux_case,(order,'jphi'),timeout=120)
    assert results is not None
    _torflux_check(results['phi_n_relabel'],1.E-3,2.E-3,results['psi_its'])
    assert results['phi_n_rejected']


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_ITER_torflux_bootstrap(order):
    results = mp_run(run_ITER_torflux_case,(order,'bootstrap'),timeout=600)
    assert results is not None
    _torflux_check(results['phi_n'],3.E-3,3.E-3)
    assert results['jbs_err'] < 1.E-2
    assert results['node_err'] < 1.E-12
    assert results['swb_node_err'] < 1.E-12
    assert results['jfix_err'] < 1.E-6
    assert results['swb_jfix_err'] < 1.E-6
    assert results['python_rejected']

# -----------------------------------------------------------------------
# Test: GEQDSK (g-file) reader in TokaMaker.eqdsk
#
# Expected values and per-key tolerance overrides live next to the test
# as JSON fixtures (`eqdsk_gfile_expected.json`, `eqdsk_gfile_tol.json`).
# They are produced by a one-time offline run of the OMFIT `OMFITgeqdsk`
# reader (raw parser + fluxSurfaces tracer) on the same ITER_test.eqdsk.
# The TokaMaker reader was adapted from OMFIT, so:
#   - Raw scalar/profile quantities match OMFIT to floating-point precision
#     since both parsers read the same ASCII file (default 1% tolerance
#     in validate_dict is orders of magnitude more than needed).
#   - Derived flux-surface quantities (q, geometry, li, betas) come from
#     OMFIT's fluxSurfaces; per-key tolerances are relaxed in the `_tol`
#     fixture where the TokaMaker implementation uses a slightly different
#     algorithm (e.g. OMFIT uses a 4-point triangularity definition, we
#     use the midplane-point one; edge resampling differs; etc.).
#   - j_tor_averaged (<Jt/R>/<1/R>), j_tor_averaged_direct
#     (-(p'<R> + FF'<1/R>/mu0)*(2pi)^exp_Bp, the TokaMaker default), and
#     q_profile are each sampled at five interior psi_N positions with
#     tight (0.5%) tolerance.  These are the core analytic Grad-Shafranov
#     outputs used by bootstrap and reconstruction; alignment with OMFIT
#     is ~0.004% RMS for the standard j_tor and ~0.01% for j_tor_direct
#     when the machinery is healthy.
# Regenerate with tests/physics/generate_eqdsk_expected.py in an env with
# scipy<1.14 / numpy<2.0 (e.g. the `omfit_env` conda env); the regenerator
# rewrites both JSON fixtures in place.
# -----------------------------------------------------------------------
def _load_eqdsk_fixture(filename):
    with open(os.path.join(test_dir, filename)) as fh:
        return json.load(fh)


eqdsk_gfile_expected_dict = _load_eqdsk_fixture('eqdsk_gfile_expected.json')
eqdsk_gfile_tol_dict      = _load_eqdsk_fixture('eqdsk_gfile_tol.json')


def run_eqdsk_gfile_case():
    """Read ITER_test.eqdsk with TokaMaker.eqdsk and produce a result dict
    with the same keys as eqdsk_gfile_expected_dict."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    raw = eq._raw
    results = {}

    # --- Raw scalars (pulled directly from the parsed dict) ---
    for key in ['NW','NH','RDIM','ZDIM','RCENTR','RLEFT','ZMID',
                'RMAXIS','ZMAXIS','SIMAG','SIBRY','BCENTR','CURRENT']:
        results[key] = float(raw[key])

    # --- 1-D profiles: spot samples + sum + length ---
    for key in ['FPOL','PRES','FFPRIM','PPRIME','QPSI']:
        arr = np.asarray(raw[key])
        results[f'{key}_0']   = float(arr[0])
        results[f'{key}_mid'] = float(arr[len(arr)//2])
        results[f'{key}_end'] = float(arr[-1])
        results[f'{key}_sum'] = float(np.sum(arr))
        results[f'{key}_len'] = len(arr)

    # --- 2-D PSIRZ spot values ---
    psi = np.asarray(raw['PSIRZ'])
    results['PSIRZ_min'] = float(np.min(psi))
    results['PSIRZ_max'] = float(np.max(psi))
    results['PSIRZ_sum'] = float(np.sum(psi))

    # --- Boundary and limiter ---
    for key in ['RBBBS','ZBBBS','RLIM','ZLIM']:
        arr = np.asarray(raw[key])
        results[f'{key}_min'] = float(np.min(arr))
        results[f'{key}_max'] = float(np.max(arr))
        results[f'{key}_len'] = len(arr)

    # --- Derived FSA quantities ---
    geo = eq.geometry
    avg = eq.averages
    mid = eq.midplane
    li  = eq.li
    betas = eq.betas
    results['Ip_enclosed_edge'] = float(avg['ip'][-1])
    results['q_axis']           = float(eq.q_profile[0])
    results['kappa_edge']       = float(geo['kappa'][-1])
    results['delta_edge']       = float(geo['delta'][-1])
    results['a_edge']           = float(geo['a'][-1])
    results['R_geo_edge']       = float(geo['R'][-1])
    results['vol_edge']         = float(geo['vol'][-1])
    results['cxArea_edge']      = float(geo['cxArea'][-1])
    results['li_1']             = float(li['li(1)'])
    results['li_3']             = float(li['li(3)'])
    results['R_mid_edge']       = float(mid['R'][-1])
    results['Btot_mid_edge']    = float(mid['Btot'][-1])
    results['beta_t']           = float(betas.get('beta_t', 0.0))
    results['beta_p']           = float(betas.get('beta_p', 0.0))
    results['beta_n']           = float(betas.get('beta_n', 0.0))
    results['cocos']            = int(eq.cocos)

    # --- j_tor and q profile sampling at interior psi_N positions ---
    # FSA quantities live on psi_N_levels (length nlevels), which matches the
    # raw g-file psi_N grid by default (nlevels=NW) but may differ if a
    # custom nlevels is supplied at construction.
    psi_N_lvl  = eq.psi_N_levels
    j_tor      = eq.j_tor_averaged         # <Jt/R>/<1/R>  standard convention
    j_tor_dir  = eq.j_tor_averaged_direct  # -(p'<R> + FF'<1/R>/mu0)*(2pi)^exp_Bp
    q          = eq.q_profile
    for pct in (10, 25, 50, 75, 90):
        psin_target = pct / 100.0
        results[f'j_tor_psiN_{pct:02d}']        = float(np.interp(psin_target, psi_N_lvl, j_tor))
        results[f'j_tor_direct_psiN_{pct:02d}'] = float(np.interp(psin_target, psi_N_lvl, j_tor_dir))
        results[f'q_psiN_{pct:02d}']            = float(np.interp(psin_target, psi_N_lvl, q))
    return results


def test_eqdsk_gfile():
    os.chdir(test_dir)
    results = run_eqdsk_gfile_case()
    assert validate_dict([results], eqdsk_gfile_expected_dict,
                         tol_dict=eqdsk_gfile_tol_dict)


# -----------------------------------------------------------------------
# Test: Osborne p-file reader in TokaMaker.eqdsk
#
# Expected values live next to the test as `eqdsk_pfile_expected.json`.
# They are produced by a one-time offline run of the OMFIT `OMFITpFile`
# reader on the same D3Dlike_Hmode_test.peqdsk file.  The TokaMaker reader
# was adapted from the OMFIT code, so raw profile quantities must match
# to floating-point precision.  Regenerate via
# tests/physics/generate_eqdsk_expected.py.
# -----------------------------------------------------------------------
eqdsk_pfile_expected_dict = _load_eqdsk_fixture('eqdsk_pfile_expected.json')


def run_eqdsk_pfile_case():
    """Read D3Dlike_Hmode_test.peqdsk with TokaMaker.eqdsk and produce a
    result dict matching the expected-dict keys."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_pfile
    pf = read_pfile(os.path.join(test_dir, 'D3Dlike_Hmode_test.peqdsk'))
    results = {}

    for key in pf.keys:
        if key == 'N Z A':
            nza = pf[key]
            results['N Z A_N_0'] = float(nza['N'][0])
            results['N Z A_Z_0'] = float(nza['Z'][0])
            results['N Z A_A_0'] = float(nza['A'][0])
            results['N Z A_len'] = len(nza['N'])
            continue
        d = np.asarray(pf._get_data(key))
        psin = np.asarray(pf.psinorm_for(key))
        results[f'{key}_len'] = len(d)
        results[f'{key}_0']   = float(d[0])
        results[f'{key}_mid'] = float(d[len(d)//2])
        results[f'{key}_end'] = float(d[-1])
        results[f'{key}_min'] = float(np.min(d))
        results[f'{key}_max'] = float(np.max(d))
        results[f'{key}_sum'] = float(np.sum(d))
        results[f'{key}_psinorm_end'] = float(psin[-1])
    return results


def test_eqdsk_pfile():
    os.chdir(test_dir)
    results = run_eqdsk_pfile_case()
    assert validate_dict([results], eqdsk_pfile_expected_dict)


# -----------------------------------------------------------------------
# Tests: COCOS conversion, sign flip, and round-trip serialisation
#
# ITER_test.eqdsk ships as COCOS 1.  `cocosify` applies the multiplicative
# sign / 2pi factors from Sauter & Medvedev, CPC 184 (2013) 293, Eq. 14/23
# to the raw g-file dict.  The sign table below mirrors `_cocos_params` in
# the module under test; the round-trip + per-field sign checks exercise
# every valid COCOS index (1-8 and 11-18).
# -----------------------------------------------------------------------
_COCOS_VALID = tuple(list(range(1, 9)) + list(range(11, 19)))


def _cocos_params_reference(cocos_index):
    """Reference copy of `_cocos_params` from eqdsk.py, kept locally so
    this test is an independent check of the table."""
    if cocos_index < 1 or cocos_index > 18 or cocos_index in (9, 10):
        raise ValueError(f'Invalid COCOS index: {cocos_index}')
    exp_Bp = 0 if cocos_index < 10 else 1
    base = cocos_index if cocos_index < 10 else cocos_index - 10
    return dict(
        sigma_Bp    = +1 if base in (1, 2, 5, 6) else -1,
        sigma_RpZ   = +1 if base in (1, 2, 7, 8) else -1,
        sigma_rhotp = +1 if base in (1, 3, 5, 7) else -1,
        exp_Bp      = exp_Bp,
    )


def _expected_cocos_factors(cocos_in, cocos_out):
    """Multiplicative factors cocosify() should apply when converting from
    cocos_in to cocos_out.  Returns a dict keyed by raw g-file field name."""
    cc_in  = _cocos_params_reference(cocos_in)
    cc_out = _cocos_params_reference(cocos_out)
    sBp    = cc_out['sigma_Bp']    * cc_in['sigma_Bp']
    sRpZ   = cc_out['sigma_RpZ']   * cc_in['sigma_RpZ']
    srhotp = cc_out['sigma_rhotp'] * cc_in['sigma_rhotp']
    exp_eff = cc_out['exp_Bp'] - cc_in['exp_Bp']
    twopi_exp = (2.0 * np.pi) ** exp_eff
    psi_fac  = sRpZ * sBp * twopi_exp
    dpsi_fac = sRpZ * sBp / twopi_exp
    bt_fac   = sRpZ
    ip_fac   = sRpZ
    q_fac    = srhotp
    return {
        'SIMAG':   psi_fac,
        'SIBRY':   psi_fac,
        'PSIRZ':   psi_fac,
        'PPRIME':  dpsi_fac,
        'FFPRIM':  dpsi_fac,
        'FPOL':    bt_fac,
        'BCENTR':  bt_fac,
        'CURRENT': ip_fac,
        'QPSI':    q_fac,
    }


@pytest.mark.parametrize('cocos_out', _COCOS_VALID)
def test_eqdsk_cocos_roundtrip(cocos_out):
    """cocosify(1 -> N -> 1) must reproduce the original raw g-file data."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    os.chdir(test_dir)
    eq_ref = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    eq_rt  = (eq_ref.cocosify(cocos_out, copy=True)
                    .cocosify(eq_ref.cocos, copy=True))
    assert eq_rt.cocos == eq_ref.cocos
    for key in ('SIMAG', 'SIBRY', 'BCENTR', 'CURRENT',
                'FPOL', 'PRES', 'PPRIME', 'FFPRIM', 'QPSI', 'PSIRZ'):
        assert np.allclose(eq_rt._raw[key], eq_ref._raw[key],
                           rtol=1e-12, atol=1e-12), f'{key} mismatch at COCOS {cocos_out}'


@pytest.mark.parametrize('cocos_out', _COCOS_VALID)
def test_eqdsk_cocos_signs(cocos_out):
    """After cocosify(1 -> N), each converted field equals the original
    multiplied by the sign/2pi factor predicted by the COCOS table."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    os.chdir(test_dir)
    eq_ref = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    eq_n   = eq_ref.cocosify(cocos_out, copy=True)
    assert eq_n.cocos == cocos_out
    factors = _expected_cocos_factors(eq_ref.cocos, cocos_out)
    for key, factor in factors.items():
        expected = np.asarray(eq_ref._raw[key]) * factor
        assert np.allclose(eq_n._raw[key], expected,
                           rtol=1e-12, atol=1e-12), (
            f'{key} expected scale {factor:+g} at COCOS {cocos_out}')


@pytest.mark.parametrize('bad', [0, 9, 10, 19, -1])
def test_eqdsk_cocos_invalid(bad):
    """Invalid COCOS indices must raise at construction or conversion."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    os.chdir(test_dir)
    with pytest.raises(ValueError):
        read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'), cocos=bad)
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    with pytest.raises(ValueError):
        eq.cocosify(bad, copy=True)


def test_eqdsk_flip_Bt_Ip_roundtrip():
    """flip_Bt_Ip twice is the identity on the raw g-file dict."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    os.chdir(test_dir)
    eq_ref = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    eq_rt  = eq_ref.flip_Bt_Ip(copy=True).flip_Bt_Ip(copy=True)
    for key in ('BCENTR', 'FPOL', 'CURRENT', 'SIMAG', 'SIBRY',
                'PSIRZ', 'PPRIME', 'FFPRIM'):
        assert np.allclose(eq_rt._raw[key], eq_ref._raw[key],
                           rtol=1e-12, atol=1e-12), f'{key} not invariant under flip^2'


def test_eqdsk_save_load_roundtrip(tmp_path):
    """Write ITER_test back to disk and re-read: raw scalars and 1-D/2-D
    arrays must round-trip within fixed-format ASCII precision."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    os.chdir(test_dir)
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    out_path = tmp_path / 'ITER_roundtrip.eqdsk'
    eq.save(str(out_path))
    eq_rt = read_geqdsk(str(out_path))
    for key in ('NW', 'NH', 'RDIM', 'ZDIM', 'RCENTR', 'RLEFT', 'ZMID',
                'RMAXIS', 'ZMAXIS', 'SIMAG', 'SIBRY', 'BCENTR', 'CURRENT'):
        assert np.isclose(eq._raw[key], eq_rt._raw[key], rtol=1e-6, atol=0), key
    for key in ('FPOL', 'PRES', 'PPRIME', 'FFPRIM', 'QPSI', 'PSIRZ',
                'RBBBS', 'ZBBBS', 'RLIM', 'ZLIM'):
        assert np.allclose(eq._raw[key], eq_rt._raw[key],
                           rtol=1e-6, atol=1e-10), key


def test_eqdsk_cocos_bytes_roundtrip():
    """to_bytes / from_bytes must reproduce the original equilibrium data."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk, GEQDSKEquilibrium
    os.chdir(test_dir)
    eq = read_geqdsk(os.path.join(test_dir, 'ITER_test.eqdsk'))
    raw_bytes = eq.to_bytes()
    eq_rt = GEQDSKEquilibrium.from_bytes(raw_bytes, cocos=eq.cocos)
    for key in ('FPOL', 'PRES', 'PSIRZ', 'QPSI'):
        assert np.allclose(eq._raw[key], eq_rt._raw[key],
                           rtol=1e-6, atol=1e-10), key


def test_pfile_bytes_roundtrip():
    """PFile.to_bytes / from_bytes must reproduce the original profiles."""
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_pfile, PFile
    os.chdir(test_dir)
    pf = read_pfile(os.path.join(test_dir, 'D3Dlike_Hmode_test.peqdsk'))
    pf_rt = PFile.from_bytes(pf.to_bytes())
    for key in ('ne', 'te', 'ni', 'ti'):
        d_ref = np.asarray(pf._get_data(key))
        d_rt  = np.asarray(pf_rt._get_data(key))
        assert np.allclose(d_ref, d_rt, rtol=1e-6, atol=1e-10), key


# -----------------------------------------------------------------------
# Test: X-point isoflux boundary helpers in TokaMaker.util
#
# These helpers check the generated boundary from 'create_isoflux_xpts'
# against the reference `isoflux_xpts_expected.json`.
# -----------------------------------------------------------------------
isoflux_xpts_expected = _load_eqdsk_fixture('isoflux_xpts_expected.json')


def test_xpoints_from_moments():
    """X-points are exact analytic Miller-moment coordinates."""
    r0, z0, a = 6.2, 0.3, 2.0
    kU, dU, kL, dL = 1.7, 0.33, 2.0, 0.5
    # Up-down symmetric (lower defaults to upper)
    xpts = xpoints_from_moments(r0, z0, a, kU, dU)
    assert xpts.shape == (2, 2)
    assert np.allclose(xpts[0], [r0 - a * dU, z0 + a * kU])
    assert np.allclose(xpts[1], [r0 - a * dU, z0 - a * kU])
    # Asymmetric
    xpts = xpoints_from_moments(r0, z0, a, kU, dU, kappa_lower=kL, delta_lower=dL)
    assert np.allclose(xpts[0], [r0 - a * dU, z0 + a * kU])
    assert np.allclose(xpts[1], [r0 - a * dL, z0 - a * kL])


@pytest.mark.parametrize('case', ['symmetric', 'asymmetric'])
def test_create_isoflux_xpts(case):
    """Generated boundary matches the reference contour, has exactly npts
    unique points, and passes through both X-points."""
    ref = isoflux_xpts_expected[case]
    params = ref['params']
    pts = create_isoflux_xpts(**params)
    # Shape: exactly npts points, no duplicate closing point
    assert pts.shape == (params['npts'], 2)
    closed = np.vstack([pts, pts[0]])
    gaps = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    assert np.all(gaps > 1e-9), "boundary contains a duplicate point"
    # Both X-points must appear on the contour
    xpts = np.asarray(ref['xpoints'])
    for xpt in xpts:
        assert np.min(np.linalg.norm(pts - xpt, axis=1)) < 1e-8
    # Regression against the stored reference contour
    assert np.allclose(pts, np.asarray(ref['points']), rtol=1e-10, atol=1e-12)


def test_create_isoflux_xpts_requires_min_npts():
    """npts < 3 must raise rather than allocate invalid segment counts."""
    with pytest.raises(ValueError):
        create_isoflux_xpts(2, 6.2, 0.0, 2.0, 1.7, 0.33)


# # Example of how to run single test without pytest
# if __name__ == '__main__':
#     multiprocessing.freeze_support()
#     mp_q = multiprocessing.Queue()
#     run_sph_case(0.05,2,mp_q)

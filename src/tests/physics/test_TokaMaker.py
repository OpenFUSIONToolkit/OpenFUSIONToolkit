from __future__ import print_function
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
from OpenFUSIONToolkit.TokaMaker import TokaMaker, gs_Domain, save_gs_mesh, load_gs_mesh


def mp_run(target,args,timeout=30):
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
    test_result = mp_q.get()
    p.join()
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
    mygs = TokaMaker()
    mygs.setup_mesh(mesh_pts,mesh_lc)
    mygs.settings.free_boundary = False
    mygs.setup(order=fe_order,F0=1.0,full_domain=True)
    mygs.pnorm=a
    mygs.alam=b*R*R*2.0
    mygs.set_profiles(ffp_prof={'type': 'flat'},pp_prof={'type': 'flat'})
    mygs.init_psi()
    psi_solovev_TM, _, rz_x = solovev_psi(mygs.r[:,0], mygs.r[:,1],R,a,b,c0)
    mygs.set_psi(-psi_solovev_TM)
    mygs.settings.nl_tol = 1.E-14
    mygs.update_settings()
    err_flag = mygs.solve()
    if err_flag != 0:
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


def validate_solo(results,psi_err_exp,X_err_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs((results[0]-psi_err_exp)/psi_err_exp) > 1.E-1:
        print("FAILED: psi error too high!")
        print("  Expected = {0}".format(psi_err_exp))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    if abs((results[1]-X_err_exp)/X_err_exp) > 1.E-1:
        print("FAILED: X-point error too high!")
        print("  Expected = {0}".format(X_err_exp))
        print("  Actual =   {0}".format(results[1]))
        test_result = False
    return test_result


# Test runners for Solov'ev cases
@pytest.mark.parametrize("order", (2,3,4))
def test_solo_h1(order):
    errs = [
        [3.2048631614233643e-07,0.00014929412629149645],
        [8.919954733021135e-10,4.659825491095631e-07],
        [5.084454462338564e-15,3.1842881888709607e-12]
    ]
    results = mp_run(run_solo_case,(0.015,order))
    assert validate_solo(results,errs[order-2][0],errs[order-2][1])
@pytest.mark.parametrize("order", (2,3,4))
def test_solo_h2(order):
    errs = [
        [7.725262474858205e-08,4.9243688140144384e-05],
        [1.1190059530634016e-10,2.919838380657025e-08],
        [1.0424769098635496e-14,3.434564147191569e-12]
    ]
    results = mp_run(run_solo_case,(0.015/2.0,order))
    assert validate_solo(results,errs[order-2][0],errs[order-2][1])
@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
def test_solo_h3(order):
    errs = [
        [2.0607919004158514e-08,5.955338556344096e-06],
        [1.3950375633902016e-11,1.154542061756696e-09],
        [2.0552832098467707e-14,6.859186868795993e-12]
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
    mygs = TokaMaker()
    mygs.setup_mesh(mesh_pts,mesh_lc)
    mygs.settings.free_boundary = False
    mygs.setup(order=fe_order)
    mygs.pnorm=0.0
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
    err_flag = mygs.solve()
    if err_flag != 0:
        mp_q.put(None)
        return
    psi_TM = mygs.get_psi(False)
    psi_eig_TM = spheromak_psi(mygs.r[:,0], mygs.r[:,1],1.0,1.0)
    psi_TM /= psi_TM.dot(psi_eig_TM)/psi_eig_TM.dot(psi_eig_TM)
    # Compute error in psi
    psi_err = np.linalg.norm(psi_TM-psi_eig_TM)/np.linalg.norm(psi_eig_TM)
    mp_q.put([psi_err])


def validate_sph(results,psi_err_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs((results[0]-psi_err_exp)/psi_err_exp) > 1.E-1:
        print("FAILED: psi error too high!")
        print("  Expected = {0}".format(psi_err_exp))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    return test_result


# Test runners for Spheromak cases
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
def run_coil_case(mesh_resolution,fe_order,mp_q):
    def coil_green(rc,zc,r,z,gs_obj):
        return gs_obj.eval_green(np.array([[r,z]]),np.array([rc,zc]))[0]
    def masked_err(point_mask,gs_obj,psi,sort_ind):
        bdry_points = gs_obj.r[point_mask,:]
        sort_ind = bdry_points[:,sort_ind].argsort()
        psi_bdry = psi[point_mask]
        psi_bdry = psi_bdry[sort_ind]
        bdry_points = bdry_points[sort_ind]
        green = np.zeros((bdry_points.shape[0],))
        for i in range(bdry_points.shape[0]):
            green[i], _ = dblquad(coil_green,0.75,0.85,0.75,0.85,args=(bdry_points[i,0],bdry_points[i,1],gs_obj))
        return green, psi_bdry
    # Build mesh
    gs_mesh = gs_Domain(rextent=1.0,zextents=[0.0,1.0])
    gs_mesh.define_region('air',mesh_resolution,'boundary')
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.define_region('coil',0.01,'coil')
    gs_mesh.add_rectangle(0.4,0.4,0.2,0.2,'plasma')
    gs_mesh.add_rectangle(0.8,0.8,0.1,0.1,'coil')
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    coil_dict = gs_mesh.get_coils()
    cond_dict = gs_mesh.get_conductors()
    # Run EQ
    mygs = TokaMaker()
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict)
    mygs.setup(order=fe_order)
    mygs.set_coil_currents(np.array([1.E-2]))
    err_flag = mygs.solve(True)
    if err_flag != 0:
        mp_q.put(None)
        return
    psi0 = mygs.get_psi(False)
    # Get analytic result
    green1, psi1 = masked_err(mygs.r[:,1]==1.0,mygs,psi0,0)
    green2, psi2 = masked_err(mygs.r[:,0]==1.0,mygs,psi0,1)
    green3, psi3 = masked_err(mygs.r[:,1]==-1.0,mygs,psi0,0)
    # Compute error in psi
    green_full = np.hstack((green1[1:], green2, green3[1:]))
    psi_full = np.hstack((psi1[1:], psi2, psi3[1:]))
    psi_err = np.linalg.norm(green_full+psi_full)/np.linalg.norm(green_full)
    mp_q.put([psi_err])


def validate_coil(results,psi_err_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    if abs((results[0]-psi_err_exp)/psi_err_exp) > 1.E-1:
        print("FAILED: psi error too high!")
        print("  Expected = {0}".format(psi_err_exp))
        print("  Actual =   {0}".format(results[0]))
        test_result = False
    return test_result


# Test runners for vacuum coil cases
@pytest.mark.parametrize("order", (2,3,4))
def test_coil_h1(order):
    errs = [0.010800921782063938, 0.0002851010669736233, 1.8185396736818836e-05]
    results = mp_run(run_coil_case,(0.1,order))
    assert validate_coil(results,errs[order-2])
@pytest.mark.parametrize("order", (2,3,4))
def test_coil_h2(order):
    errs = [0.0032993582771277, 2.725546769847347e-05, 8.670511127765199e-07]
    results = mp_run(run_coil_case,(0.1/2.0,order))
    assert validate_coil(results,errs[order-2])
@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
def test_coil_h3(order):
    errs = [0.0008175212508035045, 1.921137561342415e-06, 4.4282752350112954e-07]
    results = mp_run(run_coil_case,(0.1/4.0,order))
    assert validate_coil(results,errs[order-2])


#============================================================================
def run_ITER_case(mesh_resolution,fe_order,mp_q):
    def create_mesh():
        print(os.listdir())
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
            gs_mesh.define_region(key,coil_dx,'coil')
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
    mygs = TokaMaker()
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict)
    mygs.setup(order=fe_order,F0=5.3*6.2)
    vsc_signs = np.zeros((mygs.ncoils,), dtype=np.float64)
    vsc_signs[[coil_dict['VSU']['coil_id'], coil_dict['VSL']['coil_id']]] = [1.0,-1.0]
    mygs.set_coil_vsc(vsc_signs)
    coil_bounds = np.zeros((mygs.ncoils+1,2), dtype=np.float64)
    coil_bounds[:,0] = -50.E6; coil_bounds[:,1] = 50.E6
    mygs.set_coil_bounds(coil_bounds)
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
    coil_reg_mat = np.eye(mygs.ncoils+1, dtype=np.float64)
    coil_reg_weights = np.ones((mygs.ncoils+1,))
    coil_reg_targets = np.zeros((mygs.ncoils+1,))
    for key, coil in coil_dict.items():
        if key.startswith('CS'):
            if key.startswith('CS1'):
                coil_reg_weights[coil['coil_id']] = 2.E-2
            else:
                coil_reg_weights[coil['coil_id']] = 1.E-2
        elif key.startswith('PF'):
            coil_reg_weights[coil['coil_id']] = 1.E-2
        elif key.startswith('VS'):
            coil_reg_weights[coil['coil_id']] = 1.E2
    coil_reg_weights[-1] = 1.E-2
    mygs.set_coil_reg(coil_reg_mat, reg_weights=coil_reg_weights, reg_targets=coil_reg_targets)
    n_sample = 40
    psi_sample = np.linspace(0.0,1.0,n_sample)
    alpha = 1.5
    gamma = 2.0
    ffp_prof = {
        'type': 'linterp',
        'x': psi_sample,
        'y': np.power(1.0-np.power(psi_sample,alpha),gamma)
    }
    ffp_prof['y'] /= ffp_prof['y'][0]
    alpha = 4.0
    gamma = 1.0
    pp_prof = {
        'type': 'linterp',
        'x': psi_sample,
        'y': np.power(1.0-np.power(psi_sample,alpha),gamma)
    }
    pp_prof['y'] /= pp_prof['y'][0]
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
    R0 = 6.3
    Z0 = 0.5
    a = 2.0
    kappa = 1.4
    delta = 0.0
    err_flag = mygs.init_psi(R0, Z0, a, kappa, delta)
    err_flag = mygs.solve()
    if err_flag != 0:
        mp_q.put(None)
        return
    eq_info = mygs.get_stats()
    mp_q.put([eq_info])


def validate_ITER(results,dict_exp):
    if results is None:
        print("FAILED: error in solve!")
        return False
    test_result = True
    for key, exp_val in dict_exp.items():
        result_val = results[0].get(key,None)
        if result_val is None:
            print('FAILED: key "{0}" not present!'.format(key))
        else:
            if type(exp_val) is not list:
                if abs((result_val-exp_val)/exp_val) > 1.E-2:
                    print("FAILED: {0} error too high!".format(key))
                    print("  Expected = {0}".format(exp_val))
                    print("  Actual =   {0}".format(result_val))
                    test_result = False
    return test_result


# Test runners for ITER test cases
@pytest.mark.parametrize("order", (2,3))#3,4))
def test_ITER(order):
    exp_dict = {
        'Ip': 15599996.700479196,
        'Ip_centroid': [6.20274133, 0.5296048],
        'kappa': 1.86799695311941,
        'kappaU': 1.7388335731481432,
        'kappaL': 1.997160333090677,
        'delta': 0.4642130933423834,
        'deltaU': 0.3840631923067706,
        'deltaL': 0.5443629943779958,
        'vol': 820.0973897169655,
        'q_0': 0.8234473499435633,
        'q_95': 2.76048354704068,
        'P_ax': 619225.0167519478,
        'W_MHD': 242986888.67690986,
        'beta_pol': 39.73860565406112,
        'dflux': 1.5402746036620532,
        'tflux': 121.86870301036512,
        'l_i': 0.9048845463517069,
        'beta_tor': 1.768879437469196
    }
    results = mp_run(run_ITER_case,(1.0,order))
    assert validate_ITER(results,exp_dict)
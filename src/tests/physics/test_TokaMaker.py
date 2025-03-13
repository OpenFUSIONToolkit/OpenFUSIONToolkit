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
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh, load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_isoflux, eval_green, create_power_flux_fun
from OpenFUSIONToolkit.util import oftpy_dump_cov


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


def validate_dict(results,dict_exp):
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
                    if abs((result_val[i]-exp_val[i])/exp_val[i]) > 1.E-2:
                        print("FAILED: {0} ({1}) error too high!".format(key,i))
                        print("  Expected = {0}".format(exp_val[i]))
                        print("  Actual =   {0}".format(result_val[i]))
                        test_result = False
            else:
                if abs((result_val-exp_val)/exp_val) > 1.E-2:
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
    mygs.pnorm=a
    mygs.alam=b*R*R*2.0
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
    def coil_green(rc,zc,r,z):
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
    # Build mesh
    gs_mesh = gs_Domain(rextent=1.0,zextents=[0.0,1.0])
    gs_mesh.define_region('air',mesh_resolution,'boundary')
    gs_mesh.define_region('plasma',mesh_resolution,'plasma')
    gs_mesh.define_region('coil',0.01,'coil')
    gs_mesh.add_rectangle(px1,py1,pdx,pdy,'plasma')
    gs_mesh.add_rectangle(cx1,cy1,cdx,cdy,'coil')
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    coil_dict = gs_mesh.get_coils()
    cond_dict = gs_mesh.get_conductors()
    # Run EQ
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order=fe_order)
    mygs.set_coil_currents({'COIL': cdx*cdy})
    if dist is not None:
        mygs.set_coil_current_dist('COIL',dist(mygs.r[:,0],mygs.r[:,1]))
    try:
        psi0 = mygs.vac_solve()
    except ValueError:
        mp_q.put(None)
        return

    # Get analytic result
    green1, psi1 = masked_err(mygs.r[:,1]==1.0,mygs,psi0,0)
    green2, psi2 = masked_err(mygs.r[:,0]==1.0,mygs,psi0,1)
    green3, psi3 = masked_err(mygs.r[:,1]==0.0,mygs,psi0,0)
    # Compute error in psi
    green_full = np.hstack((green1[1:], green2, green3[1:]))
    psi_full = np.hstack((psi1[1:], psi2, psi3[1:]))
    psi_err = np.linalg.norm(green_full+psi_full)/np.linalg.norm(green_full)
    mp_q.put([psi_err])
    oftpy_dump_cov()


def validate_coil(results,psi_err_exp):
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


# Test runners for vacuum coil cases
def coil_dist(r,z):
    return r-z

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h1(order,dist_coil):
    if dist_coil:
        errs = np.r_[0.01840042334178343, 0.003450061648683903, 0.0008471927795560409]
        results = mp_run(run_coil_case,(0.1,order,coil_dist))
    else:
        errs = np.r_[0.010702389576304984, 0.00028240755964702516, 1.7930442330763583e-05]
        results = mp_run(run_coil_case,(0.1,order,None))
    print('Err = ',results[0])
    assert validate_coil(results,errs[order-2])
@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h2(order,dist_coil):
    if dist_coil:
        errs = np.r_[0.0036698088466649878, 0.00021755100020083888, 1.736479174843997e-05]
        results = mp_run(run_coil_case,(0.1/2.0,order,coil_dist))
    else:
        errs = np.r_[0.0032680822197860876, 2.7032824967342426e-05, 8.560758830069598e-07]
        results = mp_run(run_coil_case,(0.1/2.0,order,None))
    assert validate_coil(results,errs[order-2])
@pytest.mark.slow
@pytest.mark.parametrize("order", (2,3,4))
@pytest.mark.parametrize("dist_coil", (False, True))
def test_coil_h3(order,dist_coil):
    if dist_coil:
        errs = np.r_[0.001364423661608862, 1.5953386454257285e-05, 9.158565258919996e-07]
        results = mp_run(run_coil_case,(0.1/4.0,order,coil_dist))
    else:
        errs = np.r_[0.0008094155097004184, 1.8949323808351823e-06, 4.4169705023586007e-07]
        results = mp_run(run_coil_case,(0.1/4.0,order,None))
    assert validate_coil(results,errs[order-2])


#============================================================================
def run_ITER_case(mesh_resolution,fe_orders,eig_test,stability_test,mp_q):
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
        mygs.setup(order=fe_order,F0=5.3*6.2)
        #
        if eig_test:
            eig_vals, _ = mygs.eig_wall(10)
            mp_q.put([{'Tau_w': eig_vals[:5,0]}])
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
            mygs.solve()
        except ValueError:
            mp_q.put(None)
            return
        if stability_test:
            eig_vals, _ = mygs.eig_td(-1.E2,10,False)
            mp_q.put([{'gamma': eig_vals[:5,0]}])
            return
        mygs.save_eqdsk('test.eqdsk',lcfs_pressure=6.E4)
        eq_info = mygs.get_stats(li_normalization='ITER')
        Lmat = mygs.get_coil_Lmat()
        eq_info['LCS1'] = Lmat[mygs.coil_sets['CS1U']['id'],mygs.coil_sets['CS1U']['id']]
        eq_info['MCS1_plasma'] = Lmat[mygs.coil_sets['CS1U']['id'],-1]
        eq_info['Lplasma'] = Lmat[-1,-1]
    # Test deletion if multiple cases
    if mygs_last is not None:
        del mygs_last
    # Save final one
    mp_q.put([eq_info])
    oftpy_dump_cov()


# Test runners for ITER test cases
@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_eig(order):
    exp_dict = {
        'Tau_w': [1.51083009, 2.87431718, 3.91493237, 5.23482507, 5.61049374]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),True,False))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_stability(order):
    exp_dict = {
        'gamma': [-12.3620, 1.83981, 3.41613, 5.12470, 6.53393]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),False,True))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_eq(order):
    exp_dict = {
        'Ip': 15599996.700479196,
        'Ip_centroid': [6.20274133, 0.5296048],
        'kappa': 1.86799695311941,
        'kappaU': 1.7388335731481432,
        'kappaL': 1.997160333090677,
        'delta': 0.4642130933423834,
        # 'deltaU': 0.3840631923067706, # Disable for now
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
        'beta_tor': 1.7816206668692283,
        'beta_n': 1.1868590722509704,
        'LCS1': 2.485860941880887e-06,
        'MCS1_plasma': 8.930926092661585e-07,
        'Lplasma': 1.1899835061690724e-05
    }
    results = mp_run(run_ITER_case,(1.0,(order,),False,False))
    assert validate_dict(results,exp_dict)

def test_ITER_concurrent():
    exp_dict = {
        'Ip': 15599996.700479196,
        'Ip_centroid': [6.20274133, 0.5296048],
        'kappa': 1.86799695311941,
        'kappaU': 1.7388335731481432,
        'kappaL': 1.997160333090677,
        'delta': 0.4642130933423834,
        # 'deltaU': 0.3840631923067706, # Disable for now
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
        'beta_tor': 1.7816206668692283,
        'beta_n': 1.1868590722509704,
        'LCS1': 2.485860941880887e-06,
        'MCS1_plasma': 8.930926092661585e-07,
        'Lplasma': 1.1899835061690724e-05
    }
    results = mp_run(run_ITER_case,(1.0,(2,3),False,False))
    assert validate_dict(results,exp_dict)

#============================================================================
def run_LTX_case(fe_order,eig_test,stability_test,mp_q):
    def create_mesh():
        with open('LTX_geom.json','r') as fid:
            LTX_geom = json.load(fid)
        plasma_dx = 0.02
        coil_dx = 0.02
        vac_dx = 0.05
        gs_mesh = gs_Domain()
        #
        gs_mesh.define_region('air',vac_dx,'boundary')
        gs_mesh.define_region('plasma',plasma_dx,'plasma')
        gs_mesh.define_region('shellU',2*plasma_dx,'conductor',eta=4.E-7,noncontinuous=True)
        gs_mesh.define_region('shellL',2*plasma_dx,'conductor',eta=4.E-7,noncontinuous=True)
        for i, vv_segment in enumerate(LTX_geom['vv']):
            gs_mesh.define_region('vv{0}'.format(i),2*plasma_dx,'conductor',eta=vv_segment[1])
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
    if eig_test:
        eig_vals, _ = mygs.eig_wall(10)
        mp_q.put([{'Tau_w': eig_vals[:5,0]}])
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
    if stability_test:
        eig_vals, _ = mygs.eig_td(-1.E3,10,False)
        mp_q.put([{'gamma': eig_vals[:5,0]}])
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
        'Tau_w': [195.300148, 253.92961287, 403.74576838, 473.64151856, 550.08441557]
    }
    results = mp_run(run_LTX_case,(order,True,False))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_stability(order):
    exp_dict = {
        'gamma': [-238.0708, 216.6903, 286.1825, 394.7443, 394.7443]
    }
    results = mp_run(run_LTX_case,(order,False,True))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_LTX_eq(order):
    exp_dict = {
        'Ip': 90000.1298205169,
        'Ip_centroid': [4.05471907e-01, None],
        'kappa': 1.5213293087744595,
        'kappaU': 1.5215960005535605,
        'kappaL': 1.5210626169953587,
        # 'delta': 0.12295683642943968, # Disable for now
        # 'deltaU': 0.12289529426354395, # Disable for now
        # 'deltaL': 0.12301837859533517, # Disable for now
        'vol': 0.6511641559778095,
        'q_0': 1.3276880560032807,
        'q_95': 5.897372049554493,
        'P_ax': 1721.5000219840285,
        'W_MHD': 563.3140773902292,
        'beta_pol': 40.358078684354965,
        'dflux': 0.0009602066573419095,
        'tflux': 0.08571976036037239,
        'l_i': 1.002735427186787,
        'beta_tor': 1.9398553532544882,
        'beta_n': 1.38790732317241
    }
    results = mp_run(run_LTX_case,(order,False,False))
    assert validate_dict(results,exp_dict)

# Example of how to run single test without pytest
# if __name__ == '__main__':
#     multiprocessing.freeze_support()
#     mp_q = multiprocessing.Queue()
#     run_sph_case(0.05,2,mp_q)

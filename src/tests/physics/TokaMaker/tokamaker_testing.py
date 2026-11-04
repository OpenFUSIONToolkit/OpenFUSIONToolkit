#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Shared helpers for TokaMaker tests: process runner, result validation, and the
# ITER mesh, constraints, and profiles used by the ITER-based cases
import os
import sys
import json
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain, save_gs_mesh, load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun


def mp_run(target,args,timeout=30):
    return oft_testing.mp_run(target,args,test_dir,timeout=timeout)


def load_fixture(filename):
    with open(os.path.join(test_dir, filename)) as fh:
        return json.load(fh)


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
# ITER configuration shared by the equilibrium, reconstruction, bootstrap,
# and toroidal flux cases
#============================================================================
# LCFS isoflux points and lower X-point
ITER_isoflux_pts = np.array([
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
ITER_x_point = np.array([[5.125, -3.4],])

# Reference values for the ITER equilibrium from `solve_ITER_eq`
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


def create_ITER_mesh(mesh_resolution=1.0):
    '''Build the ITER mesh and save it to `ITER_mesh.h5`, unless that file already exists'''
    if os.path.exists('ITER_mesh.h5'):
        return
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


def setup_ITER_gs(myOFT,fe_order,F0=5.3*6.2,maxits=None):
    '''Create a TokaMaker object on the ITER mesh (see `create_ITER_mesh`)'''
    mygs = TokaMaker(myOFT)
    mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts,mesh_lc,mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    if maxits is not None:
        mygs.settings.maxits = maxits
    mygs.setup(order=fe_order,F0=F0)
    return mygs


def set_ITER_constraints(mygs,Ip_target=15.6E6,P0_target=6.2E5):
    '''Set the vertical stability coil, coil bounds, global targets, LCFS shape, and coil regularization'''
    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.E6, 50.E6] for key in mygs.coil_sets})
    mygs.set_targets(Ip=Ip_target, pax=P0_target)
    mygs.set_isoflux(np.vstack((ITER_isoflux_pts,ITER_x_point)))
    mygs.set_saddles(ITER_x_point)
    # Set regularization weights (all ITER coil sets are CS, PF or VS coils)
    regularization_terms = []
    for name in mygs.coil_sets:
        if name.startswith('CS1'):
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=2.E-2))
        else:
            regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=0.0,weight=1.E-2))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E2))
    mygs.set_coil_reg(reg_terms=regularization_terms)


def solve_ITER_eq(mygs):
    '''Initialize and solve the reference ITER equilibrium, returning the equilibrium and iteration count'''
    ffp_prof = create_power_flux_fun(40,1.5,2.0)
    pp_prof = create_power_flux_fun(40,4.0,1.0)
    mygs.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)
    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
    return mygs.solve(return_its=True)


def ITER_Hmode_profiles(n_sample):
    '''H-mode electron density [m^-3] and temperature [eV] profiles on `n_sample` uniform psi_N points'''
    from OpenFUSIONToolkit.TokaMaker.bootstrap import Hmode_profiles
    xphalf = 0.965
    ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n_sample,
                        expin=1.6, expout=1.6, widthp=0.35, xphalf=xphalf) * 1e20
    Te = Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=n_sample,
                        expin=1.3, expout=1.7, widthp=0.1, xphalf=xphalf)
    return ne, Te


def setup_ITER_bootstrap_internal(mesh_resolution, fe_order):
    '''ITER mesh, H-mode profiles and an initial no-bootstrap solve for the internal bootstrap cases; None on failure'''
    try:
        create_ITER_mesh(mesh_resolution)
    except Exception as e:
        print(e)
        return None

    # --- Kinetic and current profiles (ITER H-mode) ---
    n_sample = 257
    psi_sample = np.linspace(0.0, 1.0, n_sample)
    Ip_target = 13.0e6
    Zeff_val = 1.5
    ne, Te = ITER_Hmode_profiles(n_sample)
    ni = ne.copy()
    Ti = Te.copy()
    inductive_jphi = create_power_flux_fun(n_sample, 2.25, 2.5)['y']

    # --- Set up GS solver ---
    myOFT = OFT_env(nthreads=-1)
    mygs = setup_ITER_gs(myOFT, fe_order, maxits=100)
    set_ITER_constraints(mygs, Ip_target=Ip_target)

    # Initial no-bootstrap solve
    mygs.settings.pm = False
    mygs.update_settings()
    try:
        mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
        mygs.solve()
    except ValueError:
        return None
    return dict(mygs=mygs, myOFT=myOFT, psi_sample=psi_sample, n_sample=n_sample, Ip_target=Ip_target,
                Zeff_val=Zeff_val, ne=ne, Te=Te, ni=ni, Ti=Ti, inductive_jphi=inductive_jphi)

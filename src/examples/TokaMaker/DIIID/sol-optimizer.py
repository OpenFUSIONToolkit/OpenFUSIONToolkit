from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, read_mhdin, read_kfile
import numpy as np
# import matplotlib.pyplot as plt
import sys
import traceback
import os
import json
from collections import OrderedDict

def convert_Mirnov_coordinates(sensor, angle):
    R, z, _ = sensor
    position = np.array([R, z])

    angle_rads = angle * np.pi / 180
    orientation = np.array([np.cos(angle_rads), 0, np.sin(angle_rads)])
    return position, orientation

def read_points_flux(file_path,  sensor_names = None):
    points = OrderedDict()
    with open(file_path, 'r') as file:
        for line in file:
            if sensor_names is not None:
                name, r, z = line.strip().split()  # Split and ignore the first field (name)
                if name in sensor_names:
                    points[name] = [float(r),float(z)]
            else:
                name, r, z = line.strip().split()  # Split and ignore the first field (name)
                points[name] = [float(r),float(z)]
    return points

def read_points_mag(file_path, sensor_names = None):
    points = OrderedDict()
    with open(file_path, 'r') as file:
        for line in file:
            if sensor_names is not None: 
                name, r, z, pol, orient, _, _, _ = line.strip().split()  # Split and ignore the first field (name)
                if name in sensor_names:
                    points[name] = [float(r), float(z), float(orient)]
            else: 
                name, r, z, pol, orient, _, _, _ = line.strip().split()  # Split and ignore the first field (name)
                points[name] = [float(r), float(z), float(orient)]  # Convert r and z to floats and append
    return points

def load_shot(mygs, shot_tag):
    eqdsk = read_eqdsk('DIIID_files/g'+shot_tag)
    e_coil_names = ['ECOILA','ECOILB','E567UP','E567DN','E89DN','E89UP']
    f_coil_names = ['F1A', 'F2A', 'F3A', 'F4A', 'F5A', 'F6A', 'F7A', 'F8A', 'F9A', 'F1B', 'F2B', 'F3B', 'F4B', 'F5B', 'F6B', 'F7B', 'F8B', 'F9B']
    machine_dict = read_mhdin('DIIID_files/mhdin.dat', e_coil_names, f_coil_names)
    probes_dict, loops_dict, e_coil_dict, f_coil_dict = read_kfile('DIIID_files/k'+shot_tag, e_coil_names, f_coil_names, machine_dict)

    plasma_dx = 0.04
    coil_dx = 0.03
    vv_dx = 0.04
    vac_dx = 0.10

    with open('DIIID_files/DIIID_geom.json','r') as fid:
        DIIID_geom = json.load(fid)

    # Create a G-S domain
    gs_mesh = gs_Domain()
    # Define region information for mesh
    gs_mesh.define_region('air',vac_dx,'boundary')                     # Define the bounding region
    gs_mesh.define_region('plasma',plasma_dx,'plasma')                 # Define the plasma region and resolution
    gs_mesh.define_region('vacuum',vv_dx,'vacuum',allow_xpoints=True)  # Define the vacuum inside the VV
    # Define regions for VV
    for i, vv_segment in enumerate(DIIID_geom["vv"]):
        gs_mesh.define_region('vv{0}'.format(i),vv_dx,'conductor',eta=vv_segment[1])

    # Define geometry
    gs_mesh.add_polygon(DIIID_geom['limiter'],'plasma',parent_name='vacuum')  # Define the shape of the limiter
    gs_mesh.add_enclosed([1.75,1.25],'vacuum')
    # Define regions for VV
    for i, vv_segment in enumerate(DIIID_geom["vv"]):
        gs_mesh.add_polygon(vv_segment[0],'vv{0}'.format(i),parent_name='air')

    # Define F coils
    for key, coil in DIIID_geom["coils"].items():
        if key.startswith('ECOIL'):
            continue
        gs_mesh.define_region(key,coil_dx,'coil',nTurns=machine_dict['FCOIL'][key][4])
        gs_mesh.add_polygon(coil["pts"],key,parent_name='air')

    # Define E Coils
    for coil_set_name, coil_set in machine_dict['ECOIL'].items():
        i = 0
        for coil in coil_set:
            coil_name = "{}_{}".format(coil_set_name, i)
            i = i + 1
            r, z, w, h = coil
            gs_mesh.define_region(coil_name,coil_dx,'coil',coil_set=coil_set_name,nTurns=1) # Check nTurns
            gs_mesh.add_rectangle(r, z, w, h, coil_name, parent_name='air')
    
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    coil_dict = gs_mesh.get_coils()
    cond_dict = gs_mesh.get_conductors()

    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
    mygs.setup(order = 2, F0 = eqdsk['rcentr']*eqdsk['bcentr'])

    mygs.update_settings()

    vsc_signs = {key: 0 for key in mygs.coil_sets}
    vsc_signs['F9A'] = 1.0
    vsc_signs['F9B'] = -1.0
    mygs.set_coil_vsc(vsc_signs)
    return probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk

def load_mirnov(myrecon):
    magSensors322 = read_points_mag('DIIID_files/3222magsensorloc.txt')
    B_vals = [probes_dict[key][0] for key in magSensors322 if key in probes_dict]
    B_mean = np.mean(np.positive(B_vals))
    mirnov_names = []

    for key, mag in magSensors322.items():
        mirnov_names.append(key)
        position, orientation = convert_Mirnov_coordinates(mag, machine_dict['AMP2'][key])
        if key not in probes_dict:
            print("No data for probe {}".format(key))
            continue
        B_meas, _ = probes_dict[key]
        myrecon.add_Mirnov(position, orientation, B_meas, err=0.1*max(B_mean, abs(B_meas)))
    
    return magSensors322

def load_flux_loop(mygs, myrecon):    
    fluxLoops_onf = read_points_flux('DIIID_files/fcoil_fluxloops.txt')
    fluxLoops_vv = read_points_flux('DIIID_files/vv_fluxloops.txt')
    flux_locs = []
    flux_vals = []

    psi_vals_1 = [loops_dict[key][0] for key in fluxLoops_onf if key in loops_dict]
    psi_mean_1 = np.mean(np.positive(psi_vals_1))
    psi_vals_2 = [loops_dict[key][0] for key in fluxLoops_vv if key in loops_dict]
    psi_mean_2 = np.mean(np.positive(psi_vals_2))

    flux_names = []

    for key, fl in fluxLoops_onf.items():
        if key not in loops_dict:
            print("No data for Loop {}".format(key))
            continue
        flux_names.append(key)
        B_tmp, _ = loops_dict[key]
        flux_locs.append(fl)
        flux_vals.append(B_tmp)
        psi_val = B_tmp*2.0*np.pi
        myrecon.add_flux_loop(fl, psi_val, err=0.1*(max(psi_mean_1, abs(psi_val))))

    for key, fl in fluxLoops_vv.items():
        if key not in loops_dict:
            print("No data for Loop {}".format(key))
            continue
        flux_names.append(key)
        B_tmp, _ = loops_dict[key]
        flux_locs.append(fl)
        flux_vals.append(B_tmp)
        psi_val = B_tmp*2.0*np.pi
        myrecon.add_flux_loop(fl, psi_val, err=0.1*(max(psi_mean_2, abs(psi_val))))
    
    mygs.set_flux(np.array(flux_locs), np.array(flux_vals))

    loop_coords = {**fluxLoops_onf, **fluxLoops_vv}
    return loop_coords, np.array(flux_locs), np.array(flux_vals)

def load_coil(mygs):
    target_currents = {}
    target_errs = {}

    for key, value in e_coil_dict.items():
        target_currents[key] = value[0]
        target_errs[key] = abs(value[1])

    for key, value in f_coil_dict.items():
        target_currents[key] = f_coil_dict[key][0] / machine_dict['FCOIL'][key][4]
        target_errs[key] = abs(f_coil_dict[key][1] / machine_dict['FCOIL'][key][4])

    # Set coil regularization to weakly track measured coil currents
    regularization_terms = []
    for key in e_coil_dict:
        regularization_terms.append(mygs.coil_reg_term({key: 1.0},
                                                    target=target_currents[key],
                                                    weight=1.0E1))
    for key in f_coil_dict:
        weight = 1.0
        if key.startswith('F5'):
            weight = 1.0E2
        regularization_terms.append(mygs.coil_reg_term({key: 1.0},
                                                    target=target_currents[key],
                                                    weight=1.0E2))

    # Set zero target current and small weight on virtual VSC to allow up-down adjustment
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-2))

    # Pass regularization terms to TokaMaker
    mygs.set_coil_reg(reg_terms=regularization_terms)    

def eq_reconstruct(mygs, myrecon, ffp_prof, probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk, probe_coords, loop_coords):
    Ip_target=abs(eqdsk['ip'])
    P0_target=eqdsk['pres'][0]
    mygs.set_targets(Ip=Ip_target, pax=P0_target)

    # Set reconstruction settings
    myrecon.settings.fitPnorm = False
    myrecon.settings.fitR0 = True
    myrecon.settings.fitV0 = True
    myrecon.settings.fitCoils = True
    myrecon.settings.pm = False

    pprime = [1.0, 0.0, 0.0]
    psi_pprime = [0.0, 1.0, 1.5]
    psi_sample = np.linspace(0.0,1.5,150)
    psi_prof = psi_sample.copy()
    pp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_pprime,pprime)))).copy()
    mygs.set_profiles(pp_prof={'type': 'linterp', 'y': pp_prof[:,1], 'x': psi_sample}, f_SOL=True)
    mygs.set_profiles(ffp_prof={'type': 'linterp', 'y': ffp_prof[:,1], 'x': psi_sample},pp_prof={'type': 'linterp', 'y': pp_prof[:,1], 'x': psi_sample},f_SOL=True)

    # Initial equilibrium with very rough guess
    R0 = eqdsk['rcentr']
    Z0 = 0.0
    a = 0.6
    kappa = 1.13
    delta = .1

    err_flag = mygs.init_psi(R0, Z0, a, kappa, delta)
    if err_flag:
        print("Initialzie PSI Failed.")
        sys.exit(1)

    print("Initial Solve...")
    mygs.solve()

    # Remove all shape constraints
    mygs.set_isoflux(None)
    mygs.set_flux(None,None)
    mygs.set_saddles(None)

    # Set Recon Plasma Current
    myrecon.set_Ip(Ip_target, err=0.1*Ip_target)

    # Set initial position targets from current values
    mygs.set_targets(R0=mygs.o_point[0],V0=mygs.o_point[1])

    # Perform reconstructions
    print("Reconstructing...")
    err_flag = myrecon.reconstruct()

    if err_flag:
        raise(Exception("Reconstruction Failed."))
    
    # psi_recons = mygs.get_psi()
    # fig, ax = plt.subplots(1,1)
    # mygs.plot_machine(fig,ax,coil_colormap='seismic',coil_scale=1.E-6,coil_clabel=r'$I_C$ [MA]',coil_symmap=True)
    # mygs.plot_psi(fig,ax,plasma_nlevels=5,vacuum_nlevels=5)
    # mygs.plot_psi(fig,ax,psi_recons,plasma_levels=[1.0,],plasma_color='red',vacuum_nlevels=0,plasma_linestyles='dashed')
    # plt.savefig('out/LCFS-{}-i={}.png'.format(shot_tag, i))

    coil_currents, _ = mygs.get_coil_currents()
    coil_dict = {**e_coil_dict, **f_coil_dict}
    current_err = 0.0
    for key in coil_dict:
        coil_uncty = coil_dict[key][1]
        if key.startswith('F'):
            coil_uncty = np.abs(f_coil_dict[key][1] / machine_dict['FCOIL'][key][4])
        sq_diff = ((coil_currents[key] - coil_dict[key][0]) * coil_uncty) ** 2 
        current_err += sq_diff

    B_eval = mygs.get_field_eval('B')
    probe_err = 0.0
    for key in probes_dict:
        if key not in probe_coords:
            print("Probe {} not recognized. Skipping analysis.".format(key))
            continue
        coords = probe_coords[key]
        position, orientation = convert_Mirnov_coordinates(coords, machine_dict['AMP2'][key])
        B_field = np.dot(B_eval.eval(position), orientation)
        sq_diff = ((B_field - probes_dict[key][0]) * probes_dict[key][1]) ** 2
        probe_err += sq_diff

    psi_eval = mygs.get_field_eval('PSI')
    loop_err = 0.0
    for key in loops_dict:
        if key not in loop_coords:
            print("Loop {} not recognized. Skipping analysis.".format(key))
            continue
        position = loop_coords[key]
        psi_realized = psi_eval.eval(position)[0] * 2 * np.pi
        sq_diff = ((psi_realized - loops_dict[key][0]) * loops_dict[key][1]) ** 2
        loop_err += sq_diff

    return current_err, probe_err, loop_err

def gen_linear_profile(params):
    p1, p2 = params
    ffprime = [1.0, p1, p2, 0.0]
    psi_def = [0.0, 0.5, 1.0, 1.5]
    psi_sample = np.linspace(0.0, 1.5, 150)
    psi_prof = psi_sample.copy()
    ffp_prof = np.transpose(np.vstack((psi_prof, np.interp(psi_sample, psi_def, ffprime))))
    return ffp_prof.copy()

def gen_bump_profile(params):
    delta, psi_sol = params
    def simple_bump(psi_n):
        return np.exp(-1.0 * delta ** 2 / (psi_sol**2 - psi_n**2))
    x = np.linspace(0.0, 1.5, 150)
    y = [simple_bump(x_i) for x_i in x]
    ffp_prof = np.array([[x[i], y[i]] for i in range(len(x))])
    return ffp_prof

def total_err(curr_err, probe_err, loop_err):
    return curr_err + probe_err + loop_err

# Initialize and Perform Optimization
if len(sys.argv) < 2:
    print("USAGE: python sol-optimizer.py <SHOT_TAG> <LINEAR/BUMP>")
    sys.exit(0)
shot_tag = sys.argv[1]

use_bump = False
if len(sys.argv) > 2:
    mode = sys.argv[2].upper()
    if mode not in ['LINEAR', 'BUMP']:
        print("USAGE: python sol-optimizer.py <SHOT_TAG> <LINEAR/BUMP>")
        sys.exit(0)
    use_bump = mode == 'BUMP'
    print("Using {} mode.".format(mode))
else:
    print("Mode not specified. Defaulting to linear.")

print("=== LOADING SHOT FILE {} ===".format(shot_tag))
myOFT = OFT_env(nthreads=2)
mygs = TokaMaker(myOFT)
mygs.settings.maxits=1000
myrecon = reconstruction(mygs)
probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk = load_shot(mygs, shot_tag)
mygs.update_settings()

probe_coords = load_mirnov(myrecon)
loop_coords, flux_locs, flux_vals = load_flux_loop(mygs, myrecon)
load_coil(mygs)

print("=== BEGINNING OPTIMIZATION ===")
paramter_range = np.linspace(0.0, 1.0, 4)
parameter_space = [(p1, p2) for p1 in paramter_range for p2 in paramter_range]
if use_bump:
    delta_space= np.linspace(0.5, 3.0, 5)
    psi_sol_space = np.linspace(1.5, 2.0, 5)
    parameter_space = [(delta, psi_sol) for delta in delta_space for psi_sol in psi_sol_space]

min_err = float('inf')
min_params = None
i = 0
for params in parameter_space:
    try:
        ffp_prof = gen_linear_profile(params)
        if use_bump:
            ffp_prof = gen_bump_profile(params)
        print("RUNNING OPTIIMIZATION [{}/{}]".format(i + 1, len(parameter_space)))
        # Reset Flux Loops
        mygs.set_flux(flux_locs, flux_vals)
        # Run Reconstruction
        curr_err, probe_err, loop_err = eq_reconstruct(mygs, myrecon, ffp_prof, probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk, probe_coords, loop_coords)
        # Check Error
        err = total_err(curr_err, probe_err, loop_err)
        if err < min_err:
            min_err = err
            min_params = params
        print("Last Err = {}".format(err))
        print("Last Parameters = {}".format(params))
        print("Min Err = {}".format(min_err))
        print("Min Parameters = {}".format(min_params))
        print()
    except Exception as e:
        print("Solve failed for step={} parameters={}".format(i + 1, params))
        traceback.print_exc()
    finally:
        i = i + 1

print("=== OPTIMIZATION COMPLETE ===")
print("Min Err = {}".format(min_err))
print("Min Parameters = {}".format(min_params))
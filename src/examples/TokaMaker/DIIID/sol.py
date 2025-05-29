from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, read_mhdin, read_kfile
import numpy as np
# import matplotlib.pyplot as plt
import sys
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
    machine_dict, _ = read_mhdin('DIIID_files/mhdin.dat', e_coil_names, f_coil_names)
    probes_dict, loops_dict, e_coil_dict, f_coil_dict, _ = read_kfile('DIIID_files/k'+shot_tag, machine_dict, e_coil_names, f_coil_names)

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

def load_mirnov(myrecon, probes_dict, machine_dict):
    magSensors322 = read_points_mag('DIIID_files/3222magsensorloc.txt')
    B_vals = [probes_dict[key][0] for key in magSensors322 if key in probes_dict]
    B_mean = np.mean(np.positive(B_vals))
    mirnov_names = []

    for key, mag in magSensors322.items():
        if probes_dict[key][1] == 0:
            print("No data for probe {}".format(key))
            continue
        position, orientation = convert_Mirnov_coordinates(mag, machine_dict['AMP2'][key])
        mirnov_names.append(key)
        B_meas, _ = probes_dict[key]
        myrecon.add_Mirnov(position, orientation, B_meas, err=0.1*max(B_mean, abs(B_meas)))
    
    return magSensors322

def load_flux_loop(mygs, myrecon, loops_dict):    
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
        if loops_dict[key][1] == 0:
            print("No data for loop {}".format(key))
            continue
        flux_names.append(key)
        B_tmp, _ = loops_dict[key]
        flux_locs.append(fl)
        flux_vals.append(B_tmp)
        psi_val = B_tmp*2.0*np.pi
        myrecon.add_flux_loop(fl, psi_val, err=0.1*(max(psi_mean_1, abs(psi_val))))

    for key, fl in fluxLoops_vv.items():
        if loops_dict[key][1] == 0:
            print("No data for loop {}".format(key))
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

def load_coil(mygs, e_coil_dict, f_coil_dict, machine_dict):
    target_currents = {}
    target_errs = {}

    for key in e_coil_dict:
        if e_coil_dict[key][1] == 0:
            continue
        target_currents[key] = e_coil_dict[key][0]
        target_errs[key] = abs(e_coil_dict[key][1])

    for key in f_coil_dict:
        if f_coil_dict[key][1] == 0:
            continue
        target_currents[key] = f_coil_dict[key][0] / machine_dict['FCOIL'][key][4]
        target_errs[key] = abs(f_coil_dict[key][1] / machine_dict['FCOIL'][key][4])

    # Set coil regularization to weakly track measured coil currents
    regularization_terms = []
    for key in e_coil_dict:
        if e_coil_dict[key][1] == 0:
            continue
        regularization_terms.append(mygs.coil_reg_term({key: 1.0},
                                                    target=target_currents[key],
                                                    weight=1.0E1))
    for key in f_coil_dict:
        if f_coil_dict[key][1] == 0:
            continue
        regularization_terms.append(mygs.coil_reg_term({key: 1.0},
                                                    target=target_currents[key],
                                                    weight=1.0E2))

    # Set zero target current and small weight on virtual VSC to allow up-down adjustment
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-2))

    # Pass regularization terms to TokaMaker
    mygs.set_coil_reg(reg_terms=regularization_terms)    

def eq_reconstruct(mygs, myrecon, ffp_prof, eqdsk, f_SOL=True):
    Ip_target=abs(eqdsk['ip'])
    P0_target=eqdsk['pres'][0]
    mygs.set_targets(Ip=Ip_target, pax=P0_target)

    # Set reconstruction settings
    myrecon.settings.fitPnorm = False
    myrecon.settings.fitR0 = True
    myrecon.settings.fitV0 = True
    myrecon.settings.fitCoils = True
    myrecon.settings.pm = False

    if f_SOL:
        pprime = [1.0, 0.0, 0.0]
        psi_pprime = [0.0, 1.0, 1.5]
        psi_sample = np.linspace(0.0,1.5,100)
        psi_prof = psi_sample.copy()
    else:
        pprime = [1.0, 0.0]
        psi_pprime = [0.0, 1.0]
        psi_sample = np.linspace(0.0,1.0,100)
        psi_prof = psi_sample.copy()
    pp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_pprime,pprime)))).copy()
    mygs.set_profiles(pp_prof={'type': 'linterp', 'y': pp_prof[:,1], 'x': psi_sample}, f_SOL=False)
    mygs.set_profiles(ffp_prof={'type': 'linterp', 'y': ffp_prof[:,1], 'x': psi_sample},pp_prof={'type': 'linterp', 'y': pp_prof[:,1], 'x': psi_sample},f_SOL=False)

    # Initial equilibrium with very rough guess
    R0 = eqdsk['rcentr']
    Z0 = 0.0
    a = 0.6
    kappa = 1.13
    delta = .1

    err_flag = mygs.init_psi(R0, Z0, a, kappa, delta)
    if err_flag:
        raise(Exception("Initialize PSI Failed."))

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
    
    chi_sq = 0.0
    with open('fit.out') as fit:
        for line in fit:
            chi = float(line.split()[1])
            chi_sq += chi ** 2
    os.remove('fit.out')

    return chi_sq

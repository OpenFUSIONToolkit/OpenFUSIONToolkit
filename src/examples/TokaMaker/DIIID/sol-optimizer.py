from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
import numpy as np
import sys
import traceback
from collections import OrderedDict
from sol import *

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

probe_coords = load_mirnov(myrecon, probes_dict,machine_dict)
loop_coords, flux_locs, flux_vals = load_flux_loop(mygs, myrecon, loops_dict)
load_coil(mygs, e_coil_dict, f_coil_dict, machine_dict)

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
        err = eq_reconstruct(mygs, myrecon, ffp_prof, probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk, probe_coords, loop_coords)
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
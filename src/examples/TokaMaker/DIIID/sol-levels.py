from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
import numpy as np
import matplotlib.pyplot as plt
import sys
import traceback
from collections import OrderedDict
from sol import *

def gen_modified_eqdsk(sol_height, eqdsk, graph=False):
    ffprim = np.array(eqdsk['ffprim'])
    ffprim /= ffprim[0]
    ffprim -= ffprim[-1]
    psi_eqdsk = np.linspace(0.0,1.0,np.size(ffprim))
    ffprim = np.r_[ffprim, ffprim[-1]]
    ffprim += sol_height
    psi_eqdsk = np.r_[psi_eqdsk, 1.5]

    psi_sample = np.linspace(0.0,1.5,150)
    psi_prof = np.copy(psi_sample)
    ffp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_eqdsk,ffprim)))).copy()

    if graph:
        plt.plot(ffp_prof[:,0], ffp_prof[:,1])
        plt.show()

    return ffp_prof

# Initialize and Perform Optimization
if len(sys.argv) < 2:
    print("USAGE: python sol-levels.py <SHOT_TAG>")
    sys.exit(0)
shot_tag = sys.argv[1]

print("=== LOADING SHOT FILE {} ===".format(shot_tag))
myOFT = OFT_env(nthreads=2)
mygs = TokaMaker(myOFT)
mygs.settings.maxits=1000
myrecon = reconstruction(mygs)
probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk = load_shot(mygs, shot_tag)
mygs.update_settings()

probe_coords = load_mirnov(myrecon, probes_dict, machine_dict)
loop_coords, flux_locs, flux_vals = load_flux_loop(mygs, myrecon, loops_dict)
load_coil(mygs, e_coil_dict, f_coil_dict, machine_dict)

print("=== BEGINNING OPTIMIZATION ===")
parameter_space = np.linspace(0.0, 1.0, 6)

min_err = float('inf')
min_params = None
i = 0
for param in parameter_space:
    try:
        ffp_prof = gen_modified_eqdsk(param, eqdsk, graph=False)
        print("RUNNING OPTIIMIZATION [{}/{}]".format(i + 1, len(parameter_space)))
        # Reset Flux Loops
        mygs.set_flux(flux_locs, flux_vals)
        # Run Reconstruction
        recons_err = eq_reconstruct(mygs, myrecon, ffp_prof, probes_dict, loops_dict, e_coil_dict, f_coil_dict, machine_dict, eqdsk, probe_coords, loop_coords)
        # Check Error
        if recons_err < min_err:
            min_err = recons_err
            min_params = param
        print("Last Err = {}".format(recons_err))
        print("Last Parameters = {}".format(param))
        print("Min Err = {}".format(min_err))
        print("Min Parameters = {}".format(min_params))
        print()
    except Exception as e:
        print("Solve failed for step={} parameters={}".format(i + 1, param))
        traceback.print_exc()
    finally:
        i = i + 1

print("=== OPTIMIZATION COMPLETE ===")
print("Min Err = {}".format(min_err))
print("Min Parameters = {}".format(min_params))
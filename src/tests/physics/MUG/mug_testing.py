#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Shared helpers for MUG tests: validation of error norms and plot output
import math
import h5py


def validate_errors(results_file, expected, factor):
    '''Check errors written by a test program, one per line, against their expected values

    `expected` is a sequence of (name, expected error) pairs in file order, and each error
    may exceed its expected value by at most `factor`. Non-finite errors always fail.
    '''
    retval = True
    with open(results_file, 'r') as fid:
        for name, err_exp in expected:
            err_test = float(fid.readline())
            if not math.isfinite(err_test) or err_test > factor*err_exp:
                print("FAILED: {0} error too high!".format(name))
                print("  Expected = {0:.8E}".format(err_exp))
                print("  Actual =   {0:.8E}".format(err_test))
                retval = False
    return retval


def validate_mug_output(steps_exp=11, linear=False, two_temp=False):
    '''Check the number of time steps and fields in the MUG plot group of `oft_xdmf.0001.h5`

    Linear runs must also output the equilibrium field (`B0`) and two-temperature runs the
    electron temperature (`Te`).
    '''
    retval = True
    step_count = 0
    B0_found = False
    Te_found = False
    with h5py.File("{0}.{1}.h5".format('oft_xdmf',str(1).zfill(4)),'r') as h5_file:
        if 'mug' in h5_file:
            for _, mesh_obj in h5_file['mug'].items():
                if mesh_obj['TYPE'][()] > 30:
                    for i in range(9999):
                        timestep = mesh_obj.get('{0:04d}'.format(i),None)
                        if timestep is None:
                            break
                        step_count += 1
                        if 'B0' in timestep:
                            B0_found = True
                        if 'Te' in timestep:
                            Te_found = True
        else:
            print('FAILED: "mug" plot group not found in output file')
            retval = False
    if step_count != steps_exp:
        print("FAILED: Incorrect number of time steps!")
        print("  Expected = {0}".format(steps_exp))
        print("  Actual =   {0}".format(step_count))
        retval = False
    if linear and (not B0_found):
        print("FAILED: Linearization plot fields missing for linear run!")
        retval = False
    if two_temp and (not Te_found):
        print("FAILED: Te plot field missing for two-temperature run!")
        retval = False
    return retval

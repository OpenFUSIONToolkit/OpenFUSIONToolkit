'''! Helper interfaces for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''

## @file util.py
#
# Helper interfaces for Open FUSION Toolkit (OFT) Python interfaces
import os
import subprocess
import numpy
import h5py
from ._interface import oftpy_dump_cov

# Common parameters
## Vacuum magnetic permeability
mu0 = numpy.pi*4.E-7
## Electron charge
eC = 1.60217663e-19


def run_shell_command(command, timeout=10, env_vars={}):
    '''! Run a shell command

    @param command Command and arguments to run
    @param timeout Timeout for command to complete [seconds]
    @param env_vars Modifications to runtime environment
    @returns `result` STDOUT, `errcode` Error code from command
    '''
    # Run shell command
    my_env = os.environ.copy()
    for key, val in env_vars.items():
        my_env[key] = val
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=my_env)
    # Wait for process to complete or timeout
    try:
        outs, _ = pid.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        pid.kill()
        outs, _ = pid.communicate()
        print("WARNING: Command timeout")
    errcode = pid.poll()
    result = outs.decode("utf-8")
    return result, errcode


def write_native_mesh(filename, r, lc, reg, nodesets=[], sidesets=[], ho_info=None, periodic_info=None):
    r'''Create a native HDF5 mesh file for OFT from the given mesh information

    @param filename Filename for mesh file
    @param r Points list [np,3]
    @param lc Cell list [nc,3] (1-based)
    @param reg Region list [nc]
    @param nodesets List of node sets
    @param sidesets List of side sets
    @param ho_info High-order grid information
    @param periodic_info Information for mesh periodicity
    '''
    print()
    print("Saving mesh: {0}".format(filename))
    with h5py.File(filename, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
        if len(nodesets) > 0:
            h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(nodesets),], dtype='i4')
            for i, node_set in enumerate(nodesets):
                h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set, dtype='i4')
        if len(sidesets) > 0:
            h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(sidesets),], dtype='i4')
            for i, side_set in enumerate(sidesets):
                h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set, dtype='i4')
        if ho_info is not None:
            h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
            h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1], dtype='i4')
            if ho_info[2] is not None:
                h5_file.create_dataset('mesh/ho_info/LF', data=ho_info[2], dtype='i4')
        if periodic_info is not None:
            h5_file.create_dataset('mesh/periodicity/nodes', data=periodic_info, dtype='i4')
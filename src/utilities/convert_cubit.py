from __future__ import print_function
import argparse
import os
import numpy as np
import h5py
import netCDF4

def read_mesh(filename):
    print()
    print("Reading mesh: {0}".format(filename))
    ncdf_file = netCDF4.Dataset(filename, "r")
    # Read points
    if 'coord' in ncdf_file.variables:
        r = np.asarray(ncdf_file.variables['coord']).transpose().copy()
    else:
        r = []
        r.append(np.asarray(ncdf_file.variables['coordx']))
        r.append(np.asarray(ncdf_file.variables['coordy']))
        if 'coordz' in ncdf_file.variables:
            rz = np.asarray(ncdf_file.variables['coordz'])
            if np.max(abs(rz)) > 0.0:
                r.append(rz)
        r = np.transpose(np.vstack(r)).copy()
    print("  Dimension    = {0}".format(r.shape[1]))
    print("  # of points  = {0}".format(r.shape[0]))
    # Read regions
    lc = []
    reg = []
    node_sets = []
    side_sets = []
    for varname in ncdf_file.variables:
        if varname.startswith('connect'):
            lc_tmp = np.asarray(ncdf_file.variables[varname])
            lc.append(lc_tmp)
            nReg = len(lc)
            reg.append(nReg*np.ones((lc_tmp.shape[0],), dtype=np.int32))
        elif varname.startswith('node_ns'):
            node_sets.append(np.asarray(ncdf_file.variables[varname]))
        elif varname.startswith('elem_ss'):
            side_sets.append(np.asarray(ncdf_file.variables[varname]))
    nReg = len(lc)
    lc = np.vstack(lc)
    reg = np.hstack(reg)
    print("  # of cells   = {0}".format(lc.shape[0]))
    print("  # of regions = {0}".format(nReg))
    print("  # of nSets   = {0}".format(len(node_sets)))
    print("  # of sSets   = {0}".format(len(side_sets)))
    #
    return r, lc, reg, node_sets, side_sets

def write_file(filename, r, lc, reg, node_sets=[], side_sets=[]):
    print()
    print("Saving mesh: {0}".format(filename))
    h5_file = h5py.File(filename, 'w')
    h5_file.create_dataset('mesh/R', data=r, dtype='f8')
    h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
    h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
    if len(node_sets) > 0:
        h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(node_sets),], dtype='i4')
        for i, node_set in enumerate(node_sets):
            h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set, dtype='i4')
    if len(side_sets) > 0:
        h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(side_sets),], dtype='i4')
        for i, side_set in enumerate(side_sets):
            h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set, dtype='i4')


parser = argparse.ArgumentParser()
parser.description = "Pre-processing script for mesh files"
parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
parser.add_argument("--out_file", type=str, default=None, help="Ouput mesh file")
options = parser.parse_args()

out_file = options.out_file
if out_file is None:
    out_file = os.path.splitext(options.in_file)[0] + ".h5"

r, lc, reg, node_sets, side_sets = read_mesh(options.in_file)
write_file(out_file, r, lc, reg, node_sets, side_sets)
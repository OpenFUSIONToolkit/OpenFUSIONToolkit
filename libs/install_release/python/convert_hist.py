#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! @file convert_hist.py

 Script to convert OFT history files to MATLAB or HDF5'
'''
from __future__ import print_function
import argparse
from OpenFUSIONToolkit.io import histfile


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Pre-processing script for mesh files"
    parser.add_argument("--files", type=str, default=None, nargs='+', required=True, help="Files to view or convert")
    parser.add_argument("--convert_hdf5", action="store_true", default=False, help="Convert files to HDF5? (default: False)")
    parser.add_argument("--convert_matlab", action="store_true", default=False, help="Convert files to MATLAB? (default: False)")
    options = parser.parse_args()

    for file in options.files:
        hist_file = histfile(file)
        file_prefix = file.split('.')[0]
        if options.convert_hdf5:
            hist_file.save_to_hdf5(file_prefix + ".h5")
        if options.convert_matlab:
            hist_file.save_to_matlab(file_prefix + ".mat")
        if not (options.convert_hdf5 or options.convert_matlab):
            print(hist_file)
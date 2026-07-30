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
from OpenFUSIONToolkit.io import _convert_hist_cli


if __name__ == "__main__":
    _convert_hist_cli()
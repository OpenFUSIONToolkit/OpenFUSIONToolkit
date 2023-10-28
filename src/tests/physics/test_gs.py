from __future__ import print_function
import os
import sys
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='square'
 cad_type=92
 nlevels=1
 nbase=1
 grid_order=1
/

&cube_options
 mesh_type=1
 ni=5,5,0
 shift=0.0,-0.5,0.0
/

&tokamaker_options
 order=2
 pm=T
 alam=4.9
 pnorm={0}
 urf=0.1
 maxits=50
 has_plasma=T
 nl_tol=1.E-11
/
"""

def gs_setup(pnorm,f_prof="flat\n",p_prof="flat\n"):
    """
    Common setup and run operations for thin-wall physics module test cases
    """
    # Create main input file from template
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(pnorm))
    # Create F profile specification
    with open('f_prof.in','w+') as fid:
        fid.write(f_prof)
    # Create P profile specification
    with open('p_prof.in','w+') as fid:
        fid.write(p_prof)
    # Run gs driver
    return run_OFT("../../tokamaker_gs oft.in", 1, 180, return_stdout=True)

def validate_gs(run_output, expected, tol=1.E-2):
    """
    Helper function to validate G-S run output for test case.
    """
    map_dict = {
        "q_0, q_95, q_a": ("q_0", "q_95"),
        "Toroidal Current": ("Ip",),
        "Magnetic Axis": ("R_0", "Z_0"),
        "Elongation": ("kappa",),
        "Triangularity": ("delta",),
        "Stored Energy": ("W_MHD",),
        "<Beta_pol>": ("beta_p",)
    }
    gs_data = {key: None for key in expected}
    found_stats=False
    for line in run_output.splitlines():
        line_stripped = line.strip()
        if line_stripped.startswith("Equilibrium Statistics:"):
            found_stats=True
            continue
        if found_stats:
            for key in map_dict:
                if line_stripped.startswith(key):
                    vals = line_stripped.split("=")[1].split()
                    for (i, key_inner) in enumerate(map_dict[key]):
                        gs_data[key_inner] = float(vals[i])
    retval = True
    for key in expected:
        if gs_data[key] is None:
            print('FAILED: Variable "{0}" not found in output!'.format(key))
            retval = False
            continue
        if abs(gs_data[key]-expected[key])/max(1.E-10,abs(expected[key])) > tol:
            print('FAILED: Variable "{0}" incorrect!'.format(key))
            print("  Expected = {0}".format(expected[key]))
            print("  Actual =   {0}".format(gs_data[key]))
            retval = False
    return retval

#============================================================================
# Test runners for time-dependent cases
def test_taylor():
    expected_values = {
        "q_0": 0.648, "q_95": 0.444,
        "Ip": 1.405E+07,
        "R_0": 0.627, "Z_0": 0.0,
        "kappa": 1.341, "delta": -0.168,
        "W_MHD": 0.0, "beta_p": 0.0
    }
    succ_flag, run_output = gs_setup(0.0)
    assert succ_flag
    assert validate_gs(run_output, expected_values)
def test_sph_pressure():
    expected_values = {
        "q_0": 0.640, "q_95": 0.438,
        "Ip": 1.415E+07,
        "R_0": 0.631, "Z_0": 0.0,
        "kappa": 1.354, "delta": -0.168,
        "W_MHD": 1.561E+06, "beta_p": 3.281
    }
    succ_flag, run_output = gs_setup(1.0)
    assert succ_flag
    assert validate_gs(run_output, expected_values)
def test_sph_linlam():
    expected_values = {
        "q_0": 0.527, "q_95": 0.357,
        "Ip": 1.174E+07,
        "R_0": 0.636, "Z_0": 0.0,
        "kappa": 1.373, "delta": -0.164,
        "W_MHD": 0.0, "beta_p": 0.0
    }
    succ_flag, run_output = gs_setup(0.0,f_prof="linear\n 1\n 0.2\n")
    assert succ_flag
    assert validate_gs(run_output, expected_values)
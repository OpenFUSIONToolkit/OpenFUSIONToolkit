#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: vacuum field in a torus driven by one toroidal flux jump plane through the Python interface
# (curved mesh), checked against regression values of |B| at a single point
import pytest
from marklin_testing import mp_run, clear_cache, run_marklin_vac


def marklin_torus_vac(order,mag_val):
    if order == 1:
        clear_cache()
    results = mp_run(run_marklin_vac,('torus_test.h5',order,2))
    assert results is not None
    assert abs((results-mag_val)/mag_val) < 1.E-4

@pytest.mark.coverage
def test_marklin_vac_g2_p1():
    '''Linear elements on curved (quadratic geometry) cells.'''
    field_mag = 1.145501
    marklin_torus_vac(1,field_mag)

@pytest.mark.coverage
def test_marklin_vac_g2_p2():
    '''Quadratic elements on curved (quadratic geometry) cells.'''
    field_mag = 1.188539
    marklin_torus_vac(2,field_mag)

def test_marklin_vac_g2_p3():
    '''Cubic elements on curved (quadratic geometry) cells.'''
    field_mag = 1.188352
    marklin_torus_vac(3,field_mag)

def test_marklin_vac_g2_p4():
    '''Quartic elements on curved (quadratic geometry) cells.'''
    field_mag = 1.188126
    marklin_torus_vac(4,field_mag)

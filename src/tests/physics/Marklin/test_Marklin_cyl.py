#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: Taylor states in a cylinder on straight (linear) and curved (quadratic) meshes
import pytest
from marklin_testing import mp_run, clear_cache, run_marklin, validate_eigs


# Test runners for cylinder taylor states (straight mesh)
def marklin_lin_cell(order):
    eigs = [5.029946,4.958005,4.957746,4.957744]
    if order == 1:
        clear_cache()
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,1,None))
    assert validate_eigs(results,eigs[order-1])

@pytest.mark.coverage
def test_marklin_g1_p1():
    marklin_lin_cell(1)

@pytest.mark.coverage
def test_marklin_g1_p2():
    marklin_lin_cell(2)

def test_marklin_g1_p3():
    marklin_lin_cell(3)

def test_marklin_g1_p4():
    marklin_lin_cell(4)

# Test runners for cylinder taylor states (curved mesh)
def marklin_quad_cell(order):
    eigs = [5.027066,4.955207,4.954955,4.954955]
    if order == 1:
        clear_cache()
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,2,None))
    assert validate_eigs(results,eigs[order-1])

@pytest.mark.coverage
def test_marklin_g2_p1():
    marklin_quad_cell(1)

@pytest.mark.coverage
def test_marklin_g2_p2():
    marklin_quad_cell(2)

def test_marklin_g2_p3():
    marklin_quad_cell(3)

def test_marklin_g2_p4():
    marklin_quad_cell(4)

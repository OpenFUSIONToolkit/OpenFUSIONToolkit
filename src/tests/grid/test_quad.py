from __future__ import print_function
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels=1
 nbase=1
 grid_order=1
/

&cube_options
 mesh_type={1}
/

&test_quad_options
 order={0}
/
"""

def quad_case(order,grid_type):
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(order,grid_type))
    result = run_OFT("./test_quad", 1, 60)
    if result:
        #
        errs = [0.,0.,0.]
        for i in range(3):
            fid = open('quad_{0}d.results'.format(i+1),'r')
            for line in fid:
                errs[i] = max(errs[i],abs(float(line)))
            fid.close()
        return errs
    else:
        return [99., 99., 99.]

@pytest.mark.parametrize("grid_type", range(1,3))
@pytest.mark.parametrize("order", range(1,12))
def test_quadrature(order,grid_type):
    tol = 1.e-12
    errs = quad_case(order,grid_type)
    result = True
    for (i, err) in enumerate(errs):
        if abs(err) >= tol:
            print("FAILED: Quadrature dimension = {0}, err = {1}".format(i+1, err))
            result = False
    assert result

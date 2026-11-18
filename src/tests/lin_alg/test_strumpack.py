import os
import sys
import pytest
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

&test_strumpack_options
 nsize={0}
 comp_type={1}
 leaf_size={2}
 min_sep_size={3}
 comp_rel_tol={4}
 use_xml={5}
/
"""

# Template for XML input file
xml_template = """
<oft>
  <solver type="lu">
    <package>strumpack</package>
    <strumpack>
      <comp_type>{0}</comp_type>
      <comp_rel_tol>{1}</comp_rel_tol>
      <leaf_size>{2}</leaf_size>
      <min_sep_size>{3}</min_sep_size>
    </strumpack>
  </solver>
</oft>
"""

# Common setup function and process handling
def strumpack_setup(nsize, comp_type, leaf_size=-1, min_sep_size=-1, comp_rel_tol=-1.0, xml_comp_type=None):
    os.chdir(test_dir)
    use_xml = 'F'
    command = "./test_strumpack"
    if xml_comp_type is not None:
        use_xml = 'T'
        command = "./test_strumpack oft.in oft_in.xml"
        with open('oft_in.xml','w+') as fid:
            fid.write(xml_template.format(xml_comp_type, comp_rel_tol, leaf_size, min_sep_size))
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nsize, comp_type, leaf_size, min_sep_size, comp_rel_tol, use_xml))
    return run_OFT(command, 1, 60)

def validate_result(max_error, xml_options=None):
    with open('strumpack.results','r') as fid:
        errors = [float(fid.readline()) for _ in range(4)]
        if xml_options is not None:
            package = fid.readline().strip()
            options = [int(fid.readline()), float(fid.readline()), int(fid.readline()), int(fid.readline())]
    for i, error in enumerate(errors):
        if error > max_error:
            print("FAILED: Solution error too large (solve {0})!".format(i+1))
            print("  Expected < {0}".format(max_error))
            print("  Actual =   {0}".format(error))
            return False
    if xml_options is not None:
        if package != 'strumpack':
            print("FAILED: Incorrect package from XML ({0})".format(package))
            return False
        if options != xml_options:
            print("FAILED: Incorrect STRUMPACK options from XML!")
            print("  Expected = {0}".format(xml_options))
            print("  Actual =   {0}".format(options))
            return False
    return True

#============================================================================
# Test runner for sparse solves (0 -> none, 1 -> HSS, 2 -> BLR)
@pytest.mark.coverage
@pytest.mark.parametrize("comp_type", (0, 1, 2))
def test_sparse(comp_type):
    assert strumpack_setup(60, comp_type, leaf_size=16, min_sep_size=16, comp_rel_tol=1.E-10)
    assert validate_result(1.E-8)

# Test runner for options read from XML
@pytest.mark.coverage
@pytest.mark.parametrize("comp_type", (("none", 0), ("hss", 1), ("blr", 2)))
def test_sparse_xml(comp_type):
    assert strumpack_setup(60, -1, leaf_size=16, min_sep_size=16, comp_rel_tol=1.E-10, xml_comp_type=comp_type[0])
    assert validate_result(1.E-8, xml_options=[comp_type[1], 1.E-10, 16, 16])

import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn={2}
 debug=0
 test_run=T
 use_petsc={0}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels=3
 nbase={1}
 grid_order=1
/

&cube_options
 mesh_type=1
/
"""
#
xml_template = """
<oft>
  {0}
</oft>
"""

# Common setup function and process handling
def xml_setup(xml_file, petsc=False, nprocs=1, ppn=1, cell_parts=0):
    os.chdir(test_dir)
    petsc_flag=('T' if petsc else 'F')
    nbase = 3
    if nprocs > 1:
        nbase = 2 # Distribute finest level for MPI runs
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(petsc_flag,nbase,ppn))
        if cell_parts > 0:
            fid.write("\n&test_solver_xml_options\n cell_parts={0}\n/\n".format(cell_parts))
    with open('oft_in.xml','w+') as fid:
        fid.write(xml_template.format(xml_file))
    return run_OFT("./test_solver_xml oft.in oft_in.xml", nprocs, 60)

def validate_result(iteration_count,converged_error):
    #
    fid = open('xml.results','r')
    its_test = int(fid.readline())
    if iteration_count is not None:
        if abs(iteration_count-its_test) >= max(3,0.05*iteration_count):
            print("FAILED: Iteration count incorrect!")
            print("  Expected = {0}".format(iteration_count))
            print("  Actual =   {0}".format(its_test))
            return False
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-6:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_gmres(petsc_flag):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
    </solver>
    """
    iteration_count = 530
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_gmres_jacobi(petsc_flag):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="jacobi"></pre>
    </solver>
    """
    iteration_count = 492
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_cg(petsc_flag):
    xml_file = """
    <solver type="cg">
        <its>-3</its>
    </solver>
    """
    iteration_count = 145
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_cg_jacobi(petsc_flag):
    xml_file = """
    <solver type="cg">
        <its>-3</its>
        <pre type="jacobi"></pre>
    </solver>
    """
    iteration_count = 154
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# PETSc additive Schwarz/block Jacobi with explicitly defined local subdomains
# (nlocal=-1 uses field-based subdomains)
@pytest.mark.parametrize("nlocal,overlap,iteration_count", (
    (-1, 0, 1),
    (2, 0, 100),
    pytest.param(2, 1, 1, marks=pytest.mark.coverage),
))
def test_gmres_petsc_asm(nlocal, overlap, iteration_count):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>{0}</nlocal>
            <overlap>{1}</overlap>
        </pre>
    </solver>
    """.format(nlocal, overlap)
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=True)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# PETSc additive Schwarz/block Jacobi with explicitly defined local subdomains (MPI)
# (full_node=1 uses node subdomains with PCGASM, but with ppn=1 each node is a single process)
@pytest.mark.mpi
@pytest.mark.parametrize("nlocal,overlap,full_node,iteration_count", (
    (-1, 0, 0, 38),
    (1, 1, 0, 8),
    (2, 0, 0, 61),
    (2, 1, 0, 15),
    (-1, 0, 1, 38),
    (1, 1, 1, 8),
))
def test_gmres_petsc_asm_mpi(nlocal, overlap, full_node, iteration_count):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>{0}</nlocal>
            <overlap>{1}</overlap>
            <full_node>{2}</full_node>
        </pre>
    </solver>
    """.format(nlocal, overlap, full_node)
    converged_error = 0.6394036978170452
    assert xml_setup(xml_file,petsc=True,nprocs=2)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# PETSc GASM with node subdomains spanning both processes (ppn=2), using a
# symmetric Jacobi sub-solver as the subdomain matrices are parallel
# (nlocal=-1 uses field-based subdomains; the XML input sets overlap=0 for nlocal=-1)
@pytest.mark.mpi
@pytest.mark.parametrize("nlocal,overlap,iteration_count", (
    (-1, 0, 42),
    pytest.param(1, 1, 39, marks=pytest.mark.coverage),
))
def test_gmres_petsc_gasm_node_mpi(nlocal, overlap, iteration_count):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>{0}</nlocal>
            <overlap>{1}</overlap>
            <full_node>1</full_node>
            <solver type="sjacobi"><its>2</its></solver>
        </pre>
    </solver>
    """.format(nlocal, overlap)
    converged_error = 0.6394036978170452
    assert xml_setup(xml_file,petsc=True,nprocs=2,ppn=2)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Native overlapping Schwarz selected from XML: blocks within the process grown by one layer
# (plain additive combination counts overlap rows twice, hence more iterations)
@pytest.mark.parametrize("schwarz_type", ("restricted", "additive"))
def test_gmres_schwarz(schwarz_type):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>4</nlocal>
            <overlap>1</overlap>
            <type>{0}</type>
            <pre type="ilu"></pre>
        </pre>
    </solver>
    """.format(schwarz_type)
    iteration_count = {"restricted": 58, "additive": 102}[schwarz_type]
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=False)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Native Schwarz in cell mode selected from XML: 4 blocks (from the cell overlap builder
# attached to the matrix) grown by one layer of face-adjacent cells (with ILU sub-solves the
# counts are close to graph mode; the overlap trends are checked in test_cell_schwarz)
@pytest.mark.parametrize("schwarz_type", ("restricted", "additive"))
def test_gmres_schwarz_cell(schwarz_type):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <overlap>1</overlap>
            <overlap_mode>cell</overlap_mode>
            <adjacency>face</adjacency>
            <base_layers>2</base_layers>
            <type>{0}</type>
            <pre type="ilu"></pre>
        </pre>
    </solver>
    """.format(schwarz_type)
    iteration_count = {"restricted": 58, "additive": 100}[schwarz_type]
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=False,cell_parts=4)
    assert validate_result(iteration_count,converged_error)

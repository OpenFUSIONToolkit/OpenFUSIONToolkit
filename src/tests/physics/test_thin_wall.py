from __future__ import print_function
import os
import sys
import pytest
import numpy as np
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT
from oft_io import oft_histfile

mu0 = np.pi*4.E-7

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 cad_type=0
/

&native_mesh_options
 filename="{0}"
/

&cubit_options
 filename="{0}"
 lf_file=T
/

&thincurr_td_options
 {4}
 {5}
 dt=2.e-5
 nsteps=200
 nplot=10
 direct={1}
 cg_tol=1.E-8
 save_L=F
 save_Mcoil=F
 plot_run=F
/

&thincurr_eig_options
 direct={1}
 plot_run=F
 neigs={6}
 reduce_model={7}
/

&thincurr_fr_options
 direct={1}
 freq={2}
 fr_limit={3}
/
"""

oft_in_xml_template_template = """
<oft>
  <thincurr>
    <eta>{1}</eta>
    {0}
  </thincurr>
</oft>
"""

def thin_wall_setup(meshfile,run_type,direct_flag,freq=0.0,fr_limit=0,eta=10.0,
                    icoils=None,vcoils=None,floops=None,curr_waveform=None,volt_waveform=None):
    """
    Common setup and run operations for thin-wall physics module test cases
    """
    nPhi = 180
    phi_fac = 2.0*np.pi/(nPhi-1)
    # Create main input file from template
    os.chdir(test_dir)
    if curr_waveform is None:
        coil_file_line=""
    else:
        coil_file_line='curr_file="curr.drive"' 
    if volt_waveform is None:
        volt_file_line=""
    else:
        volt_file_line='volt_file="volt.drive"' 
    neigs = 4
    reduce_model_flag = 'F'
    if run_type == 4:
        neigs = 10
        reduce_model_flag = 'T'
        run_type = 2
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(meshfile,direct_flag,freq,fr_limit,coil_file_line,volt_file_line,neigs,reduce_model_flag))
    # Create XML input file for coils
    coil_string = ""
    if icoils is not None:
        coil_string += '<icoils><coil_set type="2">\n'
        for icoil in icoils:
            coil_string += '<coil npts="{0}">'.format(nPhi)
            R = icoil[0]; Z = icoil[1]
            for i in range(nPhi):
                phi = i*phi_fac
                coil_string += '{0:.12E} {1:.12E} {2:.12E}\n'.format(R*np.cos(phi), R*np.sin(phi), Z)
            coil_string += "</coil>\n"
        coil_string += "</coil_set></icoils>"
    if vcoils is not None:
        coil_string += '<vcoils>\n'
        for pcoil in vcoils:
            coil_string += '<coil_set type="2" res_per_len="1.256637E-5" radius="1.E-2"><coil npts="{0}">\n'.format(nPhi)
            R = pcoil[0]; Z = pcoil[1]
            for i in range(nPhi):
                phi = i*phi_fac
                coil_string += '{0:.12E} {1:.12E} {2:.12E}\n'.format(R*np.cos(phi), R*np.sin(phi), Z)
            coil_string += "</coil></coil_set>\n"
        coil_string += "</vcoils>"
    with open('oft_in.xml','w+') as fid:
        fid.write(oft_in_xml_template_template.format(coil_string, eta*mu0))
    # Create flux loop definition file for sensors
    if floops is not None:
        with open('floops.loc', 'w+') as fid:
            fid.write('{0}\n'.format(len(floops)))
            #
            for (k, floop) in enumerate(floops):
                fid.write('\n{0} 1.0 FLOOP_{1}\n'.format(nPhi, k))
                R = floop[0]; Z = floop[1]
                for i in range(nPhi):
                    phi = i*phi_fac
                    fid.write('{0:.6E} {1:.6E} {2:.6E}\n'.format(R*np.cos(phi), R*np.sin(phi), Z))
    # Create coil drive waveforms
    if curr_waveform is not None:
        n = len(curr_waveform)
        m = len(curr_waveform[0])
        with open('curr.drive', 'w+') as fid:
            fid.write('{0} {1}\n'.format(m,n))
            for i in range(n):
                fid.write(' '.join(['{0}'.format(curr_waveform[i][0])] + ['{0}'.format(val/mu0) for val in curr_waveform[i][1:]])+'\n')
    if volt_waveform is not None:
        n = len(volt_waveform)
        m = len(volt_waveform[0])
        with open('volt.drive', 'w+') as fid:
            fid.write('{0} {1}\n'.format(m,n))
            for i in range(n):
                fid.write(' '.join(['{0}'.format(val) for val in volt_waveform[i]])+'\n')
    # Run thin-wall model
    if run_type == 1:
        return run_OFT("../../bin/thincurr_td oft.in oft_in.xml", 1, 180)
    elif run_type == 2:
        return run_OFT("../../bin/thincurr_eig oft.in oft_in.xml", 1, 180)
    elif run_type == 3:
        return run_OFT("../../bin/thincurr_fr oft.in oft_in.xml", 1, 180)

def validate_eigs(eigs, tols=(1.E-5, 1.E-9)):
    """
    Helper function to validate eigenvalues against test case.
    """
    eigs_run_real = []
    eigs_run_imag = []
    with open('thincurr_eigs.dat', 'r') as fid:
        for line in fid:
            vals = line.split()
            eigs_run_real.append(float(vals[0]))
            eigs_run_imag.append(float(vals[1]))
    if not len(eigs_run_real) == len(eigs):
        print("FAILED: Number of eigenvalues does not match")
        return False
    retval = True
    for (i, val) in enumerate(eigs):
        if abs((val-eigs_run_real[i])/val) > tols[0]:
            print("FAILED: Eigenvalue {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(eigs_run_real[i]))
            retval = False
        if abs(eigs_run_imag[i]) > tols[1]:
            print("FAILED: Imaginary eigenvalue detected!")
            print("  Value =    {0}".format(eigs_run_real[i]))
            retval = False
    return retval

def validate_fr(fr_real, fr_imag, tols=(1.E-4, 1.E-4)):
    """
    Helper function to validate frequency-response results against test case.
    """
    hist_file = oft_histfile('thincurr_fr.dat')
    fr_run_real = [hist_file.data[field][0] for field in hist_file.field_tags]
    fr_run_imag = [hist_file.data[field][1] for field in hist_file.field_tags]
    if not len(fr_run_real) == len(fr_real):
        print("FAILED: Number of sensors does not match")
        return False
    retval = True
    for (i, val) in enumerate(fr_real):
        if abs((val-fr_run_real[i])/val) > tols[0]:
            print("FAILED: Real response {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(fr_run_real[i]))
            retval = False
    for (i, val) in enumerate(fr_imag):
        if abs((val-fr_run_imag[i])/val) > tols[1]:
            print("FAILED: Imaginary response {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(fr_run_imag[i]))
            retval = False
    return retval

def validate_td(sigs_final, tols=(1.E-8, 1.E-3)):
    """
    Helper function to validate time-dependent results against test case.
    """
    hist_file = oft_histfile('floops.hist')
    td_sigs_final = [hist_file.data[field][-1] for field in hist_file.field_tags]
    if not len(td_sigs_final) == len(sigs_final):
        print("FAILED: Number of sensors does not match")
        return False
    retval = True
    if abs(sigs_final[0]-td_sigs_final[0]) > tols[0]:
        print("FAILED: Final time incorrect!")
        print("  Expected = {0}".format(sigs_final[0]))
        print("  Actual =   {0}".format(td_sigs_final[0]))
        retval = False
    for (i, val) in enumerate(sigs_final[1:]):
        if abs((val-td_sigs_final[i+1])/val) > tols[1]:
            print("FAILED: Signal {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(td_sigs_final[i+1]))
            retval = False
    return retval

def validate_model_red(eigs, tols=(1.E-5, 1.E-9)):
    """
    Helper function to validate eigenvalues against test case.
    """
    import h5py
    with h5py.File('tCurr_reduced.h5','r') as file:
        L = np.asarray(file['L'])
        R = np.asarray(file['R'])
    eigs_run_real, _ = np.linalg.eig(np.dot(np.linalg.inv(R),L))
    eigs_run_real = -np.sort(-eigs_run_real)
    retval = True
    for (i, val) in enumerate(eigs):
        if abs((val-eigs_run_real[i])/val) > tols[0]:
            print("FAILED: Reduced model eigenvalue {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(eigs_run_real[i]))
            retval = False
    return retval

#============================================================================
# Test runners for time-dependent cases
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_plate(direct_flag):
    sigs_final = (4.E-3, 8.459371E-4, 7.130923E-4)
    assert thin_wall_setup("tw_test-plate.h5",1,direct_flag, 
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_plate_volt(direct_flag):
    sigs_final = (4.E-3, 2.383774E+1, 2.005698E+1)
    assert thin_wall_setup("tw_test-plate.h5",1,direct_flag, 
                           vcoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           volt_waveform=((0.0, 1.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_cyl(direct_flag):
    sigs_final = (4.E-3, 7.254196E-4, 6.151460E-4)
    assert thin_wall_setup("tw_test-cyl.h5",1,direct_flag, 
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_cyl_volt(direct_flag):
    sigs_final = (4.E-3, 8.559103E0, 7.268509E0)
    assert thin_wall_setup("tw_test-cyl.h5",1,direct_flag, 
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_torus(direct_flag):
    sigs_final = (4.E-3, 4.935683E-4, 3.729159E-5)
    assert thin_wall_setup("tw_test-torus.h5",1,direct_flag, 
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_torus_volt(direct_flag):
    sigs_final = (4.E-3, 3.249705E0, 2.3204651E-1)
    assert thin_wall_setup("tw_test-torus.h5",1,direct_flag, 
                           vcoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_passive(direct_flag):
   sigs_final = (4.E-3, 7.706778E-4, 7.903190E-4)
   assert thin_wall_setup("tw_test-passive.h5",1,direct_flag,eta=1.E4, 
                          icoils=((0.5, 0.1),),
                          vcoils=((0.5, 0.0),),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)))
   assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_passive_volt(direct_flag):
   sigs_final = (4.E-3, 2.114789E+1, 2.170003E+1)
   assert thin_wall_setup("tw_test-passive.h5",1,direct_flag,eta=1.E4, 
                          vcoils=((0.5, 0.0), (0.5, 0.1)),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          volt_waveform=((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)))
   assert validate_td(sigs_final)

#============================================================================
# Test runners for eigenvalue cases
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_plate(direct_flag):
    eigs = (9.735667E-3, 6.532314E-3, 6.532201E-3, 5.251598E-3)
    assert thin_wall_setup("tw_test-plate.h5",2,direct_flag)
    assert validate_eigs(eigs)
    
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_cyl(direct_flag):
    eigs = (2.657195E-2, 1.248071E-2, 1.247103E-2, 1.200566E-2)
    assert thin_wall_setup("tw_test-cyl.h5",2,direct_flag)
    assert validate_eigs(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_torus(direct_flag):
    eigs = (4.751344E-2, 2.564491E-2, 2.555695E-2, 2.285850E-2)
    assert thin_wall_setup("tw_test-torus.h5",2,direct_flag)
    assert validate_eigs(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_passive(direct_flag):
    eigs = (1.483589E-1, 6.207849E-2, 2.942791E-2, 2.693574E-2)
    assert thin_wall_setup("tw_test-passive.h5",2,direct_flag,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)))
    assert validate_eigs(eigs)

#============================================================================
# Test runners for eigenvalue-based model reduction
def test_mred_plate():
    eigs = (9.735667E-3, 6.532314E-3, 6.532201E-3, 5.251598E-3)
    assert thin_wall_setup("tw_test-plate.h5",4,True)
    assert validate_model_red(eigs)
    
def test_mred_cyl():
    eigs = (2.657195E-2, 1.248071E-2, 1.247103E-2, 1.200566E-2)
    assert thin_wall_setup("tw_test-cyl.h5",4,True)
    assert validate_model_red(eigs)

@pytest.mark.coverage
def test_mred_torus():
    eigs = (4.751344E-2, 2.564491E-2, 2.555695E-2, 2.285850E-2)
    assert thin_wall_setup("tw_test-torus.h5",4,True)
    assert validate_model_red(eigs)

@pytest.mark.coverage
def test_mred_passive():
    eigs = (1.483589E-1, 6.207849E-2, 2.942791E-2, 2.693574E-2)
    assert thin_wall_setup("tw_test-passive.h5",4,True,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)))
    assert validate_model_red(eigs)

#============================================================================
# Test runners for frequency-response cases
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_plate(direct_flag):
    fr_real = (6.807649E-2, 7.207748E-2)
    fr_imag = (-3.011666E-3, -2.177010E-3)
    assert thin_wall_setup("tw_test-plate.h5",3,direct_flag,freq=5.E3,fr_limit=0, 
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_cyl(direct_flag):
    fr_real = (6.118337E-2, 4.356188E-3)
    fr_imag = (-1.911861E-3, -2.283493E-3)
    assert thin_wall_setup("tw_test-cyl.h5",3,direct_flag,freq=5.E3,fr_limit=0, 
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_torus(direct_flag):
    fr_real = (-2.806665E-3, -1.194625E-4)
    fr_imag = (-1.869726E-3, -1.248644E-4)
    assert thin_wall_setup("tw_test-torus.h5",3,direct_flag,freq=5.E3,fr_limit=0, 
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_passive(direct_flag):
    fr_real = (1.777366E-1, 1.868689E-1)
    fr_imag = (-2.331033E-4, -1.671967E-4)
    assert thin_wall_setup("tw_test-passive.h5",3,direct_flag,eta=1.E4,freq=5.E3,fr_limit=0, 
                           icoils=((0.5, 0.1),),
                           vcoils=((0.5, 0.0),),
                           floops=((0.5, -0.05), (0.5, -0.1)))
    assert validate_fr(fr_real, fr_imag)
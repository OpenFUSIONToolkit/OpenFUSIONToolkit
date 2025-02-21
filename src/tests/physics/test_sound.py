from __future__ import print_function
import os
import sys
import h5py
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
 use_petsc={10}
/

&mesh_options
 meshname='cube'
 cad_type={11}
 nlevels={1}
 nbase={0}
 grid_order=1
/

&t3d_options
 filename='cube.t3d'
 inpname='cube.inp'
 reflect='xyz'
 ref_per=T,T,T
/

&cube_options
 mesh_type=2
 ni=1,1,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_sound_options
 order={2}
 minlev={3}
 delta=1.e-3
 linear={5}
 two_temp={12}
/

&xmhd_options
 xmhd_jcb=F
 xmhd_advec=F
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_visc_heat=F
 bbc='bc'
 vbc='all'
 dt={6}
 eta=1.
 nu_par=1.e-3
 nsteps={7}
 rst_freq=25
 lin_tol={8}
 nl_tol=1.E-5
 nu_xmhd={4}
 nclean=500
 rst_ind=0
 maxextrap=2
 kappa_par=1.d-6
 kappa_perp=1.d-6
 ittarget=100
 xmhd_mfnk={9}
/

&xmhd_plot_options
 t0=.9e-5
 dt=1.e-5
 rst_start=0
 rst_end={7}
/
"""

# Common setup function and process handling
def sound_setup(nbase,nlevels,order,minlev,nu_mhd='-1',linear=False,mf=False,petsc='F',hex_mesh=False,two_temp=False):
    os.chdir(test_dir)
    mesh_type=1
    if hex_mesh:
        mesh_type=92
    dt='1.e-6'
    its='100'
    tol='1.E-12'
    mf_flag='F'
    tt_flag='F'
    if two_temp:
        tt_flag='T'
    if linear:
        lin_flag='T'
    else:
        lin_flag='F'
        if mf:
            tol='1.E-9'
            mf_flag='T'
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nbase,nlevels,order,minlev,nu_mhd,
                                       lin_flag,dt,its,tol,mf_flag,petsc,mesh_type,tt_flag))
    return run_OFT("./test_sound", nproc, 2000)

def validate_result(nerr_exp,terr_exp,verr_exp,steps_exp=11,linear=False,two_temp=False):
    retval = True
    with open('sound.results','r') as fid:
        nerr_test = float(fid.readline())
        if nerr_test > 1.01*nerr_exp:
            print("FAILED: Density error too high!")
            print("  Expected = {0}".format(nerr_exp))
            print("  Actual =   {0}".format(nerr_test))
            retval = False
        terr_test = float(fid.readline())
        if terr_test > 1.01*terr_exp:
            print("FAILED: Temperature error too high!")
            print("  Expected = {0}".format(terr_exp))
            print("  Actual =   {0}".format(terr_test))
            retval = False
        verr_test = float(fid.readline())
        if verr_test > 1.01*verr_exp:
            print("FAILED: Velocity error too high!")
            print("  Expected = {0}".format(verr_exp))
            print("  Actual =   {0}".format(verr_test))
            retval = False
    #
    step_count=0
    B0_found = False
    Te_found = False
    with h5py.File("{0}.{1}.h5".format('oft_xdmf',str(1).zfill(4)),'r') as h5_file:
        if 'mug' in h5_file:
            for _, mesh_obj in h5_file['mug'].items():
                if mesh_obj['TYPE'][()] > 30:
                    for i in range(9999):
                        timestep = mesh_obj.get('{0:04d}'.format(i),None)
                        if timestep is None:
                            break
                        step_count += 1
                        if 'B0' in timestep:
                            B0_found = True
                        if 'Te' in timestep:
                            Te_found = True
        else:
            print('FAILED: "mug" plot group not found in output file')
            retval = False
    if step_count != steps_exp:
        print("FAILED: Incorrect number of time steps!")
        print("  Expected = {0}".format(steps_exp))
        print("  Actual =   {0}".format(step_count))
        retval = False
    if linear and (not B0_found):
        print("FAILED: Linearization plot fields missing for linear run!")
        retval = False
    if two_temp and (not Te_found):
        print("FAILED: Te plot field missing for two-temperature run!")
        retval = False
    return retval

#============================================================================
# Non-Linear test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p2(petsc_flag):
    nerr_exp = 5.9119040036454996E-002
    terr_exp = 5.9231719732675436E-002
    verr_exp = 3.6153775372848231E-002
    assert sound_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p2_mpi(mf,two_temp,petsc_flag):
    nerr_exp = 5.9119040036454996E-002
    terr_exp = 5.9231719732675436E-002
    verr_exp = 3.6153775372848231E-002
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p3(petsc_flag):
    nerr_exp = 3.4138360614025419E-003
    terr_exp = 3.4144365953483555E-003
    verr_exp = 3.3335039251361225E-003
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p3_mpi(mf,two_temp,petsc_flag):
    nerr_exp = 3.4138360614025419E-003
    terr_exp = 3.4144365953483555E-003
    verr_exp = 3.3335039251361225E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p4(petsc_flag):
    nerr_exp = 2.7069013754645772E-003
    terr_exp = 2.7074485164640862E-003
    verr_exp = 2.7098032859103993E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_nl_r1_p4_mpi(mf,two_temp,petsc_flag):
    nerr_exp = 2.7069013754645772E-003
    terr_exp = 2.7074485164640862E-003
    verr_exp = 2.7098032859103993E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p2(petsc_flag):
    nerr_exp = 5.9421723477144917E-002
    terr_exp = 5.9530267790148306E-002
    verr_exp = 3.6808185835839932E-002
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p2_mpi(two_temp,petsc_flag):
    nerr_exp = 5.9421723477144917E-002
    terr_exp = 5.9530267790148306E-002
    verr_exp = 3.6808185835839932E-002
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p3(petsc_flag):
    nerr_exp = 9.1038624886568657E-003
    terr_exp = 9.1039170745524662E-003
    verr_exp = 9.0864241920798547E-003
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p3_mpi(two_temp,petsc_flag):
    nerr_exp = 9.1038624886568657E-003
    terr_exp = 9.1039170745524662E-003
    verr_exp = 9.0864241920798547E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p4(petsc_flag):
    nerr_exp = 8.8610147661516616E-003
    terr_exp = 8.8613658988151399E-003
    verr_exp = 8.8604718664127627E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_lin_r1_p4_mpi(two_temp,petsc_flag):
    nerr_exp = 8.8610147661516616E-003
    terr_exp = 8.8613658988151399E-003
    verr_exp = 8.8604718664127627E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p2(petsc_flag):
    nerr_exp = 2.1360749524494939E-002
    terr_exp = 2.1362460079360893E-002
    verr_exp = 2.1373898556061369E-002
    assert sound_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p2_mpi(mf,two_temp,petsc_flag):
    nerr_exp = 2.1360749524494939E-002
    terr_exp = 2.1362460079360893E-002
    verr_exp = 2.1373898556061369E-002
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p3(petsc_flag):
    nerr_exp = 4.5268508665425016E-003
    terr_exp = 4.5266729911366405E-003
    verr_exp = 4.5140710240899129E-003
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p3_mpi(mf,two_temp,petsc_flag):
    nerr_exp = 4.5268508665425016E-003
    terr_exp = 4.5266729911366405E-003
    verr_exp = 4.5140710240899129E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p4(petsc_flag):
    nerr_exp = 3.2905133020702569E-003
    terr_exp = 3.2907964458522860E-003
    verr_exp = 3.2936707864742649E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_nl_r1_p4_mpi(mf,two_temp,petsc_flag):
    nerr_exp = 3.2905133020702569E-003
    terr_exp = 3.2907964458522860E-003
    verr_exp = 3.2936707864742649E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p2(petsc_flag):
    nerr_exp = 2.1008321863506766E-002
    terr_exp = 2.1008318352739700E-002
    verr_exp = 2.0998273206287381E-002
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p2_mpi(two_temp,petsc_flag):
    nerr_exp = 2.1008321863506766E-002
    terr_exp = 2.1008318352739700E-002
    verr_exp = 2.0998273206287381E-002
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p3(petsc_flag):
    nerr_exp = 9.2415214063574797E-003
    terr_exp = 9.2415774088205841E-003
    verr_exp = 9.2421965458169622E-003
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p3_mpi(two_temp,petsc_flag):
    nerr_exp = 9.2415214063574797E-003
    terr_exp = 9.2415774088205841E-003
    verr_exp = 9.2421965458169622E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p4(petsc_flag):
    nerr_exp = 8.8588190950046782E-003
    terr_exp = 8.8582074491383549E-003
    verr_exp = 8.8591147973098201E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_lin_r1_p4_mpi(two_temp,petsc_flag):
    nerr_exp = 8.8588190950046782E-003
    terr_exp = 8.8582074491383549E-003
    verr_exp = 8.8591147973098201E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

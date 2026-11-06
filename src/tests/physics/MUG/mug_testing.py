#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Shared helpers for MUG tests: input templates and runners for the 3-D wave cases, validation of error norms and
# plot output
import os
import sys
import math
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT


def validate_errors(results_file, expected, factor):
    '''Check errors written by a test program, one per line, against their expected values

    `expected` is a sequence of (name, expected error) pairs in file order, and each error
    may exceed its expected value by at most `factor`. Non-finite errors always fail.
    '''
    retval = True
    with open(results_file, 'r') as fid:
        for name, err_exp in expected:
            err_test = float(fid.readline())
            if not math.isfinite(err_test) or err_test > factor*err_exp:
                print("FAILED: {0} error too high!".format(name))
                print("  Expected = {0:.8E}".format(err_exp))
                print("  Actual =   {0:.8E}".format(err_test))
                retval = False
    return retval


def validate_mug_output(steps_exp=11, linear=False, two_temp=False):
    '''Check the number of time steps and fields in the MUG plot group of `oft_xdmf.0001.h5`

    Linear runs must also output the equilibrium field (`B0`) and two-temperature runs the
    electron temperature (`Te`).
    '''
    retval = True
    step_count = 0
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


# Input for a traveling sound wave in a triply periodic box (test_sound.F90, test_sound_lag.F90)
SOUND_TEMPLATE = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={petsc}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={nlevels}
 nbase={nbase}
 grid_order=1
/

&cube_options
 mesh_type={mesh_type}
 ni=1,1,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_sound_options
 order={order}
 minlev={minlev}
 delta=1.e-3
 linear={linear}
 two_temp={two_temp}
/

&xmhd_options
 xmhd_jcb=F
 xmhd_advec=F
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_visc_heat=F
{rep_options} vbc='all'
 dt={dt}
 eta=1.
 nu_par=1.e-3
 nsteps={nsteps}
 rst_freq=25
 lin_tol={lin_tol}
 nl_tol=1.E-5
 nu_xmhd={nu_mhd}
 nclean=500
 rst_ind=0
 maxextrap=2
 kappa_par=1.d-6
 kappa_perp=1.d-6
 ittarget=100
 xmhd_mfnk={mfnk}
/

&xmhd_plot_options
 t0=.9e-5
 dt=1.e-5
 rst_start=0
 rst_end={nsteps}
/
"""

# Input for a traveling shear Alfven wave in a triply periodic box (test_alfven.F90, test_alfven_lag.F90)
ALFVEN_TEMPLATE = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={petsc}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={nlevels}
 nbase={nbase}
 grid_order=1
/

&cube_options
 mesh_type={mesh_type}
 ni=1,2,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_alfven_options
 order={order}
 minlev={minlev}
 delta=1.e-3
 linear={linear}
/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=F
 xmhd_hall=F
 xmhd_adv_temp=F
 xmhd_adv_den=F
{rep_options} vbc='all'
 dt={dt}
 eta=1.e-3
 nu_par=1.e-3
 nsteps={nsteps}
 rst_freq=125
 lin_tol={lin_tol}
 nl_tol=1.E-5
 nu_xmhd={nu_mhd}
 nclean=5000
 rst_ind=0
 maxextrap=2
 ittarget=100
 xmhd_mfnk={mfnk}
/

&xmhd_plot_options
 t0=.9e-5
 dt=1.e-5
 rst_start=0
 rst_end={nsteps}
/
"""

# Magnetic field representation of the 3-D xMHD modules: the mixed H(curl)/Lagrange module (xmhd) sets a Dirichlet
# condition on B at the wall; the nodal Lagrange module (xmhd_lag) relies on divergence cleaning instead
REPRESENTATIONS = ('hcurl', 'lagrange')


def namelist_flag(value):
    return 'T' if value else 'F'


def mpi_ranks(nbase, nlevels):
    '''Refinement past the base level runs distributed over two ranks'''
    return 2 if nbase != nlevels else 1


def sound_setup(representation,nbase,nlevels,order,minlev,nu_mhd='-1',linear=False,mf=False,petsc=False,hex_mesh=False,two_temp=False):
    '''Write the input for a traveling sound wave and run test_sound (`representation` 'hcurl') or test_sound_lag
    ('lagrange')

    `nu_mhd` gives the multigrid smoother iterations per level. The matrix-free Newton-Krylov solver (`mf`, nonlinear
    runs only) is run with a looser linear tolerance.
    '''
    exe, rep_options, timeout = {'hcurl': ('./test_sound', " bbc='bc'\n", 2000),
                                 'lagrange': ('./test_sound_lag', " k_db=1.\n", 1000)}[representation]
    lin_tol = '1.E-12'
    mfnk = False
    if (not linear) and mf:
        lin_tol = '1.E-9'
        mfnk = True
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(SOUND_TEMPLATE.format(nbase=nbase, nlevels=nlevels, order=order, minlev=minlev, nu_mhd=nu_mhd,
                                        linear=namelist_flag(linear), dt='1.e-6', nsteps='100', lin_tol=lin_tol,
                                        mfnk=namelist_flag(mfnk), petsc=namelist_flag(petsc),
                                        mesh_type=(2 if hex_mesh else 1), two_temp=namelist_flag(two_temp),
                                        rep_options=rep_options))
    return run_OFT(exe, mpi_ranks(nbase, nlevels), timeout)


def validate_sound(nerr_exp,terr_exp,verr_exp,steps_exp=11,linear=False,two_temp=False):
    '''Density, temperature and velocity errors of a sound wave run, and its plot output'''
    errors_ok = validate_errors('sound.results', (('Density', nerr_exp), ('Temperature', terr_exp), ('Velocity', verr_exp)), 1.01)
    output_ok = validate_mug_output(steps_exp, linear=linear, two_temp=two_temp)
    return errors_ok and output_ok


def alfven_setup(representation,nbase,nlevels,order,minlev,nu_mhd='-1',linear=False,mf=False,petsc=False,hex_mesh=False):
    '''Write the input for a traveling shear Alfven wave and run test_alfven (`representation` 'hcurl') or
    test_alfven_lag ('lagrange')

    `nu_mhd` gives the multigrid smoother iterations per level; `mf` selects the matrix-free Newton-Krylov solver
    (nonlinear runs only). Order 4 needs a tighter linear tolerance to reach the expected error.
    '''
    exe, rep_options, timeout = {'hcurl': ('./test_alfven', " bbc='bc'\n", 4000),
                                 'lagrange': ('./test_alfven_lag', "", 2000)}[representation]
    lin_tol = '1.E-12'
    mfnk = False
    if not linear:
        if order == 4:
            lin_tol = '1.E-13'
        mfnk = mf
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(ALFVEN_TEMPLATE.format(nbase=nbase, nlevels=nlevels, order=order, minlev=minlev, nu_mhd=nu_mhd,
                                         linear=namelist_flag(linear), dt='4.E-7', nsteps='250', lin_tol=lin_tol,
                                         mfnk=namelist_flag(mfnk), petsc=namelist_flag(petsc),
                                         mesh_type=(2 if hex_mesh else 1), rep_options=rep_options))
    return run_OFT(exe, mpi_ranks(nbase, nlevels), timeout)


def validate_alfven(berr_exp,verr_exp,steps_exp=11,linear=False):
    '''Magnetic field and velocity errors of an Alfven wave run, and its plot output'''
    errors_ok = validate_errors('alfven.results', (('Magnetic', berr_exp), ('Velocity', verr_exp)), 1.01)
    output_ok = validate_mug_output(steps_exp, linear=linear)
    return errors_ok and output_ok

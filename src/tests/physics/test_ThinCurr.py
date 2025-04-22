from __future__ import print_function
import os
import sys
import time
import multiprocessing
import pytest
import numpy as np
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..','..','python')))
from oft_testing import run_OFT
from OpenFUSIONToolkit.io import histfile
from OpenFUSIONToolkit.util import oftpy_dump_cov

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
 cg_tol={10}
 save_L=T
 save_Mcoil=F
 plot_run=F
 jumper_start={11}
/

&thincurr_eig_options
 direct={1}
 plot_run=F
 neigs={6}
 reduce_model={9}
 jumper_start={11}
/

&thincurr_fr_options
 direct={1}
 freq={2}
 fr_limit={3}
 jumper_start={11}
/

&thincurr_hodlr_options
 target_size=200
 L_svd_tol={7}
 L_aca_rel_tol=0.05
 B_svd_tol={8}
 B_aca_rel_tol=0.05
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


def mp_run(target,args,timeout=180):
    if os.environ.get('OFT_DEBUG_TEST', 0):
        timeout *= 4
    os.chdir(test_dir)
    mp_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=target, args=args + (mp_q,))
    p.start()
    start = time.time()
    while time.time() - start <= timeout:
        if not p.is_alive():
            break
        time.sleep(.5)
    else: # Reached timeout
        print("Timeout reached")
        p.terminate()
        p.join()
        return False
    # Completed successfully?
    try:
        test_result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return False
    p.join()
    return test_result


def run_td(meshfile,direct_flag,use_aca,floops,curr_waveform,volt_waveform,lin_tol,jumper_start,run_reduced,mp_q):
    result = True
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        if meshfile is None:
            from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy
            r_dummy, lc_dummy = build_ThinCurr_dummy([0.0,0.0,10.0],size=0.25)
            tw_model.setup_model(r=r_dummy,lc=lc_dummy,xml_filename='oft_in.xml',jumper_start=jumper_start)
            tw_model.set_eta_values(np.r_[1.E4*mu0])
        else:
            tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
        tw_model.setup_io()
        if floops is not None:
            from OpenFUSIONToolkit.ThinCurr.sensor import circular_flux_loop, save_sensors
            sensors = []
            for (k, floop) in enumerate(floops):
                sensors.append(circular_flux_loop(floop[0],floop[1],'FLOOP_{0}'.format(k),npts=180))
            save_sensors(sensors)
            _, _, sensor_obj = tw_model.compute_Msensor('floops.loc')
        else:
            sensor_obj = None
        if curr_waveform is not None:
            curr_waveform = np.array(curr_waveform)
            curr_waveform[:,1:] /= mu0
        if volt_waveform is not None:
            volt_waveform = np.array(volt_waveform)
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat(use_hodlr=use_aca,cache_file='Lmat.save')
        tw_model.compute_Rmat()
        tw_model.run_td(2.E-5,200,direct=(direct_flag == 'T'),lin_tol=lin_tol,coil_currs=curr_waveform,coil_volts=volt_waveform,sensor_obj=sensor_obj)
        tw_model.plot_td(200,sensor_obj=sensor_obj)
        if run_reduced:
            eig_vals, eig_vecs = tw_model.get_eigs(30,direct=(direct_flag == 'T'))
            tw_reduced = tw_model.build_reduced_model(eig_vecs, compute_B=True, sensor_obj=sensor_obj, filename='tCurr_reduced.h5')
            eig_vals_r, _ = tw_reduced.get_eigs()
            if np.linalg.norm(eig_vals_r[:10]-eig_vals[:10])/np.linalg.norm(eig_vals[:10]) > 1.E-2:
                print("Reduced model eigenvalue error too high ",np.linalg.norm(eig_vals_r[:10]-eig_vals[:10])/np.linalg.norm(eig_vals[:10]))
                result = False
            sensors, currents = tw_reduced.run_td(2.E-5,200,coil_currs=curr_waveform)
            Jreduced = tw_model.reconstruct_current(tw_reduced.reconstruct_potential(currents['curr'][10]),centering='vertex')
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


def run_eig(meshfile,direct_flag,use_aca,jumper_start,mp_q):
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        if meshfile is None:
            from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy
            r_dummy, lc_dummy = build_ThinCurr_dummy([0.0,0.0,10.0],size=0.25,nsplit=1)
            tw_model.setup_model(r=r_dummy,lc=lc_dummy,xml_filename='oft_in.xml',jumper_start=jumper_start)
            tw_model.set_eta_values(np.r_[1.E4*mu0])
        else:
            tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
        tw_model.setup_io()
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat(use_hodlr=use_aca,cache_file='Lmat.save')
        tw_model.compute_Rmat()
        eig_vals, _ = tw_model.get_eigs(4,direct=(direct_flag == 'T'))
        eig_file = '\n'.join(['{0} 0.0'.format(eig_val) for eig_val in eig_vals])
        with open('thincurr_eigs.dat','w+') as fid:
            fid.write(eig_file)
        result = True
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


def run_fr(meshfile,direct_flag,use_aca,freq,fr_limit,floops,jumper_start,mp_q):
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        if meshfile is None:
            from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy
            r_dummy, lc_dummy = build_ThinCurr_dummy([0.0,0.0,10.0],size=0.25,nsplit=1)
            tw_model.setup_model(r=r_dummy,lc=lc_dummy,xml_filename='oft_in.xml',jumper_start=jumper_start)
            tw_model.set_eta_values(np.r_[1.E4*mu0])
        else:
            tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
        tw_model.setup_io()
        if floops is not None:
            from OpenFUSIONToolkit.ThinCurr.sensor import circular_flux_loop, save_sensors
            sensors = []
            for (k, floop) in enumerate(floops):
                sensors.append(circular_flux_loop(floop[0],floop[1],'FLOOP_{0}'.format(k),npts=180))
            save_sensors(sensors)
            Msensor, Msc, _ = tw_model.compute_Msensor('floops.loc')
        Mcoil = tw_model.compute_Mcoil()
        tw_model.compute_Lmat(use_hodlr=use_aca,cache_file='Lmat.save')
        tw_model.compute_Rmat()
        driver_current = 1.0/mu0 # Current is 1.0/mu0 [A]
        driver = np.zeros((2,tw_model.nelems))
        driver[0,:] = Mcoil[0,:]*driver_current
        result = tw_model.compute_freq_response(driver,fr_limit=fr_limit,freq=freq,direct=(direct_flag == 'T'))
        probe_signals = np.dot(result,Msensor)
        probe_signals[0,:] += np.dot(np.r_[driver_current],Msc)
        fr_file = '\n'.join(['{0} {1}'.format(*probe_signals[:,i]) for i in range(probe_signals.shape[1])])
        with open('thincurr_fr.dat','w+') as fid:
            fid.write(fr_file)
        result = True
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


def run_mode(meshfile,freq,mp_q):
    result = True
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        from OpenFUSIONToolkit.ThinCurr.meshing import build_torus_bnorm_grid, ThinCurr_periodic_toroid
        myOFT = OFT_env(nthreads=-1)
        ntheta = 20
        nphi = 20
        r_grid, bnorm, nfp = build_torus_bnorm_grid('tCurr_test_mode.dat',ntheta,nphi,resample_type='theta',use_spline=False)
        plasma_mode = ThinCurr_periodic_toroid(r_grid,nfp,ntheta,nphi)
        plasma_mode.write_to_file('thincurr_mode.h5')
        #
        tw_mode = ThinCurr(myOFT)
        tw_mode.setup_model(mesh_file='thincurr_mode.h5')
        tw_mode.setup_io(basepath='plasma/')
        tw_mode.compute_Lmat()
        Lmat_new = plasma_mode.condense_matrix(tw_mode.Lmat)
        Linv = np.linalg.inv(Lmat_new)
        #
        bnorm_flat = bnorm.reshape((2,bnorm.shape[1]*bnorm.shape[2]))
        flux_flat = bnorm_flat.copy()
        flux_flat[0,plasma_mode.r_map] = tw_mode.scale_va(bnorm_flat[0,plasma_mode.r_map])
        flux_flat[1,plasma_mode.r_map] = tw_mode.scale_va(bnorm_flat[1,plasma_mode.r_map])
        mode_drive = np.zeros((2,tw_mode.nelems))
        for j in range(2):
            output_unique = np.dot(Linv,plasma_mode.nodes_to_unique(flux_flat[j,:]))
            mode_drive[j,:] = plasma_mode.expand_vector(output_unique)
        with h5py.File('thincurr_mode.h5', 'r+') as h5_file:
            h5_file.create_dataset('thincurr/driver', data=mode_drive, dtype='f8')
        #
        tw_torus = ThinCurr(myOFT)
        tw_torus.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml')
        tw_torus.setup_io(basepath='wall/')
        tw_torus.compute_Lmat()
        tw_torus.compute_Rmat()
        mode_driver = tw_mode.cross_eval(tw_torus,mode_drive)
        fr_result = tw_torus.compute_freq_response(fdriver=mode_driver,freq=freq)
        with open('thincurr_mode.dat','w+') as fid:
            fid.write('{0} {1}\n'.format(*np.linalg.norm(mode_drive,axis=1)))
            fid.write('{0} {1}\n'.format(*np.linalg.norm(fr_result,axis=1)))
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)
    

def ThinCurr_setup(meshfile,run_type,direct_flag,freq=0.0,fr_limit=0,eta=10.0,use_aca=False,
                    icoils=None,vcoils=None,floops=None,curr_waveform=None,volt_waveform=None,
                    python=False,lin_tol=1.E-9,jumper_start=0,run_reduced=False):
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
    if use_aca:
        L_svd_tol = 1.E-8
        B_svd_tol = 1.E-3
    else:
        L_svd_tol = -1.0
        B_svd_tol = -1.0
    oft_in_meshfile = meshfile
    if meshfile is None:
        oft_in_meshfile = "tw_test-passive.h5"
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(
            oft_in_meshfile,direct_flag,freq,fr_limit,coil_file_line,volt_file_line,
            neigs,L_svd_tol,B_svd_tol,reduce_model_flag,lin_tol,jumper_start
        ))
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
        coil_string += "</coil_set></icoils>\n"
    if vcoils is not None:
        coil_string += '<vcoils>\n'
        with h5py.File('test_vcoils.h5','w') as h5_file:
            for j, pcoil in enumerate(vcoils):
                coil_string += '<coil_set type="2" res_per_len="1.256637E-5" radius="1.E-2"><coil path="test_vcoils.h5:vcoil_{0:04d}"/></coil_set>\n'.format(j+1)
                R = pcoil[0]; Z = pcoil[1]
                vcoil_pts = []
                for i in range(nPhi):
                    phi = i*phi_fac
                    vcoil_pts.append([R*np.cos(phi), R*np.sin(phi), Z])
                h5_file.create_dataset('vcoil_{0:04d}'.format(j+1), data=np.array(vcoil_pts), dtype='f8')
        coil_string += "</vcoils>\n"
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
        if python:
            return mp_run(run_td,(meshfile,direct_flag,use_aca,floops,curr_waveform,volt_waveform,lin_tol,jumper_start,run_reduced))
        else:
            return run_OFT("../../bin/thincurr_td oft.in oft_in.xml", 1, 180)
    elif run_type == 2:
        if python:
            return mp_run(run_eig,(meshfile,direct_flag,use_aca,jumper_start))
        else:
            return run_OFT("../../bin/thincurr_eig oft.in oft_in.xml", 1, 180)
    elif run_type == 3:
        if python:
            return mp_run(run_fr,(meshfile,direct_flag,use_aca,freq,fr_limit,floops,jumper_start))
        else:
            return run_OFT("../../bin/thincurr_fr oft.in oft_in.xml", 1, 180)
    elif run_type == 5:
        return mp_run(run_mode,(meshfile,freq))

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
    if len(eigs_run_real) < len(eigs):
        print("FAILED: Number of eigenvalues is not correct")
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

def validate_fr(fr_real, fr_imag, tols=(1.E-4, 1.E-4),python=False):
    """
    Helper function to validate frequency-response results against test case.
    """
    try:
        if python:
            fr_run_real = []
            fr_run_imag = []
            with open('thincurr_fr.dat', 'r') as fid:
                for line in fid:
                    vals = line.split()
                    fr_run_real.append(float(vals[0]))
                    fr_run_imag.append(float(vals[1]))
        else:
            hist_file = histfile('thincurr_fr.hist')
            fr_run_real = [hist_file[field][0] for field in hist_file]
            fr_run_imag = [hist_file[field][1] for field in hist_file]
    except BaseException as e:
        print(e)
        return False
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

def validate_td(sigs_final, jumpers_final=None, tols=(1.E-8, 1.E-3)):
    """
    Helper function to validate time-dependent results against test case.
    """
    try:
        hist_file = histfile('floops.hist')
        td_sigs_final = [hist_file[field][-1] for field in hist_file]
    except BaseException as e:
        print(e)
        return False
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
    #
    if jumpers_final is not None:
        try:
            hist_file = histfile('jumpers.hist')
            td_sigs_final = [hist_file[field][-1] for field in hist_file]
        except BaseException as e:
            print(e)
            return False
        if not len(td_sigs_final) == len(jumpers_final):
            print("FAILED: Number of jumpers does not match")
            return False
        retval = True
        if abs(jumpers_final[0]-td_sigs_final[0]) > tols[0]:
            print("FAILED: Final time incorrect!")
            print("  Expected = {0}".format(jumpers_final[0]))
            print("  Actual =   {0}".format(td_sigs_final[0]))
            retval = False
        for (i, val) in enumerate(jumpers_final[1:]):
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


def validate_mode(drive_exp,result_exp):
    try:
        with open('thincurr_mode.dat','r') as fid:
            drive_amps = [float(val) for val in fid.readline().split()]
            result_amps = [float(val) for val in fid.readline().split()]
    except BaseException as e:
        print(e)
        return False
    result_val = True
    for i in range(2):
        if abs((drive_exp[i]-drive_amps[i])/drive_exp[i]) > 1.E-4:
            print("FAILED: driver {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(drive_exp[i]))
            print("  Actual =   {0}".format(drive_amps[i]))
            result_val = False
        if abs((result_exp[i]-result_amps[i])/result_exp[i]) > 1.E-4:
            print("FAILED: result {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(result_exp[i]))
            print("  Actual =   {0}".format(result_amps[i]))
            result_val = False
    return result_val

#============================================================================
# Test runners for plate
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_eig_plate(direct_flag,python):
    eigs = (9.735667E-3, 6.532314E-3, 6.532201E-3, 5.251598E-3)
    assert ThinCurr_setup("tw_test-plate.h5",2 if python else 4,direct_flag,python=python)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_plate(direct_flag,python):
    sigs_final = (4.E-3, 8.459371E-4, 7.130923E-4)
    assert ThinCurr_setup("tw_test-plate.h5",1,direct_flag,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           python=python)
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_fr_plate(direct_flag,python):
    fr_real = (6.807649E-2, 7.207748E-2)
    fr_imag = (-3.011666E-3, -2.177010E-3)
    assert ThinCurr_setup("tw_test-plate.h5",3,direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           python=python)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_plate_volt(direct_flag,python):
    sigs_final = (4.E-3, 4.580643E-4, 3.854292E-4)
    assert ThinCurr_setup("tw_test-plate.h5",1,direct_flag,
                           vcoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           volt_waveform=((0.0, 1.0), (1.0, 1.0)),
                           python=python)
    assert validate_td(sigs_final)

#============================================================================
# Test runners for cylinder
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_eig_cyl(direct_flag,python):
    eigs = (2.657195E-2, 1.248071E-2, 1.247103E-2, 1.200566E-2)
    assert ThinCurr_setup("tw_test-cyl.h5",2 if python else 4,direct_flag,python=python,jumper_start=2)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_cyl(direct_flag,python):
    sigs_final = (4.E-3, 7.254196E-4, 6.151460E-4)
    jumpers_final = (4.E-3, 5.445469E3)
    assert ThinCurr_setup("tw_test-cyl.h5",1,direct_flag,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           python=python,jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_fr_cyl(direct_flag,python):
    fr_real = (6.118337E-2, 4.356188E-3)
    fr_imag = (-1.911861E-3, -2.283493E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",3,direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           python=python,jumper_start=2)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_cyl_volt(direct_flag,python):
    sigs_final = (4.E-3, 1.504279E-4, 1.276624E-4)
    jumpers_final = (4.E-3, 1.1203960E3)
    assert ThinCurr_setup("tw_test-cyl.h5",1,direct_flag,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           python=python,jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

#============================================================================
# Test runners for torus
@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_eig_torus(direct_flag,python):
    eigs = (4.751344E-2, 2.564491E-2, 2.555695E-2, 2.285850E-2)
    assert ThinCurr_setup("tw_test-torus.h5",2 if python else 4,direct_flag,python=python)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_torus(direct_flag,python):
    sigs_final = (4.E-3, 4.935683E-4, 3.729159E-5)
    assert ThinCurr_setup("tw_test-torus.h5",1,direct_flag,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           lin_tol=1.E-10,
                           python=python)
    assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_fr_torus(direct_flag,python):
    fr_real = (-2.807955E-3, -1.196091E-4)
    fr_imag = (-1.869732E-3, -1.248642E-4)
    assert ThinCurr_setup("tw_test-torus.h5",3,direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           python=python)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.coverage
def test_mode_torus():
    drive_exp = (73.91361257364348, 48.430633246949554)
    result_exp = (58.713811707231145, 41.5238351334917)
    assert ThinCurr_setup("tw_test-torus.h5",5,False,freq=1.E3)
    assert validate_mode(drive_exp,result_exp)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_torus_volt(direct_flag,python):
    sigs_final = (4.E-3, 5.653338E-5, 4.035387E-6)
    assert ThinCurr_setup("tw_test-torus.h5",1,direct_flag,
                           vcoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           lin_tol=1.E-11,
                           python=python)
    assert validate_td(sigs_final)

#============================================================================
# Test runners for filament model
@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_eig_passive(direct_flag,python):
    eigs = (1.504155E-1, 6.423383E-2, 3.190175E-2, 2.942398E-2)
    assert ThinCurr_setup(None,2 if python else 4,direct_flag,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)),python=python)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_passive(direct_flag,python):
   sigs_final = (4.E-3, 8.349309E-4, 8.364054E-4)
   assert ThinCurr_setup(None,1,direct_flag,eta=1.E4,
                          icoils=((0.5, 0.1),),
                          vcoils=((0.5, 0.0),),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                          python=python)
   assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_fr_passive(direct_flag,python):
    fr_real = (1.947713E-1, 1.990873E-1)
    fr_imag = (-2.174952E-4, -1.560016E-4)
    assert ThinCurr_setup(None,3,direct_flag,eta=1.E4,freq=5.E3,fr_limit=0,
                           icoils=((0.5, 0.1),),
                           vcoils=((0.5, 0.0),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           python=python)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (False, True))
def test_td_passive_volt(direct_flag,python):
   sigs_final = (4.E-3, 4.379235E-4, 4.389248E-4)
   assert ThinCurr_setup(None,1,direct_flag,eta=1.E4,
                          vcoils=((0.5, 0.0), (0.5, 0.1)),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          volt_waveform=((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)),
                          python=python)
   assert validate_td(sigs_final)

#============================================================================
# Test runners for large cylinder (w/ ACA+)
@pytest.mark.coverage
@pytest.mark.parametrize("python", (False, True))
def test_eig_aca(python):
    eigs = (2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",2 if python else 4,'F',use_aca=True,python=python,jumper_start=2)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("python", (False, True))
def test_td_aca(python):
    eigs = (2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2)
    sigs_final = (4.E-3, 7.280671E-4, 6.211245E-4)
    jumpers_final = (4.E-3, 5.447048E3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",1,'F',use_aca=True,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           python=python,jumper_start=2,run_reduced=True)
    if python:
        assert validate_model_red(eigs)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
@pytest.mark.parametrize("python", (False, True))
def test_fr_aca(python):
    fr_real = (5.888736E-2, 4.881440E-3)
    fr_imag = (-2.017045E-3, -2.313881E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",3,'F',use_aca=True,freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           python=python,jumper_start=2)
    assert validate_fr(fr_real, fr_imag, tols=(1.E-3, 1.E-3))

@pytest.mark.coverage
@pytest.mark.parametrize("python", (False, True))
def test_td_volt_aca(python):
    sigs_final = (4.E-3, 1.512679E-4, 1.291681E-4)
    jumpers_final = (4.E-3, 1.122550E3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",1,'F',use_aca=True,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           python=python,jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)
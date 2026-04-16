from __future__ import print_function
import os
import sys
import time
import glob
import shutil
import multiprocessing
import pytest
import numpy as np
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..','..','python')))
from oft_testing import run_OFT
from OpenFUSIONToolkit.io import histfile, write_oft_xml
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_Vcoil, ThinCurr_XML



mu0 = np.pi*4.E-7
_oft_env_singleton = None

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
            tw_model.set_eta_values(eta_values=np.r_[1.E4*mu0])
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
            tw_model.set_eta_values(eta_values=np.r_[1.E4*mu0])
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
            tw_model.set_eta_values(eta_values=np.r_[1.E4*mu0])
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
    

def run_td_for_Mirnov(meshfile,direct_flag,curr_waveform,lin_tol,mp_q):
    result = True
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in_remc.xml')
        tw_model.setup_io()
        _, _, sensor_obj = tw_model.compute_Msensor('floops.loc')
        if curr_waveform is not None:
            curr_waveform = np.array(curr_waveform)
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat()
        tw_model.compute_Rmat()
        tw_model.run_td(2.E-4,200,direct=(direct_flag == 'T'),lin_tol=lin_tol,coil_currs=curr_waveform,sensor_obj=sensor_obj)
        tw_model.plot_td(200,sensor_obj=sensor_obj)
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
    xml = ThinCurr_XML()
    xml.set_eta([eta*mu0])
    if icoils is not None:
        test_Icoil = ThinCurr_Icoil(name="Test")
        for icoil in icoils:
            R = icoil[0]; Z = icoil[1]
            coil_pts = []
            for i in range(nPhi):
                phi = i*phi_fac
                coil_pts.append([R*np.cos(phi), R*np.sin(phi), Z])
            test_Icoil.add_subcoil(pts=coil_pts)
        xml.add_Icoil(test_Icoil)
    if vcoils is not None:
        with h5py.File('test_vcoils.h5','w') as h5_file:
            for j, vcoil in enumerate(vcoils):
                test_Vcoil = ThinCurr_Vcoil(name="Test_{0:d}".format(j+1), resistivity_per_length=1.256637E-5, radius=1.E-2)
                R = vcoil[0]; Z = vcoil[1]
                vcoil_pts = []
                for i in range(nPhi):
                    phi = i*phi_fac
                    vcoil_pts.append([R*np.cos(phi), R*np.sin(phi), Z])
                h5_file.create_dataset('vcoil_{0:04d}'.format(j+1), data=np.array(vcoil_pts), dtype='f8')
                test_Vcoil.add_subcoil(hdf5_path='test_vcoils.h5:vcoil_{0:04d}'.format(j+1))
                xml.add_Vcoil(test_Vcoil)
    write_oft_xml([xml], "oft_in.xml")
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
    elif run_type == 6:
        return mp_run(run_td_for_Mirnov,(meshfile,direct_flag,curr_waveform,lin_tol))


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
        for (i, val) in enumerate(td_sigs_final[1:]):
            print(val)
        if len(td_sigs_final) != len(jumpers_final):
            print("FAILED: Number of jumpers does not match")
            return False
        retval = True
        if abs(jumpers_final[0]-td_sigs_final[0]) > tols[0]:
            print("FAILED: Final time incorrect!")
            print("  Expected = {0}".format(jumpers_final[0]))
            print("  Actual =   {0}".format(td_sigs_final[0]))
            retval = False
        for (i, val) in enumerate(jumpers_final[1:]):
            if val is None:
                continue
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


def _write_thincurr_xml(xml_filename, eta_values=None, thickness_values=None, eta_vol_values=None):
    from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_XML

    thincurr_xml = ThinCurr_XML()
    if eta_values is not None:
        thincurr_xml.set_eta(eta_values)
    if eta_vol_values is not None:
        thincurr_xml.set_eta_vol(eta_vol_values)
    if thickness_values is not None:
        thincurr_xml.set_thickness(thickness_values)
    write_oft_xml([thincurr_xml], xml_filename)


def _build_dummy_model(xml_filename):
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.ThinCurr import ThinCurr
    from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy

    global _oft_env_singleton
    if _oft_env_singleton is None:
        _oft_env_singleton = OFT_env(nthreads=-1)

    r_dummy, lc_dummy = build_ThinCurr_dummy([0.0, 0.0, 10.0], size=0.25)
    model = ThinCurr(_oft_env_singleton)
    model.setup_model(r=r_dummy, lc=lc_dummy, xml_filename=xml_filename)
    return model


def validate_torus_fourier_sensor(interface,sigs_nmodes_1D_PEST,sigs_nmodes_1D_Hamada,sigs_mnmodes_2D_PEST,sigs_mnmodes_2D_Hamada,t,delta_phi,tol=1.E-6):
    try:
        interface.load_histfile()
    except BaseException as e:
        print(e)
        return False
    
    result_val = True
    import matplotlib.pyplot as plt
    _,ax=plt.subplots(1,1,figsize=(8,6))
    data = interface.plot_1D_fourier_amplitude(t,1,ax,toroidal_harmonics=True,hamada_dphi=None,part='r')[1]
    if np.linalg.norm(abs(sigs_nmodes_1D_PEST-data),np.inf)>tol:
        print(f"FAILED: 1D PEST toroidal Fourier transform at t = {t} and helicity = {interface.helicity} incorrect!")
        np.save('sigs_nmodes_1D_PEST-new.npy',data)
        result_val = False
    data = interface.plot_1D_fourier_amplitude(t,1,ax,toroidal_harmonics=True,hamada_dphi=delta_phi,part='r')[1]
    if np.linalg.norm(abs(sigs_nmodes_1D_Hamada-data),np.inf)>tol:
        print(f"FAILED: 1D Hamada toroidal Fourier transform at t = {t} and helicity = {interface.helicity} incorrect!")
        np.save('sigs_nmodes_1D_Hamada-new.npy',data)
        result_val = False
    data = interface.fft2(interface.get_B_mesh(t),hamada_dphi=None)[0]
    if np.linalg.norm(abs(sigs_mnmodes_2D_PEST-data),np.inf)>tol:
        print(f"FAILED: 2D PEST Fourier transform at t = {t} and helicity = {interface.helicity} incorrect!")
        np.save('sigs_mnmodes_2D_PEST-new.npy',data)
        result_val = False
    data = interface.fft2(interface.get_B_mesh(t),hamada_dphi=delta_phi)[0]
    if np.linalg.norm(abs(sigs_mnmodes_2D_Hamada-data),np.inf)>tol:
        print(f"FAILED: 2D Hamada Fourier transform at t = {t} and helicity = {interface.helicity} incorrect!")
        np.save('sigs_mnmodes_2D_Hamada-new.npy',data)
        result_val = False
    run_files = [f for f in os.listdir('.') if f.endswith('.rst') or f.endswith('.xmf') or f.endswith('.loc')]
    for file in run_files:
        os.remove(file)
    save_files = [f for f in os.listdir('.') if f.startswith('mesh') or f.startswith('vector') or f.startswith('scalar') or f.startswith('dump.dat')]
    for file in save_files:
        os.remove(file)  
    return result_val


@pytest.mark.coverage
def test_thickness_api_roundtrip_and_validation():
    assert mp_run(run_thickness_api_roundtrip_and_validation, ())


def run_thickness_api_roundtrip_and_validation(mp_q):
    result = True
    try:
        import warnings
        os.chdir(test_dir)
        xml_filename = 'oft_in_thickness_api.xml'
        eta_values = np.r_[1.E4*mu0]
        thickness_values = np.r_[2.5E-3]
        eta_vol_values = eta_values * thickness_values
        _write_thincurr_xml(xml_filename, eta_values=eta_values, thickness_values=thickness_values)

        tw_model = _build_dummy_model(xml_filename)

        tw_model.set_eta_values(eta_surf=eta_values, thickness=thickness_values)
        
        # Verify get_eta_values() without flag shows deprecation warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            eta_returned = tw_model.get_eta_values()
            if len(w) == 0 or not issubclass(w[-1].category, DeprecationWarning):
                print("WARNING: Expected DeprecationWarning for get_eta_values() not found")
        
        if not np.allclose(eta_returned, eta_values):
            result = False
        eta_surf_out, eta_vol_out = tw_model.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_out, eta_values):
            result = False
        if not np.allclose(eta_vol_out, eta_vol_values):
            result = False
        if not np.allclose(tw_model.get_thickness(), thickness_values):
            result = False

        tw_model_pair = _build_dummy_model(xml_filename)
        tw_model_pair.set_eta_values(eta_surf=eta_values, eta_vol=eta_vol_values)
        eta_surf_pair, eta_vol_pair = tw_model_pair.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_pair, eta_values):
            result = False
        if not np.allclose(eta_vol_pair, eta_vol_values):
            result = False
        if not np.allclose(tw_model_pair.get_thickness(), thickness_values):
            result = False

        tw_model_alias = _build_dummy_model(xml_filename)
        tw_model_alias.set_eta_values(eta_values=eta_values, thickness=thickness_values)
        eta_surf_alias, eta_vol_alias = tw_model_alias.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_alias, eta_values):
            result = False
        if not np.allclose(eta_vol_alias, eta_vol_values):
            result = False

        tw_model_eta_only = _build_dummy_model(xml_filename)
        tw_model_eta_only.set_eta_values(eta_surf=eta_values)
        eta_surf_only, eta_vol_only = tw_model_eta_only.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_only, eta_values):
            result = False
        if eta_vol_only is not None:
            result = False
        if not np.allclose(tw_model_eta_only.get_thickness(), -1.0):
            result = False

        eta_surf_override = eta_values * 1.5
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            tw_model_three = _build_dummy_model(xml_filename)
            tw_model_three.set_eta_values(eta_surf=eta_surf_override, eta_vol=eta_vol_values, thickness=thickness_values)
        if len(w) == 0 or not any(issubclass(item.category, UserWarning) for item in w):
            result = False
        eta_surf_three, eta_vol_three = tw_model_three.get_eta_values(include_eta_vol=True)
        if not np.allclose(eta_surf_three, eta_values):
            result = False
        if not np.allclose(eta_vol_three, eta_vol_values):
            result = False
        if not np.allclose(tw_model_three.get_thickness(), thickness_values):
            result = False

        # Test updating thickness with a new set_eta_values call
        new_thickness = np.r_[3.5E-3]
        tw_model.set_eta_values(eta_surf=eta_values, thickness=new_thickness)
        if not np.allclose(tw_model.get_thickness(), new_thickness):
            result = False

        with pytest.raises(ValueError):
            tw_model.set_eta_values()
        with pytest.raises(ValueError):
            tw_model.set_eta_values(eta_values=np.r_[-1.0])
        with pytest.raises(ValueError):
            tw_model.set_eta_values(thickness=thickness_values)
        with pytest.raises(ValueError):
            tw_model.set_eta_values(eta_surf=eta_values, thickness=np.r_[-1.0])
        # eta_vol requires pre-existing thickness when thickness is not passed.
        # Use a fresh model with no thickness XML to exercise this error path.
        xml_no_thickness = 'oft_in_no_thickness.xml'
        _write_thincurr_xml(xml_no_thickness, eta_values=eta_values)
        tw_model_no_thickness = _build_dummy_model(xml_no_thickness)
        with pytest.raises(Exception):
            tw_model_no_thickness.set_eta_values(eta_vol=np.r_[1.0e-6])
        with pytest.raises(IndexError):
            tw_model.set_eta_values(eta_surf=eta_values, thickness=np.r_[1.E-3, 2.E-3])
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_plot_td_compute_jvol_outputs_fields():
    assert mp_run(run_plot_td_compute_jvol_outputs_fields, ())


def run_plot_td_compute_jvol_outputs_fields(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_compute_jvol.xml'
        io_basepath = 'td_jvol_regression'
        eta_values = np.r_[1.E4*mu0]
        thickness_values = np.r_[2.0E-3]
        _write_thincurr_xml(xml_filename, eta_values=eta_values, thickness_values=thickness_values)

        if os.path.isdir(io_basepath):
            shutil.rmtree(io_basepath)

        tw_model = _build_dummy_model(xml_filename)
        tw_model.setup_io(basepath=io_basepath)
        tw_model.compute_Mcoil()
        tw_model.compute_Lmat()
        tw_model.compute_Rmat()
        tw_model.run_td(2.E-5, 5, direct=True, plot_freq=1)
        tw_model.plot_td(5, plot_freq=1)
        tw_model.build_XDMF()

        xmf_files = sorted(glob.glob(os.path.join(io_basepath, '*.xmf')))
        if len(xmf_files) == 0:
            result = False
        xmf_text = ''
        for xmf_file in xmf_files:
            with open(xmf_file, 'r') as fid:
                xmf_text += fid.read()
        if 'J_vol' not in xmf_text or 'thickness' not in xmf_text:
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_eta_only_matches_surface_resistivity_with_thickness():
    assert mp_run(run_eta_only_matches_surface_resistivity_with_thickness, ())


def run_eta_only_matches_surface_resistivity_with_thickness(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_compat.xml'
        eta_surface = np.r_[1.E4*mu0]
        thickness_values = np.r_[2.5E-3]
        eta_bulk = eta_surface*thickness_values
        _write_thincurr_xml(xml_filename, eta_values=eta_surface)

        model_surface = _build_dummy_model(xml_filename)
        model_surface.set_eta_values(eta_values=eta_surface)
        model_surface.compute_Rmat(copy_out=True)
        R_surface = model_surface.Rmat

        model_bulk = _build_dummy_model(xml_filename)
        model_bulk.set_eta_values(eta_vol=eta_bulk, thickness=thickness_values)
        model_bulk.compute_Rmat(copy_out=True)
        R_bulk = model_bulk.Rmat

        if not np.allclose(R_surface, R_bulk, rtol=1.E-10, atol=1.E-12):
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_eta_vol_without_thickness_warns_and_loads():
    assert mp_run(run_eta_vol_without_thickness_warns_and_loads, ())


def run_eta_vol_without_thickness_warns_and_loads(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_vol_only_invalid.xml'
        eta_vol_values = np.r_[2.5E-3]
        _write_thincurr_xml(xml_filename, eta_vol_values=eta_vol_values)

        # eta_vol-only XML should load with warning; eta_surf cannot be inferred.
        tw_model = _build_dummy_model(xml_filename)
        eta_surf = tw_model.get_eta_values()
        if np.all(eta_surf > 0.0):
            result = False
    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_get_eta_values_deprecation_warning():
    """Test that get_eta_values() without include_eta_vol flag shows deprecation warning."""
    assert mp_run(run_get_eta_values_deprecation_warning, ())


def run_get_eta_values_deprecation_warning(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_deprecation.xml'
        eta_values = np.r_[1.E4*mu0]
        _write_thincurr_xml(xml_filename, eta_values=eta_values)

        tw_model = _build_dummy_model(xml_filename)
        tw_model.set_eta_values(eta_surf=eta_values)

        import warnings
        # Check that calling get_eta_values() without flag triggers deprecation warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            eta = tw_model.get_eta_values()
            
            # Should get a DeprecationWarning
            if len(w) == 0 or not issubclass(w[-1].category, DeprecationWarning):
                print("ERROR: Expected DeprecationWarning not found")
                result = False
            elif "include_eta_vol" not in str(w[-1].message):
                print(f"ERROR: Expected 'include_eta_vol' in warning message, got: {w[-1].message}")
                result = False
            elif not np.allclose(eta, eta_values):
                print("ERROR: eta_surf values don't match")
                result = False

    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)


@pytest.mark.coverage
def test_get_eta_values_with_flag_returns_tuple():
    """Test that get_eta_values(include_eta_vol=True) returns tuple with eta_vol."""
    assert mp_run(run_get_eta_values_with_flag_returns_tuple, ())


def run_get_eta_values_with_flag_returns_tuple(mp_q):
    result = True
    try:
        os.chdir(test_dir)
        xml_filename = 'oft_in_eta_vol_tuple.xml'
        eta_surf = np.r_[1.E4*mu0]
        eta_vol = np.r_[2.5E-2]
        thickness = np.r_[2.5E-3]
        _write_thincurr_xml(xml_filename, eta_values=eta_surf, thickness_values=thickness)

        tw_model = _build_dummy_model(xml_filename)
        tw_model.set_eta_values(eta_vol=eta_vol, thickness=thickness)

        import warnings
        # Call with include_eta_vol=True should NOT trigger warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result_tuple = tw_model.get_eta_values(include_eta_vol=True)
            
            # Should not get DeprecationWarning
            for warning in w:
                if issubclass(warning.category, DeprecationWarning):
                    print(f"ERROR: Unexpected DeprecationWarning when include_eta_vol=True: {warning.message}")
                    result = False
                    break

        # Check that result is a tuple with (eta_surf, eta_vol)
        if not isinstance(result_tuple, tuple) or len(result_tuple) != 2:
            print(f"ERROR: Expected tuple of length 2, got {type(result_tuple)} of length {len(result_tuple) if isinstance(result_tuple, tuple) else 'N/A'}")
            result = False
        else:
            eta_surf_ret, eta_vol_ret = result_tuple
            if eta_vol_ret is None:
                print("ERROR: eta_vol should not be None when it was explicitly set")
                result = False
            elif not np.allclose(eta_vol_ret, eta_vol):
                print(f"ERROR: eta_vol values don't match. Expected {eta_vol}, got {eta_vol_ret}")
                result = False

    except BaseException as e:
        print(e)
        result = False
    oftpy_dump_cov()
    mp_q.put(result)








#============================================================================
# Test runners for plate
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_eig_plate(direct_flag,python):
    eigs = (9.735667E-3, 6.532314E-3, 6.532201E-3, 5.251598E-3)
    assert ThinCurr_setup("tw_test-plate.h5",2 if python else 4,direct_flag,python=python)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_td_plate(direct_flag,python):
    sigs_final = (4.E-3, 8.459371E-4, 7.130923E-4)
    assert ThinCurr_setup("tw_test-plate.h5",1,direct_flag,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           python=python)
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_fr_plate(direct_flag,python):
    fr_real = (6.807649E-2, 7.207748E-2)
    fr_imag = (-3.011666E-3, -2.177010E-3)
    assert ThinCurr_setup("tw_test-plate.h5",3,direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           python=python)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_td_plate_volt(direct_flag,python):
    sigs_final = (4.E-3, 4.580643E-4, 3.854292E-4)
    jumpers_final = (4.E-3, 1697.895)
    assert ThinCurr_setup("tw_test-plate.h5",1,direct_flag,
                           vcoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           volt_waveform=((0.0, 1.0), (1.0, 1.0)),
                           python=python)
    assert validate_td(sigs_final,jumpers_final)

#============================================================================
# Test runners for cylinder
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_eig_cyl(direct_flag,python):
    eigs = (2.657195E-2, 1.248071E-2, 1.247103E-2, 1.200566E-2)
    assert ThinCurr_setup("tw_test-cyl.h5",2 if python else 4,direct_flag,python=python,jumper_start=2)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_td_cyl(direct_flag,python):
    sigs_final = (4.E-3, 7.254196E-4, 6.151460E-4)
    jumpers_final = (4.E-3, 5.445469E3, 5445.469)
    assert ThinCurr_setup("tw_test-cyl.h5",1,direct_flag,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           python=python,jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_fr_cyl(direct_flag,python):
    fr_real = (6.118337E-2, 4.356188E-3)
    fr_imag = (-1.911861E-3, -2.283493E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",3,direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           python=python,jumper_start=2)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_td_cyl_volt(direct_flag,python):
    sigs_final = (4.E-3, 1.504279E-4, 1.276624E-4)
    jumpers_final = (4.E-3, 1.1203960E3, 1120.396, 655.853, 655.850)
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
@pytest.mark.parametrize("python", (True,))
def test_eig_torus(direct_flag,python):
    eigs = (4.751344E-2, 2.564491E-2, 2.555695E-2, 2.285850E-2)
    assert ThinCurr_setup("tw_test-torus.h5",2 if python else 4,direct_flag,python=python)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
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
@pytest.mark.parametrize("python", (True,))
def test_fr_torus(direct_flag,python):
    fr_real = (-2.807955E-3, -1.196091E-4)
    fr_imag = (-1.869732E-3, -1.248642E-4)
    assert ThinCurr_setup("tw_test-torus.h5",3,direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           python=python)
    assert validate_fr(fr_real, fr_imag, python=python)

# @pytest.mark.coverage
# def test_mode_torus():
#     drive_exp = (73.91361257364348, 48.430633246949554)
#     result_exp = (58.713811707231145, 41.5238351334917)
#     assert ThinCurr_setup("tw_test-torus.h5",5,False,freq=1.E3)
#     assert validate_mode(drive_exp,result_exp)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_td_torus_volt(direct_flag,python):
    sigs_final = (4.E-3, 5.653338E-5, 4.035387E-6)
    jumpers_final = (4.E-3, None, -597.6068, 371.74769, 371.74780)
    assert ThinCurr_setup("tw_test-torus.h5",1,direct_flag,
                           vcoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           lin_tol=1.E-11,
                           python=python)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_torus_fourier_sensor(direct_flag):
    from OpenFUSIONToolkit.ThinCurr.util import torus_fourier_sensor
    import xarray as xr
    R_0 = 1.0
    ds = xr.open_dataset("torus_gpec_control_output_n1_nc.nc")
    R_gpec=ds.R.to_dataframe().values[:,0][:-1]
    Z_gpec=ds.z.to_dataframe().values[:,0][:-1]
    delta_phi = ds.delta_phi.to_dataframe().values[:,0][:-1]
    interface_h1 = torus_fourier_sensor(R_gpec,Z_gpec,R_0,1)
    interface_h1.place_normal_sensors(nphi=15,filename='floops.loc')
    interface_hminus1 = torus_fourier_sensor(R_gpec,Z_gpec,R_0,-1)
    interface_hminus1.place_normal_sensors(nphi=15,filename='floops.loc')
    t = 80
    assert ThinCurr_setup("tw_test-torus.h5",6,direct_flag,
                           curr_waveform=((0.0, 1.E6), (4.E-3, 0.0), (1.0, 0.0)),
                           lin_tol=1.E-10,
                           python=True)
    sigs_nmodes_1D_PEST = np.load('sigs_nmodes_1D_PEST-h1.npy')
    sigs_nmodes_1D_Hamada = np.load('sigs_nmodes_1D_Hamada-h1.npy')
    sigs_mnmodes_2D_PEST = np.load('sigs_mnmodes_2D_PEST-h1.npy')
    sigs_mnmodes_2D_Hamada = np.load('sigs_mnmodes_2D_Hamada-h1.npy')
    assert validate_torus_fourier_sensor(interface_h1,sigs_nmodes_1D_PEST,sigs_nmodes_1D_Hamada,sigs_mnmodes_2D_PEST,sigs_mnmodes_2D_Hamada,t,delta_phi)
    sigs_nmodes_1D_PEST = np.load('sigs_nmodes_1D_PEST-hminus1.npy')
    sigs_nmodes_1D_Hamada = np.load('sigs_nmodes_1D_Hamada-hminus1.npy')
    sigs_mnmodes_2D_PEST = np.load('sigs_mnmodes_2D_PEST-hminus1.npy')
    sigs_mnmodes_2D_Hamada = np.load('sigs_mnmodes_2D_Hamada-hminus1.npy')
    assert validate_torus_fourier_sensor(interface_hminus1,sigs_nmodes_1D_PEST,sigs_nmodes_1D_Hamada,sigs_mnmodes_2D_PEST,sigs_mnmodes_2D_Hamada,t,delta_phi)

#============================================================================
# Test runners for filament model
@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_eig_passive(direct_flag,python):
    eigs = (1.503561E-1, 6.420533E-2, 3.188782E-2, 2.941118E-2)
    assert ThinCurr_setup(None,2 if python else 4,direct_flag,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)),python=python)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
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
@pytest.mark.parametrize("python", (True,))
def test_fr_passive(direct_flag,python):
    fr_real = (1.947713E-1, 1.990873E-1)
    fr_imag = (-2.175942E-4, -1.560726E-4)
    assert ThinCurr_setup(None,3,direct_flag,eta=1.E4,freq=5.E3,fr_limit=0,
                           icoils=((0.5, 0.1),),
                           vcoils=((0.5, 0.0),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           python=python)
    assert validate_fr(fr_real, fr_imag, python=python)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("python", (True,))
def test_td_passive_volt(direct_flag,python):
   sigs_final = (4.E-3, 4.379235E-4, 4.389248E-4)
   jumpers_final = (4.E-3, -641.4736, 1673.2893)
   assert ThinCurr_setup(None,1,direct_flag,eta=1.E4,
                          vcoils=((0.5, 0.0), (0.5, 0.1)),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          volt_waveform=((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)),
                          python=python)
   assert validate_td(sigs_final, jumpers_final)

#============================================================================
# Test runners for large cylinder (w/ ACA+)
@pytest.mark.coverage
@pytest.mark.parametrize("python", (True,))
def test_eig_aca(python):
    eigs = (2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",2 if python else 4,'F',use_aca=True,python=python,jumper_start=2)
    assert validate_eigs(eigs)
    if not python:
        assert validate_model_red(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("python", (True,))
def test_td_aca(python):
    eigs = (2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2)
    sigs_final = (4.E-3, 7.280671E-4, 6.211245E-4)
    jumpers_final = (4.E-3, 5.447048E3, 5447.048)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",1,'F',use_aca=True,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           python=python,jumper_start=2,run_reduced=True)
    if python:
        assert validate_model_red(eigs)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
@pytest.mark.parametrize("python", (True,))
def test_fr_aca(python):
    fr_real = (5.888736E-2, 4.881440E-3)
    fr_imag = (-2.017045E-3, -2.313881E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",3,'F',use_aca=True,freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           python=python,jumper_start=2)
    assert validate_fr(fr_real, fr_imag, python=python, tols=(1.E-3, 1.E-3))

@pytest.mark.coverage
@pytest.mark.parametrize("python", (True,))
def test_td_volt_aca(python):
    sigs_final = (4.E-3, 1.512679E-4, 1.291681E-4)
    jumpers_final = (4.E-3, 1.122550E3, 1122.550, 656.9544, 656.9797)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",1,'F',use_aca=True,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           python=python,jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)
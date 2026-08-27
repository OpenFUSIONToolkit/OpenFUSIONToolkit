import os
import sys
import time
import glob
import shutil
import pathlib
import subprocess
import multiprocessing
import pytest
import numpy as np
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
python_dir = os.path.abspath(os.path.join(test_dir, '..','..','python'))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
sys.path.append(python_dir)
from oft_testing import approx_list
from OpenFUSIONToolkit.io import histfile, write_oft_xml
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_Vcoil, ThinCurr_XML
from OpenFUSIONToolkit.ThinCurr.sensor import save_sensors, circular_flux_loop



mu0 = np.pi*4.E-7
_oft_env_singleton = None


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


def run_td(meshfile,direct_flag,use_aca,floops,curr_waveform,volt_waveform,lin_tol,jumper_start,run_reduced,basepath,mp_q):
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
            if pathlib.Path("oft_in.xml").is_file():
                tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
            else:
                tw_model.setup_model(mesh_file=meshfile)
        tw_model.setup_io(basepath=basepath)
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
        rst_before = set(glob.glob('pThinCurr_*.rst'))
        hist_before = os.path.exists('floops.hist')
        tw_model.run_td(2.E-5,200,direct=(direct_flag == 'T'),lin_tol=lin_tol,coil_currs=curr_waveform,coil_volts=volt_waveform,sensor_obj=sensor_obj)
        tw_model.plot_td(200,sensor_obj=sensor_obj)
        if basepath is not None:
            # Verify output files landed in basepath, not CWD
            rst_in_outdir = glob.glob(os.path.join(basepath, 'pThinCurr_*.rst'))
            new_rst_in_cwd = set(glob.glob('pThinCurr_*.rst')) - rst_before
            if len(rst_in_outdir) == 0:
                print("FAILED: No .rst files found in output directory")
                result = False
            if len(new_rst_in_cwd) > 0:
                print("FAILED: {0} .rst file(s) written to CWD instead of output directory".format(len(new_rst_in_cwd)))
                result = False
            if floops is not None:
                if not os.path.exists(os.path.join(basepath, 'floops.hist')):
                    print("FAILED: floops.hist not found in output directory")
                    result = False
                if not hist_before and os.path.exists('floops.hist'):
                    print("FAILED: floops.hist written to CWD instead of output directory")
                    result = False
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
    finally:
        if basepath is not None and os.path.isdir(basepath):
            shutil.rmtree(basepath)
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
            if pathlib.Path("oft_in.xml").is_file():
                tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
            else:
                tw_model.setup_model(mesh_file=meshfile)
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
            if pathlib.Path("oft_in.xml").is_file():
                tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
            else:
                tw_model.setup_model(mesh_file=meshfile)
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
        if pathlib.Path("oft_in.xml").is_file():
            tw_torus.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml')
        else:
            tw_torus.setup_model(mesh_file=meshfile)
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
                    lin_tol=1.E-9,jumper_start=0,run_reduced=False,basepath=None,convert_xml=False):
    """
    Common setup and run operations for thin-wall physics module test cases
    """
    nPhi = 180
    phi_fac = 2.0*np.pi/(nPhi-1)
    # Create main input file from template
    os.chdir(test_dir)
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
    if convert_xml:
        mesh_tool_args = ["--in_file={0}".format(meshfile), "--out_file=test.h5", "--eta_from_xml=oft_in.xml"]
        if (icoils is not None) or (vcoils is not None):
            mesh_tool_args += ["--coils_from_xml=oft_in.xml"]
        if jumper_start > 0:
            mesh_tool_args += ["--jumper_range", "{0}".format(jumper_start-1)]
        result = subprocess.run([sys.executable, os.path.join(python_dir, "OFT_ThinCurr_mesh_tool.py"), "modify"] + mesh_tool_args,
                                capture_output=True, text=True)
        if result.returncode != 0:
            print("Error running mesh tool:")
            print(result.stdout)
            print(result.stderr)
            raise RuntimeError("Mesh tool failed with return code {}".format(result.returncode))
        else:
            os.remove("oft_in.xml")
            meshfile = "test.h5"
    # Create flux loop definition file for sensors
    if floops is not None:
        sensors = []
        for (k, floop) in enumerate(floops):
            sensors.append(circular_flux_loop(floop[0],floop[1],'FLOOP_{0}'.format(k),npts=180))
        save_sensors(sensors, 'floops.loc')
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
       return mp_run(run_td,(meshfile,direct_flag,use_aca,floops,curr_waveform,volt_waveform,lin_tol,jumper_start,run_reduced,basepath))
    elif run_type == 2:
        return mp_run(run_eig,(meshfile,direct_flag,use_aca,jumper_start))
    elif run_type == 3:
        return mp_run(run_fr,(meshfile,direct_flag,use_aca,freq,fr_limit,floops,jumper_start))
    elif run_type == 5:
        return mp_run(run_mode,(meshfile,freq))
    elif run_type == 6:
        return mp_run(run_td_for_Mirnov,(meshfile,direct_flag,curr_waveform,lin_tol))


def validate_eigs(eigs, imag_tol=1.E-9):
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
        if val != eigs_run_real[i]:
            print("FAILED: Eigenvalue {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(eigs_run_real[i]))
            retval = False
        if abs(eigs_run_imag[i]) > imag_tol:
            print("FAILED: Imaginary eigenvalue detected!")
            print("  Value =    {0}".format(eigs_run_imag[i]))
            retval = False
    return retval


def validate_fr(fr_real, fr_imag):
    """
    Helper function to validate frequency-response results against test case.
    """
    try:
        fr_run_real = []
        fr_run_imag = []
        with open('thincurr_fr.dat', 'r') as fid:
            for line in fid:
                vals = line.split()
                fr_run_real.append(float(vals[0]))
                fr_run_imag.append(float(vals[1]))
    except BaseException as e:
        print(e)
        return False
    if not len(fr_run_real) == len(fr_real):
        print("FAILED: Number of sensors does not match")
        return False
    retval = True
    for (i, val) in enumerate(fr_real):
        if val != fr_run_real[i]:
            print("FAILED: Real response {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(fr_run_real[i]))
            retval = False
    for (i, val) in enumerate(fr_imag):
        if val != fr_run_imag[i]:
            print("FAILED: Imaginary response {0} incorrect!".format(i+1))
            print("  Expected = {0}".format(val))
            print("  Actual =   {0}".format(fr_run_imag[i]))
            retval = False
    return retval


def validate_td(sigs_final, jumpers_final=None):
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
    if sigs_final[0] != td_sigs_final[0]:
        print("FAILED: Final time incorrect!")
        print("  Expected = {0}".format(sigs_final[0]))
        print("  Actual =   {0}".format(td_sigs_final[0]))
        retval = False
    for (i, val) in enumerate(sigs_final[1:]):
        if val != td_sigs_final[i+1]:
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
        if len(td_sigs_final) != len(jumpers_final):
            print("FAILED: Number of jumpers does not match {0} != {1}".format(len(td_sigs_final), len(jumpers_final)))
            return False
        retval = True
        if jumpers_final[0] != td_sigs_final[0]:
            print("FAILED: Final time incorrect!")
            print("  Expected = {0}".format(jumpers_final[0]))
            print("  Actual =   {0}".format(td_sigs_final[0]))
            retval = False
        for (i, val) in enumerate(jumpers_final[1:]):
            if val is None:
                continue
            if val != td_sigs_final[i+1]:
                print("FAILED: Signal {0} incorrect!".format(i+1))
                print("  Expected = {0}".format(val))
                print("  Actual =   {0}".format(td_sigs_final[i+1]))
                retval = False
    return retval


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
        model_surface.compute_Rmat()
        R_surface = model_surface.Rmat.toarray()

        model_bulk = _build_dummy_model(xml_filename)
        model_bulk.set_eta_values(eta_vol=eta_bulk, thickness=thickness_values)
        model_bulk.compute_Rmat()
        R_bulk = model_bulk.Rmat.toarray()

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
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_plate(direct_flag,convert_xml):
    eigs = approx_list((9.735667E-3, 6.532314E-3, 6.532201E-3, 5.251598E-3), rel=1.E-5)
    assert ThinCurr_setup("tw_test-plate.h5",2,direct_flag,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_plate(direct_flag,convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((8.442110E-4, 7.117118E-4), rel=1.E-3)
    assert ThinCurr_setup("tw_test-plate.h5",1,direct_flag,convert_xml=convert_xml,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_plate(direct_flag,convert_xml):
    fr_real = approx_list((6.807649E-2, 7.207748E-2), rel=1.E-4)
    fr_imag = approx_list((-3.011666E-3, -2.177010E-3), rel=1.E-4)
    assert ThinCurr_setup("tw_test-plate.h5",3,direct_flag,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_plate_volt(direct_flag,convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.580643E-4, 3.854292E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1697.895], rel=1.E-3)
    assert ThinCurr_setup("tw_test-plate.h5",1,direct_flag,convert_xml=convert_xml,
                           vcoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           volt_waveform=((0.0, 1.0), (1.0, 1.0)))
    assert validate_td(sigs_final,jumpers_final)

#============================================================================
# Test runners for cylinder
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_cyl(direct_flag,convert_xml):
    eigs = approx_list((2.657195E-2, 1.248071E-2, 1.247103E-2, 1.200566E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-cyl.h5",2,direct_flag,jumper_start=2,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_cyl(direct_flag,convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.178595E-4, 6.040177E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.451791E3, 5.451791E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",1,direct_flag,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_cyl(direct_flag,convert_xml):
    fr_real = approx_list((6.118337E-2, 4.356188E-3), rel=1.E-4)
    fr_imag = approx_list((-1.911861E-3, -2.283493E-3), rel=1.E-4)
    assert ThinCurr_setup("tw_test-cyl.h5",3,direct_flag,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_cyl_volt(direct_flag,convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.504279E-4, 1.276624E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.1203960E3, 1120.396, 655.853, 655.850], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl.h5",1,direct_flag,convert_xml=convert_xml,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)

#============================================================================
# Test runners for torus
@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_torus(direct_flag,convert_xml):
    eigs = approx_list((4.751344E-2, 2.564491E-2, 2.555695E-2, 2.285850E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-torus.h5",2,direct_flag,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_torus(direct_flag,convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.775732E-4, 3.408042E-5), rel=1.E-3)
    assert ThinCurr_setup("tw_test-torus.h5",1,direct_flag,convert_xml=convert_xml,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)),
                           lin_tol=1.E-10)
    assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_torus(direct_flag,convert_xml):
    fr_real = approx_list((-2.807955E-3, -1.196091E-4), rel=1.E-4)
    fr_imag = approx_list((-1.869732E-3, -1.248642E-4), rel=1.E-4)
    assert ThinCurr_setup("tw_test-torus.h5",3,direct_flag,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_torus_volt(direct_flag,convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((5.653338E-5, 4.035387E-6), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([None, -597.6068, 371.74769, 371.74780], rel=1.E-3)
    assert ThinCurr_setup("tw_test-torus.h5",1,direct_flag,convert_xml=convert_xml,
                           vcoils=((1.5, 0.5), (1.5, -0.5)),
                           floops=((1.4, 0.0), (0.6, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           lin_tol=1.E-11)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_torus_fourier_sensor(direct_flag):
    from OpenFUSIONToolkit.ThinCurr.util import torus_fourier_sensor
    import netCDF4
    R_0 = 1.0
    with netCDF4.Dataset(os.path.join(test_dir,"torus_gpec_control_output_n1_nc.nc"), "r") as file:
        R_gpec = np.asarray(file.variables["R"][:-1])
        Z_gpec = np.asarray(file.variables["z"][:-1])
        delta_phi = np.asarray(file.variables["delta_phi"][:-1])
    interface_h1 = torus_fourier_sensor(R_gpec,Z_gpec,R_0,1)
    interface_h1.place_normal_sensors(nphi=15,filename='floops.loc')
    interface_hminus1 = torus_fourier_sensor(R_gpec,Z_gpec,R_0,-1)
    interface_hminus1.place_normal_sensors(nphi=15,filename='floops.loc')
    t = 80
    assert ThinCurr_setup("tw_test-torus.h5",6,direct_flag,
                           curr_waveform=((0.0, 1.E6), (4.E-3, 0.0), (1.0, 0.0)),
                           lin_tol=1.E-10)
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
    # save_spectrum: default (surfmn) output, injected-mesh path, and vac3d format
    interface_h1.save_spectrum(t,'tmp_spec_a',hamada_dphi=delta_phi)
    interface_h1.save_spectrum(t,'tmp_spec_b',hamada_dphi=delta_phi,data_type='vac3d')
    interface_h1.save_spectrum(t,'tmp_spec_c',hamada_dphi=delta_phi,sensor_mesh=interface_h1.get_B_mesh(t))
    with open('tmp_spec_a.dat','r') as fid:
        content_a = fid.read()
    with open('tmp_spec_b.dat','r') as fid:
        content_b = fid.read()
    with open('tmp_spec_c.dat','r') as fid:
        content_c = fid.read()
    assert len(content_a) > 0
    assert content_c == content_a  # injecting the mesh the file is built from is a no-op
    assert content_b != content_a  # vac3d uses a different (scientific) format
    assert 'E' in content_b.upper()
    for tmp_file in ('tmp_spec_a.dat','tmp_spec_b.dat','tmp_spec_c.dat'):
        os.remove(tmp_file)


@pytest.mark.coverage
def test_torus_fourier_sensor_from_eqdsk():
    from OpenFUSIONToolkit.ThinCurr.util import torus_fourier_sensor
    from OpenFUSIONToolkit.TokaMaker.eqdsk import read_geqdsk
    eqdsk_file = os.path.join(test_dir,'ITER_test.eqdsk')
    eq = read_geqdsk(eqdsk_file)
    theta = np.linspace(0.0,2.0*np.pi,65)
    R = 6.2+2.0*np.cos(theta)
    Z = 2.0*np.sin(theta)
    dphi = 0.3*np.sin(theta)+0.1*np.cos(2.0*theta)
    # major_radius and helicity are derived from the g-file
    interface = torus_fourier_sensor.from_eqdsk(R,Z,eqdsk_file,hamada_dphi=dphi,verbose=False)
    assert abs(interface.major_radius-float(eq.R_mag)) < 1.E-12
    assert interface.helicity == (1 if float(eq.B_center)*float(eq.Ip) > 0.0 else -1)
    # duplicate endpoint is trimmed from the surface and hamada_dphi together
    assert interface.ntheta == 64
    assert interface.hamada_dphi.shape[0] == 64
    # selectable poloidal-angle origin
    interface_c = torus_fourier_sensor.from_eqdsk(R,Z,eqdsk_file,center='rcentr',verbose=False)
    assert abs(interface_c.major_radius-float(eq.R_center)) < 1.E-12
    interface_g = torus_fourier_sensor.from_eqdsk(R,Z,eqdsk_file,center='geometric',verbose=False)
    assert abs(interface_g.major_radius-0.5*(R.max()+R.min())) < 1.E-12
    with pytest.raises(ValueError):
        torus_fourier_sensor.from_eqdsk(R,Z,eqdsk_file,center='unknown',verbose=False)

@pytest.mark.coverage
def test_torus_fourier_sensor_from_gpec_control_output():
    import netCDF4
    from OpenFUSIONToolkit.ThinCurr.util import torus_fourier_sensor
    control_file = os.path.join(test_dir,'torus_gpec_control_output_n1_nc.nc')
    with netCDF4.Dataset(control_file,'r') as ds:
        R_gpec = np.asarray(ds.variables['R'][:])
        Z_gpec = np.asarray(ds.variables['z'][:])
        delta_phi = np.asarray(ds.variables['delta_phi'][:])
        ro = float(ds.getncattr('ro'))
        expected_helicity = 1 if float(ds.getncattr('bt0'))*float(ds.getncattr('crnt')) > 0.0 else -1
    interface = torus_fourier_sensor.from_gpec_control_output(control_file,verbose=False)
    assert abs(interface.major_radius-ro) < 1.E-12
    assert interface.helicity == expected_helicity
    # identical to building positionally from the same file contents
    twin = torus_fourier_sensor(R_gpec,Z_gpec,ro,expected_helicity,hamada_dphi=delta_phi)
    assert np.allclose(interface.radial_positions,twin.radial_positions)
    assert np.allclose(interface.axial_positions,twin.axial_positions)
    assert np.allclose(interface.hamada_dphi,twin.hamada_dphi)
    # explicit overrides win over the file contents
    interface_o = torus_fourier_sensor.from_gpec_control_output(control_file,major_radius=1.5,helicity=-1,verbose=False)
    assert interface_o.major_radius == 1.5
    assert interface_o.helicity == -1
    with pytest.raises(ValueError):
        torus_fourier_sensor.from_gpec_control_output(control_file,center='unknown',verbose=False)

@pytest.mark.coverage
def test_torus_fourier_sensor_hamada_alignment():
    from OpenFUSIONToolkit.ThinCurr.util import torus_fourier_sensor
    R_0 = 3.0
    theta = np.linspace(0.0,2.0*np.pi,65)
    R = R_0+1.0*np.cos(theta)
    Z = 1.0*np.sin(theta)
    dphi = 0.3*np.sin(theta)+0.1*np.cos(2.0*theta)
    base = torus_fourier_sensor(R,Z,R_0,1,hamada_dphi=dphi)
    # the same surface supplied in a rotated point order must produce the same
    # object: hamada_dphi is reordered together with the surface points
    roll = 17
    rolled = torus_fourier_sensor(np.roll(R[:-1],roll),np.roll(Z[:-1],roll),R_0,1,
                                  hamada_dphi=np.roll(dphi[:-1],roll))
    assert np.allclose(base.radial_positions,rolled.radial_positions)
    assert np.allclose(base.axial_positions,rolled.axial_positions)
    assert np.allclose(base.hamada_dphi,rolled.hamada_dphi)
    # methods fall back to the stored hamada_dphi and give identical spectra
    B = np.random.default_rng(1).normal(size=(base.ntheta,8))
    B_n_base, _, _ = base.fft2(B)
    B_n_rolled, _, _ = rolled.fft2(B)
    assert np.allclose(B_n_base,B_n_rolled)
    B_n_expl, _, _ = base.fft2(B,hamada_dphi=base.hamada_dphi)
    assert np.allclose(B_n_base,B_n_expl)

@pytest.mark.coverage
def test_thincurr_model_prep_utils(tmp_path):
    import netCDF4
    from OpenFUSIONToolkit.ThinCurr.util import (drive_to_array, triangular_waveform,
                                                 shift_coils_in_xml, parse_coils_xml,
                                                 find_coil_current_column, add_gpec_coils_to_xml)
    # triangular_waveform: 3-point ramp holding the requested peak at twidth/2
    curr = np.array([14.0, 0.0, 2.0, -4.0])
    tri = triangular_waveform(1.E-2,0.015,curr)
    assert tri.shape == (3,4)
    assert np.allclose(tri[:,0],[0.0,5.E-3,1.E-2])
    assert np.allclose(tri[1,1:],curr[1:]*1.015)
    assert np.allclose(tri[2,1:],curr[1:])
    # drive_to_array: header is 'ncols ntimes', rows are the waveform
    drive_file = tmp_path / 'test.drive'
    drive_file.write_text('3 2\n0.0 1.0 2.0\n1.0 3.0 4.0\n')
    waveform = drive_to_array(str(drive_file))
    assert waveform.shape == (2,3)
    assert np.allclose(waveform[1],[1.0,3.0,4.0])
    # shift_coils_in_xml/parse_coils_xml round-trip
    xml_file = tmp_path / 'coils.xml'
    xml_file.write_text('<oft><thincurr><icoils><coil_set><coil npts="2">\n'
                        '1.0, 2.0, 3.0\n4.0, 5.0, 6.0\n</coil></coil_set>'
                        '</icoils></thincurr></oft>\n')
    shifted_file = tmp_path / 'coils_shifted.xml'
    nshifted = shift_coils_in_xml(str(xml_file),str(shifted_file),dx=0.5,dz=-1.0)
    assert nshifted == (1,2)
    coils = parse_coils_xml(str(shifted_file))
    assert len(coils) == 1
    assert np.allclose(coils[0],[[1.5,2.0,2.0],[4.5,5.0,5.0]])
    bad_file = tmp_path / 'coils_bad.xml'
    bad_file.write_text('<oft><thincurr><icoils><coil_set><coil>\n1.0, 2.0, 3.0\n'
                        '</coil></coil_set></icoils></thincurr></oft>\n')
    with pytest.raises(ValueError):
        shift_coils_in_xml(str(bad_file),str(tmp_path / 'unused.xml'),dx=0.5)
    # find_coil_current_column: exact, suffix-stripped, and numbered fallbacks
    header_to_idx = {'cs1u current [a/turn]': 3, 'divl_1 current [a/turn]': 7}
    assert find_coil_current_column('CS1U_feed_2',header_to_idx) == 3
    assert find_coil_current_column('divl_1',header_to_idx) == 7
    assert find_coil_current_column('unknown',header_to_idx) is None
    # add_gpec_coils_to_xml: split coil groups from a GPEC-style netCDF
    nc_file = tmp_path / 'coils.nc'
    with netCDF4.Dataset(str(nc_file),'w',format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('npts',5)
        ds.createDimension('nsub',1)
        ds.createDimension('ncoil',2)
        pts = np.linspace(0.0,1.0,5)
        for name, vals in (('divl_x',pts), ('divl_y',pts+1.0), ('divl_z',pts-1.0)):
            var = ds.createVariable(name,'f8',('npts','nsub','ncoil'))
            var[:,0,0] = vals
            var[:,0,1] = vals+10.0
        ds.createVariable('divl_nw','f8',())[...] = 2.0
        ds.createVariable('divl_current','f8',('ncoil',))[:] = [1.0,2.0]
    base_xml = tmp_path / 'base.xml'
    base_xml.write_text('<oft><thincurr><eta>1.0</eta></thincurr></oft>\n')
    out_xml = tmp_path / 'with_coils.xml'
    names = add_gpec_coils_to_xml(str(base_xml),str(nc_file),str(out_xml),verbose=False)
    assert names == ['divl_1','divl_2']
    coils = parse_coils_xml(str(out_xml))
    assert len(coils) == 2
    assert coils[0].shape == (5,3)
    with pytest.raises(ValueError):
        add_gpec_coils_to_xml(str(base_xml),str(nc_file),str(out_xml),
                              prefixes={'missing'},verbose=False)

@pytest.mark.coverage
def test_make_drive_and_xml_from_eqdsks(tmp_path):
    pytest.importorskip('shapely')
    import matplotlib
    matplotlib.use('Agg')
    from OpenFUSIONToolkit.ThinCurr.util import make_drive_and_xml_from_eqdsks, drive_to_array, parse_coils_xml
    from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk
    eqdsk_file = os.path.join(test_dir,'ITER_test.eqdsk')
    xml_file = tmp_path / 'filaments.xml'
    drive_file = tmp_path / 'filaments.drive'
    make_drive_and_xml_from_eqdsks(1.0,eqdsk_file,str(xml_file),str(drive_file),[1.E-5,2.E-5])
    coils = parse_coils_xml(str(xml_file))
    waveform = drive_to_array(str(drive_file))
    # filaments outside the LCFS are excluded and columns match the coil_sets
    assert len(coils) > 0
    assert waveform.shape == (2,len(coils)+1)
    assert np.all(waveform[1,1:] != 0.0)
    # first row holds the initial currents at t=0
    assert waveform[0,0] == 0.0
    assert np.allclose(waveform[0,1:],waveform[1,1:])
    # total filament current approximates the equilibrium plasma current
    eqdsk_obj = read_eqdsk(eqdsk_file)
    assert abs(np.sum(waveform[1,1:])-eqdsk_obj['ip'])/abs(eqdsk_obj['ip']) < 0.05

#============================================================================
# Test runners for filament model
@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_eig_passive(direct_flag):
    eigs = approx_list((1.503561E-1, 6.420533E-2, 3.188782E-2, 2.941118E-2), rel=1.E-5)
    assert ThinCurr_setup(None,2,direct_flag,eta=1.E4,
                           vcoils=((0.5, 0.1), (0.5, 0.05),
                                   (0.5, -0.05), (0.5, -0.1)))
    assert validate_eigs(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_passive(direct_flag):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((8.326557E-4, 8.347735E-4), rel=1.E-3)
    assert ThinCurr_setup(None,1,direct_flag,eta=1.E4,
                          icoils=((0.5, 0.1),),
                          vcoils=((0.5, 0.0),),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          curr_waveform=((-1.0, 0.0), (0.0, 0.0), (1.0, 1.0)))
    assert validate_td(sigs_final)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_passive(direct_flag):
    fr_real = approx_list((1.947713E-1, 1.990873E-1), rel=1.E-4)
    fr_imag = approx_list((-2.175942E-4, -1.560726E-4), rel=1.E-4)
    assert ThinCurr_setup(None,3,direct_flag,eta=1.E4,freq=5.E3,fr_limit=0,
                           icoils=((0.5, 0.1),),
                           vcoils=((0.5, 0.0),),
                           floops=((0.5, -0.05), (0.5, -0.1)))
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_td_passive_volt(direct_flag):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((4.379235E-4, 4.389248E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([-641.4736, 1673.2893], rel=1.E-3)
    assert ThinCurr_setup(None,1,direct_flag,eta=1.E4,
                          vcoils=((0.5, 0.0), (0.5, 0.1)),
                          floops=((0.5, -0.05), (0.5, -0.1)),
                          volt_waveform=((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)))
    assert validate_td(sigs_final, jumpers_final)

#============================================================================
# Test for output directory (basepath) support in run_td / plot_td
@pytest.mark.coverage
def test_td_output_dir():
    assert ThinCurr_setup("tw_test-plate.h5", 1, 'F',
                           icoils=((0.5, 0.1),),
                           floops=((0.5, -0.05), (0.5, -0.1)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           basepath='td_output_dir_test')


#============================================================================
# Test runners for large cylinder (w/ ACA+)
@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_aca(convert_xml):
    eigs = approx_list((2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",2,'F',use_aca=True,jumper_start=2,convert_xml=convert_xml)
    assert validate_eigs(eigs)

@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_aca(convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.205529E-4, 6.100309E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.453372E3, 5.453372E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",1,'F',use_aca=True,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2,run_reduced=True)
    assert validate_td(sigs_final,jumpers_final)

@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_aca(convert_xml):
    fr_real = approx_list((5.888736E-2, 4.881440E-3), rel=1.E-3)
    fr_imag = approx_list((-2.017045E-3, -2.313881E-3), rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",3,'F',use_aca=True,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)

@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_volt_aca(convert_xml):
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.512679E-4, 1.291681E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.122550E3, 1122.550, 656.9544, 656.9797], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",1,'F',use_aca=True,convert_xml=convert_xml,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)
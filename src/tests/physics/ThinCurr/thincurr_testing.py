#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Shared helpers for ThinCurr tests: case setup, process runners, and result validation
import os
import sys
import glob
import shutil
import pathlib
import subprocess
import numpy as np
import h5py
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
import oft_testing
from OpenFUSIONToolkit.io import histfile, write_oft_xml
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_Vcoil, ThinCurr_XML
from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy
from OpenFUSIONToolkit.ThinCurr.sensor import save_sensors, circular_flux_loop


mu0 = np.pi*4.E-7
# Calculations run by ThinCurr_setup()
THINCURR_RUN_TYPES = ('td', 'eig', 'fr', 'td_mirnov')
_oft_env_singleton = None

# Circular coils, defined by radius and height in the XML file ('xml') or the mesh file ('mesh'), or converted from
# point lists by OFT_ThinCurr_mesh_tool ('converted')
circular_coil_defs = {
    'xml': {'circular': True},
    'mesh': {'circular': True, 'convert_xml': True},
    'converted': {'convert_xml': True, 'convert_circular': True}
}


def mp_run(target,args,timeout=180):
    return oft_testing.mp_run(target,args,test_dir,timeout=timeout)


def setup_test_model(tw_model, meshfile, jumper_start, nsplit=0):
    '''Load `meshfile`, with the coil/resistivity XML written by ThinCurr_setup() if it exists, or build the dummy
    plate with a high surface resistivity when `meshfile` is None (filament-only cases)'''
    if meshfile is None:
        r_dummy, lc_dummy = build_ThinCurr_dummy([0.0,0.0,10.0],size=0.25,nsplit=nsplit)
        tw_model.setup_model(r=r_dummy,lc=lc_dummy,xml_filename='oft_in.xml',jumper_start=jumper_start)
        tw_model.set_eta_values(eta_values=np.r_[1.E4*mu0])
    elif pathlib.Path("oft_in.xml").is_file():
        tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml',jumper_start=jumper_start)
    else:
        tw_model.setup_model(mesh_file=meshfile)


def write_flux_loops(floops):
    '''Save circular flux loops at the (R, Z) positions `floops` to "floops.loc"'''
    sensors = [circular_flux_loop(R,Z,'FLOOP_{0}'.format(k),npts=180) for k, (R, Z) in enumerate(floops)]
    save_sensors(sensors, 'floops.loc')


def run_td(meshfile,direct_flag,use_aca,floops,curr_waveform,volt_waveform,lin_tol,jumper_start,run_reduced,basepath,hodlr_opts,mp_q):
    '''Time-dependent run driven by coil currents or voltages, with flux-loop sensors and optional model reduction'''
    result = True
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        setup_test_model(tw_model, meshfile, jumper_start)
        tw_model.setup_io(basepath=basepath)
        if floops is not None:
            write_flux_loops(floops)
            _, _, sensor_obj = tw_model.compute_Msensor('floops.loc')
        else:
            sensor_obj = None
        if curr_waveform is not None:
            curr_waveform = np.array(curr_waveform)
            curr_waveform[:,1:] /= mu0
        if volt_waveform is not None:
            volt_waveform = np.array(volt_waveform)
        tw_model.compute_Mcoil()
        if hodlr_opts is not None:
            tw_model.set_hodlr_options(**hodlr_opts)
        tw_model.compute_Lmat(use_hodlr=use_aca,cache_file='Lmat.save')
        tw_model.compute_Rmat()
        rst_before = set(glob.glob('pThinCurr_*.rst'))
        hist_before = os.path.exists('floops.hist')
        tw_model.run_td(2.E-5,200,direct=(direct_flag == 'T'),lin_tol=lin_tol,coil_currs=curr_waveform,coil_volts=volt_waveform,sensor_obj=sensor_obj)
        tw_model.plot_td(200,sensor_obj=sensor_obj)
        if basepath is not None:
            # Output must go to basepath, not the working directory
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
            tw_model.reconstruct_current(tw_reduced.reconstruct_potential(currents['curr'][10]),centering='vertex')
    except BaseException as e:
        print(e)
        result = False
    finally:
        if basepath is not None and os.path.isdir(basepath):
            shutil.rmtree(basepath)
    oftpy_dump_cov()
    mp_q.put(result)


def run_eig(meshfile,direct_flag,use_aca,jumper_start,mp_q):
    '''The four slowest L/R eigenmodes, written to "thincurr_eigs.dat"'''
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        setup_test_model(tw_model, meshfile, jumper_start, nsplit=1)
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


def run_fr(meshfile,direct_flag,use_aca,freq,fr_limit,floops,jumper_start,hodlr_opts,mp_q):
    '''Frequency response to the first coil at `freq`, as flux-loop signals written to "thincurr_fr.dat"'''
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        tw_model = ThinCurr(myOFT)
        setup_test_model(tw_model, meshfile, jumper_start, nsplit=1)
        tw_model.setup_io()
        if floops is not None:
            write_flux_loops(floops)
            Msensor, Msc, _ = tw_model.compute_Msensor('floops.loc')
        Mcoil = tw_model.compute_Mcoil()
        if hodlr_opts is not None:
            tw_model.set_hodlr_options(**hodlr_opts)
        tw_model.compute_Lmat(use_hodlr=use_aca,cache_file='Lmat.save')
        tw_model.compute_Rmat()
        # Coupling matrices omit mu0, so a 1/mu0 A drive gives unit-normalized signals
        driver_current = 1.0/mu0
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


def run_td_for_Mirnov(meshfile,direct_flag,curr_waveform,lin_tol,mp_q):
    '''Time-dependent run with the coils of "oft_in_remc.xml", recording the sensors in "floops.loc" to the history
    file for Fourier analysis'''
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
                    lin_tol=1.E-9,jumper_start=0,run_reduced=False,basepath=None,convert_xml=False,
                    hodlr_opts=None,circular=False,convert_circular=False):
    """
    Common setup and run operations for thin-wall physics module test cases

    `run_type` selects the calculation: 'td' (time-dependent), 'eig' (L/R eigenmodes), 'fr' (frequency
    response), or 'td_mirnov' (time-dependent run recording the sensors in 'floops.loc' to the history file).
    """
    if run_type not in THINCURR_RUN_TYPES:
        raise ValueError("Unknown ThinCurr run type '{0}' (expected one of {1})".format(run_type, ', '.join(THINCURR_RUN_TYPES)))
    nPhi = 180
    phi_fac = 2.0*np.pi/(nPhi-1)
    os.chdir(test_dir)
    xml = ThinCurr_XML()
    xml.set_eta([eta*mu0])
    if icoils is not None:
        test_Icoil = ThinCurr_Icoil(name="Test")
        for icoil in icoils:
            R = icoil[0]; Z = icoil[1]
            if circular:
                test_Icoil.add_circular_subcoil([R, Z], npoints=nPhi)
                continue
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
                if circular:
                    test_Vcoil.add_circular_subcoil([R, Z], npoints=nPhi)
                    xml.add_Vcoil(test_Vcoil)
                    continue
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
            if convert_circular:
                mesh_tool_args += ["--circular_coils"]
        if jumper_start > 0:
            mesh_tool_args += ["--jumper_range", "{0}".format(jumper_start-1)]
        result = subprocess.run([sys.executable, os.path.join(oft_testing.oft_python_dir, "OFT_ThinCurr_mesh_tool.py"), "modify"] + mesh_tool_args,
                                capture_output=True, text=True)
        if result.returncode != 0:
            print("Error running mesh tool:")
            print(result.stdout)
            print(result.stderr)
            raise RuntimeError("Mesh tool failed with return code {}".format(result.returncode))
        else:
            os.remove("oft_in.xml")
            meshfile = "test.h5"
    if floops is not None:
        write_flux_loops(floops)
    # Drive waveforms as (time, value...) rows; currents are given in mu0-normalized units
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
    if run_type == 'td':
        return mp_run(run_td,(meshfile,direct_flag,use_aca,floops,curr_waveform,volt_waveform,lin_tol,jumper_start,run_reduced,basepath,hodlr_opts))
    elif run_type == 'eig':
        return mp_run(run_eig,(meshfile,direct_flag,use_aca,jumper_start))
    elif run_type == 'fr':
        return mp_run(run_fr,(meshfile,direct_flag,use_aca,freq,fr_limit,floops,jumper_start,hodlr_opts))
    else:
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
        if jumpers_final[0] != td_sigs_final[0]:
            print("FAILED: Jumper final time incorrect!")
            print("  Expected = {0}".format(jumpers_final[0]))
            print("  Actual =   {0}".format(td_sigs_final[0]))
            retval = False
        for (i, val) in enumerate(jumpers_final[1:]):
            if val is None:
                continue
            if val != td_sigs_final[i+1]:
                print("FAILED: Jumper {0} incorrect!".format(i+1))
                print("  Expected = {0}".format(val))
                print("  Actual =   {0}".format(td_sigs_final[i+1]))
                retval = False
    return retval


def write_thincurr_xml(xml_filename, eta_values=None, thickness_values=None, eta_vol_values=None):
    from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_XML

    thincurr_xml = ThinCurr_XML()
    if eta_values is not None:
        thincurr_xml.set_eta(eta_values)
    if eta_vol_values is not None:
        thincurr_xml.set_eta_vol(eta_vol_values)
    if thickness_values is not None:
        thincurr_xml.set_thickness(thickness_values)
    write_oft_xml([thincurr_xml], xml_filename)


def get_oft_env():
    '''OFT environment shared by the models built in this process (created on first use)'''
    global _oft_env_singleton
    if _oft_env_singleton is None:
        from OpenFUSIONToolkit import OFT_env
        _oft_env_singleton = OFT_env(nthreads=-1)
    return _oft_env_singleton


def save_dummy_mesh(filename, eta_surf):
    '''Save the dummy mesh with its surface resistivity (no XML file needed) to `filename`'''
    from OpenFUSIONToolkit.ThinCurr._mesh_cli import ThinCurrMesh
    from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy
    r_dummy, lc_dummy = build_ThinCurr_dummy([0.0, 0.0, 10.0], size=0.25)
    mesh = ThinCurrMesh(r_dummy, lc_dummy, np.ones((lc_dummy.shape[0],), dtype=np.int32))
    mesh.set_resistivity(eta_surf=eta_surf)
    mesh.save(filename)


def build_dummy_model(xml_filename):
    from OpenFUSIONToolkit.ThinCurr import ThinCurr
    from OpenFUSIONToolkit.ThinCurr.meshing import build_ThinCurr_dummy

    r_dummy, lc_dummy = build_ThinCurr_dummy([0.0, 0.0, 10.0], size=0.25)
    model = ThinCurr(get_oft_env())
    model.setup_model(r=r_dummy, lc=lc_dummy, xml_filename=xml_filename)
    return model

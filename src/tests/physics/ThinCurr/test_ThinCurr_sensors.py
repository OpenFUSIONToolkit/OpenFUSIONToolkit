#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: Mirnov sensor geometry and couplings, and toroidal Fourier sensor arrays on a driven torus
import os
import pytest
import numpy as np
from thincurr_testing import ThinCurr_setup, mp_run, mu0
from OpenFUSIONToolkit.io import write_oft_xml
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_XML
from OpenFUSIONToolkit.ThinCurr.sensor import save_sensors, Mirnov
test_dir = os.path.abspath(os.path.dirname(__file__))


def test_mirnov_orientation_and_validation():
    '''Mirnov sensor loops are oriented along the requested normal, and zero normals or locations without three
    components are rejected'''
    for norm in ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.3, -0.5, 0.8]):
        pts = Mirnov([1.0, 0.5, -0.2], norm, 'test')._pts
        area = 0.5*sum(np.cross(pts[i], pts[i+1]) for i in range(len(pts)-1))
        n_hat = np.asarray(norm)/np.linalg.norm(norm)
        assert np.allclose(area/np.linalg.norm(area), n_hat)
    with pytest.raises(ValueError):
        Mirnov([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 'test')
    with pytest.raises(ValueError):
        Mirnov([1.0, 0.0], [0.0, 0.0, 1.0], 'test')


def validate_torus_fourier_sensor(interface,sigs_nmodes_1D_PEST,sigs_nmodes_1D_Hamada,sigs_mnmodes_2D_PEST,sigs_mnmodes_2D_Hamada,t,delta_phi,tol=1.E-6):
    '''Compare 1D toroidal and 2D poloidal/toroidal Fourier spectra (PEST and Hamada angles) at step `t` with the
    stored references, saving any mismatching spectrum as `*-new.npy` for inspection'''
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
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_torus_fourier_sensor(direct_flag):
    '''Toroidal Fourier sensor arrays placed normal to a GPEC control surface on the driven torus, for both helicities:
    1D and 2D Fourier spectra of the recorded fields in PEST and Hamada angles match stored references'''
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
    # Coil current from oft_in_remc.xml ramped down from 1 MA to zero over 4 ms
    assert ThinCurr_setup("tw_test-torus.h5",'td_mirnov',direct_flag,
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


def run_mirnov_point(meshfile,mp_q):
    """Couplings of point Mirnov sensors and of small loop Mirnov sensors at the same locations"""
    try:
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        myOFT = OFT_env(nthreads=-1)
        probes = (([1.38, 0.0, 0.0], [1.0, 0.0, 0.0]), ([1.36, 0.05, 0.1], [0.3, 0.2, 0.9]),
                  ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]))
        result = []
        # Point sensors (no dx) and small square loops of side dx at the same locations
        for dx in (None, 1.E-4):
            save_sensors([Mirnov(p, n, 'B{0}'.format(k), dx=dx) for k, (p, n) in enumerate(probes)],
                         'mirnov_point.loc')
            tw_model = ThinCurr(myOFT)
            tw_model.setup_model(mesh_file=meshfile,xml_filename='oft_in.xml')
            Ms, Msc, _ = tw_model.compute_Msensor('mirnov_point.loc')
            result.append((np.array(Ms), np.array(Msc)))
        mp_q.put(result)
    except BaseException as e:
        print(e)
        mp_q.put(None)


@pytest.mark.parametrize("coil_type", ("circular", "polygon"))
def test_mirnov_point(coil_type):
    '''Point Mirnov sensor couplings to the wall and to circular or non-planar polygon coils match those of small
    loop sensors at the same locations'''
    os.chdir(test_dir)
    xml = ThinCurr_XML()
    xml.set_eta([10.0*mu0])
    coil = ThinCurr_Icoil(name="Test")
    if coil_type == "circular":
        coil.add_circular_subcoil([1.6, 0.6])
    else:
        # Non-planar polygon coil: segment fields (point sensors) against segment-pair integrals (small loops)
        t = np.linspace(0.0, 2.0*np.pi, 181)
        coil.add_subcoil(pts=np.c_[1.6*np.cos(t), 1.6*np.sin(t), 0.6 + 0.1*np.sin(3.0*t)])
    xml.add_Icoil(coil)
    write_oft_xml([xml], "oft_in.xml")
    result = mp_run(run_mirnov_point, ("tw_test-torus.h5",))
    assert result is not None
    (Ms_point, Msc_point), (Ms_loop, Msc_loop) = result
    # A loop of side dx approaches the point value as dx^2
    assert np.abs(Ms_loop - Ms_point).max() < 1.E-5*np.abs(Ms_point).max()
    assert np.allclose(Msc_loop, Msc_point, rtol=1.E-6, atol=0.0)

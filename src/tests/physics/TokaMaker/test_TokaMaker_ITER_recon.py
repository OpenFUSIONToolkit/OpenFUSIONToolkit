#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: equilibrium reconstruction of the ITER equilibrium from noisy synthetic diagnostics
import json
import random
import pytest
import numpy as np
from tokamaker_testing import mp_run, validate_dict, ITER_eq_dict, solve_ITER_reference
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.util import mu0


def run_ITER_recon_case(mesh_resolution,fe_order,test_type,mp_q):
    from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
    # The reference equilibrium supplies the synthetic measurements
    mygs, EQ_obj = solve_ITER_reference(mesh_resolution,fe_order)
    if mygs is None:
        mp_q.put(None)
        return
    eq_info = EQ_obj.get_stats(li_normalization='ITER')
    Lmat = mygs.get_coil_Lmat()
    # Every fifth of 100 points spaced evenly along the inner vessel wall
    with open('ITER_geom.json','r') as fid:
        ITER_geom = json.load(fid)
        probe_contour = np.asarray(ITER_geom['inner_vv'][0])
    dl_contour = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(probe_contour,axis=0),axis=1))]
    probe_contour_new = np.zeros((100,2))
    probe_contour_new[:,0] = np.interp(np.linspace(0.0,dl_contour[-1],100),dl_contour,probe_contour[:,0])
    probe_contour_new[:,1] = np.interp(np.linspace(0.0,dl_contour[-1],100),dl_contour,probe_contour[:,1])
    B_locs = []
    for i, pt in enumerate(probe_contour_new):
        if i % 5 == 0:
            B_locs.append(pt)
    B_locs = np.asarray(B_locs)
    # Fixed seed so the 5% Gaussian measurement noise, and the reference values, are reproducible
    random.seed(42)
    myrecon = reconstruction(mygs)
    noise_amp = random.gauss(0.0,1.0)
    Ip_noised = eq_info['Ip']*(1.0+noise_amp*0.05)
    myrecon.set_Ip(Ip_noised, err=0.05*eq_info['Ip'])
    flux_vals = []
    field_eval = mygs.get_field_eval('PSI')
    for i in range(B_locs.shape[0]):
        B_tmp = field_eval.eval(B_locs[i,:])
        noise_amp = random.gauss(0.0,1.0)
        flux_vals.append(B_tmp[0])
        psi_val = B_tmp[0]*2.0*np.pi
        myrecon.add_flux_loop(B_locs[i,:], psi_val*(1.0 + noise_amp*0.05), err=abs(psi_val*0.05))
    field_eval = mygs.get_field_eval('B')
    for i in range(B_locs.shape[0]):
        B_tmp = field_eval.eval(B_locs[i,:])
        noise_amp = random.gauss(0.0,1.0)
        myrecon.add_Mirnov(B_locs[i,:], np.r_[1.0,0.0,0.0], B_tmp[0] + noise_amp*abs(B_tmp[0]*0.05), err=abs(B_tmp[0]*0.05))
        noise_amp = random.gauss(0.0,1.0)
        myrecon.add_Mirnov(B_locs[i,:], np.r_[0.0,0.0,1.0], B_tmp[2] + noise_amp*abs(B_tmp[2]*0.05), err=abs(B_tmp[2]*0.05))
    coil_currents, _ = mygs.get_coil_currents()
    for key in coil_currents:
        noise_amp = random.gauss(0.0,1.0)
        coil_currents[key] *= 1.0+noise_amp*0.05
    # Starting equilibrium from a circular guess, constrained by the noisy flux loops and coil currents instead
    # of the true shape
    mygs.set_isoflux(None)
    mygs.set_saddles(None)
    mygs.set_targets(Ip=Ip_noised,Ip_ratio=2.0)
    if test_type == 'recon':
        mygs.settings.ffp_target_weight=1.0/((0.05*abs(Ip_noised))*mu0)
        mygs.settings.pp_target_weight=1.0
        mygs.update_settings()
    mygs.set_flux(B_locs,np.array(flux_vals))
    regularization_terms = []
    for name in coil_currents:
        regularization_terms.append(mygs.coil_reg_term({name: 1.0},target=coil_currents[name],weight=1.E-1))
    regularization_terms.append(mygs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E2))
    mygs.set_coil_reg(reg_terms=regularization_terms)
    coil_err = {name: 0.05*abs(current) for name, current in coil_currents.items()}
    myrecon.set_coil_currents(coil_currents,coil_err)
    R0 = 6.3
    Z0 = 0.5
    a = 1.0
    kappa = 1.0
    delta = 0.0
    mygs.init_psi(R0, Z0, a, kappa, delta)
    if test_type == 'recon':
        mygs.settings.maxits=300
    elif test_type == 'recon_legacy':
        mygs.settings.maxits=100
    mygs.update_settings()
    mygs.solve()
    if test_type == 'recon':
        mygs.settings.pm = False
        mygs.update_settings()
        mygs.set_isoflux(None)
        mygs.set_saddles(None)
        mygs.set_targets(Ip=Ip_noised)
        mygs.settings.ffp_target_weight=1.0/((0.05*abs(Ip_noised))*mu0)
        myrecon.setup_constraints()
        # Fit only the pressure scale (with coil currents) using the Python least-squares driver
        myrecon.settings.fitF = False
        myrecon.settings.fitP = False
        myrecon.settings.fit_Pscale = True
        myrecon.settings.fit_FFPscale = False
        myrecon.settings.fitR0 = False
        myrecon.settings.fitZ0 = False
        myrecon.settings.fixedCentering = False
        myrecon.settings.dx = 1.E-2
        x0_nl, _ = myrecon.setup_get_opt()
        from scipy.optimize import least_squares
        least_squares(myrecon.opt_error, x0_nl, jac=myrecon.opt_error_jacobian, method='lm', ftol=1.E-3, args=(myrecon,False))
        mygs.settings.pm = True
        mygs.update_settings()
    elif test_type == 'recon_legacy':
        # Fortran reconstruction driver, fitting the magnetic axis radius and coil currents
        mygs.set_isoflux(None)
        mygs.set_flux(None,None)
        mygs.set_saddles(None)
        mygs.set_targets(R0=mygs.o_point[0],V0=mygs.o_point[1])
        myrecon.settings.fit_Pscale = False
        myrecon.settings.fitR0 = True
        myrecon.settings.fitCoils = True
        myrecon.settings.pm = False
        _ = myrecon.reconstruct()
    eq_info = mygs.get_stats(li_normalization='ITER')
    eq_info['LCS1'] = Lmat[mygs.coil_sets['CS1U']['id'],mygs.coil_sets['CS1U']['id']]
    eq_info['MCS1_plasma'] = Lmat[mygs.coil_sets['CS1U']['id'],-1]
    eq_info['Lplasma'] = Lmat[-1,-1]
    mygs.save_eqdsk('tokamaker.eqdsk',nr=64,nz=64,lcfs_pad=0.001)
    mygs.save_ifile('tokamaker.ifile',npsi=64,ntheta=64,lcfs_pad=0.001)
    mp_q.put([eq_info])
    oftpy_dump_cov()


def ITER_recon_dict_base():
    '''Global quantities of the ITER reference equilibrium (reconstruction cases do not compute flux-surface averages)'''
    return {key: val for key, val in ITER_eq_dict.items() if not key.startswith('fsa_')}


@pytest.mark.coverage
def test_ITER_recon():
    '''Equilibrium reconstruction of the ITER reference equilibrium from noisy synthetic plasma current, flux loop,
    Mirnov, and coil current measurements with the Python least-squares driver; the reconstructed equilibrium is
    compared to reference values, which stay close to the true equilibrium'''
    ITER_recon_dict = ITER_recon_dict_base()
    ITER_recon_dict['q_0'] = 0.8417344
    ITER_recon_dict['P_ax'] = 6.575446E5
    ITER_recon_dict['W_MHD'] = 2.591990E8
    ITER_recon_dict['beta_pol'] = 46.02736
    ITER_recon_dict['dflux'] = 1.440565
    ITER_recon_dict['l_i'] = 0.8872565
    ITER_recon_dict['beta_tor'] = 1.906180
    ITER_recon_dict['beta_n'] = 1.284312
    results = mp_run(run_ITER_recon_case,(1.0,2,'recon'))
    assert validate_dict(results,ITER_recon_dict)


@pytest.mark.coverage
def test_ITER_recon_legacy():
    '''Reconstruction of the same synthetic ITER measurements with the Fortran reconstruction driver'''
    ITER_recon_dict = ITER_recon_dict_base()
    ITER_recon_dict['vol'] = 8.119402E2
    ITER_recon_dict['q_95'] = 2.718779
    ITER_recon_dict['P_ax'] = 6.656184E5
    ITER_recon_dict['W_MHD'] = 2.607023E8
    ITER_recon_dict['beta_pol'] = 45.63260
    ITER_recon_dict['dflux'] = 1.462672
    ITER_recon_dict['tflux'] = 1.201402E2
    ITER_recon_dict['l_i'] = 0.8906732
    ITER_recon_dict['beta_tor'] = 1.934024
    ITER_recon_dict['beta_n'] = 1.294385
    results = mp_run(run_ITER_recon_case,(1.0,2,'recon_legacy'))
    assert validate_dict(results,ITER_recon_dict)

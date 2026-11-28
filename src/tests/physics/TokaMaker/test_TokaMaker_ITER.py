#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: ITER free-boundary equilibrium, wall modes, linear stability, and equilibrium I/O
import os
import pytest
import numpy as np
from tokamaker_testing import (mp_run, validate_dict, validate_eqdsk, validate_ifile, ifile_errors, ITER_eq_dict, test_dir,
                               create_ITER_mesh, setup_ITER_gs, set_ITER_constraints, solve_ITER_eq, ifile_q_geom, xpoint_theta)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov

# Relative error of q recomputed from X-point-packed i-file geometry (ntheta=129) against the file q, for
# surfaces with normalized flux <= 0.9 and the outermost surface (lcfs_pad=0.001); about 3x the values found
# with orders 2 and 3 (6E-5 and 1.6E-4 at order 2), where uniform angles give 8E-3 at the outermost surface
QGEOM_CORE_TOL = 2.E-4
QGEOM_EDGE_TOL = 5.E-4


def run_ITER_case(mesh_resolution,fe_orders,test_type,helicity,mp_q):
    try:
        create_ITER_mesh(mesh_resolution)
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    mygs = None
    myOFT = OFT_env(nthreads=-1)
    # Several orders build several TokaMaker objects in one OFT environment (see test_ITER_concurrent)
    for fe_order in fe_orders:
        mygs_last = mygs
        # Negative helicity reverses the toroidal field
        mygs = setup_ITER_gs(myOFT,fe_order,F0=helicity*5.3*6.2)
        if test_type.startswith('eig'):
            if test_type == 'eig_dep':
                # Deprecated interface returns 1/tau as (real, imaginary) pairs
                eig_vals, _ = mygs.eig_wall(10)
                mp_q.put([{'Tau_w': 1.0/eig_vals[:5,0]}])
            else:
                eig_vals, _ = mygs.compute_wall_modes(10)
                mp_q.put([{'Tau_w': eig_vals[:5]}])
            oftpy_dump_cov()
            return
        set_ITER_constraints(mygs)
        try:
            EQ_obj, nl_its = solve_ITER_eq(mygs)
        except ValueError:
            mp_q.put(None)
            return
        if test_type.startswith('stab'):
            if test_type == 'stab_dep':
                eig_vals, eig_modes = mygs.eig_td(-1.E2,10,False)
                growth_rates = -eig_vals[:,0]
            else:
                growth_rates, eig_modes = mygs.compute_linear_stability(1.E2,10,False)
            # Seed a brief nonlinear evolution with 1% of the unstable (vertical) mode, with its sign fixed by the
            # largest mode amplitude above the magnetic axis
            _, Z0 = mygs.o_point
            psi0 = mygs.get_psi(False)
            upper = mygs.r[:,1] > Z0
            eig_sign = eig_modes[0,upper][abs(eig_modes[0,upper]).argmax()]
            psi_ic = psi0-0.01*eig_modes[0,:]*(mygs.psi_bounds[1]-mygs.psi_bounds[0])/eig_sign
            mygs.set_psi(psi_ic)
            mygs.set_saddles(None)
            mygs.set_isoflux(None)
            # Ten steps per e-folding time of the instability, run for one e-folding time
            dt = 0.1/abs(growth_rates[0])
            mygs.setup_td(dt,1.E-13,1.E-11)
            sim_time = 0.0
            for i in range(10):
                sim_time, _, _, _, _ = mygs.step_td(sim_time,dt)
            psi1 = mygs.get_psi(False)
            mp_q.put([{'gamma': growth_rates[:5], 'nl_change': np.linalg.norm(psi1-psi0)}])
            oftpy_dump_cov()
            return
        if test_type == 'ifile_theta':
            mp_q.put([ITER_ifile_theta(mygs)])
            oftpy_dump_cov()
            return
        mygs.save_eqdsk('test.eqdsk',lcfs_pressure=6.E4)
        if test_type == 'io':
            # A reloaded equilibrium must already be converged (one nonlinear iteration)
            EQ_obj.save_TokaMaker('test_eq.h5')
            mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
            mygs.replace_eq(source_file='test_eq.h5')
            EQ_obj = mygs.copy_eq()
            _, nl_its = mygs.solve(return_its=True)
        eq_info = EQ_obj.get_stats(li_normalization='ITER')
        Lmat = mygs.get_coil_Lmat()
        eq_info['LCS1'] = Lmat[mygs.coil_sets['CS1U']['id'],mygs.coil_sets['CS1U']['id']]
        eq_info['MCS1_plasma'] = Lmat[mygs.coil_sets['CS1U']['id'],-1]
        eq_info['Lplasma'] = Lmat[-1,-1]
        eq_info['nl_its'] = nl_its
        fsa = EQ_obj.get_fsa(psi=np.r_[0.25, 0.5, 0.9])
        for key in ('q', 'F', 'dV/dPsi', '<|grad psi|>', '<|grad psi|^2>', '<Bp^2>',
                    '<1/B^2>', 'R_min', 'R_max', 'Z_min', 'Z_max', 'R_at_Zmin', 'R_at_Zmax'):
            eq_info['fsa_' + key] = fsa[key].tolist()
    # Deleting the previous object exercises TokaMaker cleanup while another instance is live
    if mygs_last is not None:
        del mygs_last
    mygs.save_eqdsk('tokamaker.eqdsk',nr=64,nz=64,lcfs_pad=0.001)
    # Uniform angles, as in the stored reference i-file (the default packs this diverted equilibrium toward the X-point)
    mygs.save_ifile('tokamaker.ifile',npsi=64,ntheta=64,lcfs_pad=0.001,theta_sampling='uniform')
    mygs.save_ifile('tokamaker_ext.ifile',npsi=64,ntheta=64,lcfs_pad=0.001,profile_derivs=True,psi_bounds=True,theta_sampling='uniform')
    eq_info.update(ifile_errors(mygs,'tokamaker_ext.ifile'))
    mygs.save_ifile('tokamaker_axis.ifile',npsi=64,ntheta=64,lcfs_pad=0.001,psi_bounds=True,profile_derivs=True,axis_pack=1.E-4)
    eq_info['ifile_psi_err_axis'] = ifile_errors(mygs,'tokamaker_axis.ifile')['ifile_psi_err']
    mp_q.put([eq_info])
    oftpy_dump_cov()


def ITER_ifile_theta(mygs):
    """i-files with uniform and X-point-packed poloidal angles and explicit radial nodes, and their checks"""
    from OpenFUSIONToolkit.TokaMaker.util import read_ifile
    info = {}
    def same_bytes(file1,file2):
        with open(file1,'rb') as fid1, open(file2,'rb') as fid2:
            return fid1.read() == fid2.read()
    # The default ('auto') packs this diverted equilibrium: same bytes as explicit X-point packing
    mygs.save_ifile('ifile_default.ifile',npsi=64,ntheta=129,lcfs_pad=0.001,psi_bounds=True,profile_derivs=True)
    mygs.save_ifile('ifile_xpt.ifile',npsi=64,ntheta=129,lcfs_pad=0.001,psi_bounds=True,profile_derivs=True,theta_sampling='xpoint')
    info['default_is_xpt'] = same_bytes('ifile_default.ifile','ifile_xpt.ifile')
    # Uniform sampling ignores the packing arguments
    mygs.save_ifile('ifile_uniform.ifile',npsi=64,ntheta=129,lcfs_pad=0.001,psi_bounds=True,profile_derivs=True,theta_sampling='uniform')
    mygs.save_ifile('ifile_uniform2.ifile',npsi=64,ntheta=129,lcfs_pad=0.001,psi_bounds=True,profile_derivs=True,
                    theta_sampling='uniform',xpoint_amp=5.0,xpoint_width=0.3)
    info['uniform_same_bytes'] = same_bytes('ifile_uniform.ifile','ifile_uniform2.ifile')
    # X-point packed angles: points on their flux surfaces at the reference packed angles
    info['xpt_psi_err'] = ifile_errors(mygs,'ifile_xpt.ifile')['ifile_psi_err']
    ifile_u = read_ifile('ifile_uniform.ifile')
    ifile_x = read_ifile('ifile_xpt.ifile')
    R0, Z0 = ifile_x['R'][0,0], ifile_x['Z'][0,0]
    theta = np.mod(np.arctan2(ifile_x['Z'][1:,:-1]-Z0,ifile_x['R'][1:,:-1]-R0),2.0*np.pi)
    theta[:,0] = np.arctan2(ifile_x['Z'][1:,0]-Z0,ifile_x['R'][1:,0]-R0) # first point at 0 (may round below)
    x_points, diverted = mygs.get_xpoints()
    thx = [np.mod(np.arctan2(x_points[-1,1]-Z0,x_points[-1,0]-R0),2.0*np.pi)]
    theta_ref = xpoint_theta(128,thx,7.0,0.10)
    info['diverted'] = bool(diverted)
    info['xpt_theta_err'] = np.max(np.abs(theta-theta_ref[None,:]))
    info['xpt_theta_monotone'] = bool(np.all(np.diff(theta,axis=1) > 0.0))
    info['xpt_periodic'] = bool(np.array_equal(ifile_x['R'][:,0],ifile_x['R'][:,-1]) and np.array_equal(ifile_x['Z'][:,0],ifile_x['Z'][:,-1]))
    info['xpt_dtheta_ratio'] = np.max(np.diff(theta_ref))/np.min(np.diff(theta_ref))
    # Flux grid, profiles and q do not depend on the poloidal sampling
    info['xpt_same_profiles'] = bool(all(np.array_equal(ifile_u[key],ifile_x[key]) for key in ('psi','f','p','ffp','pp','psi_bounds')))
    info['xpt_q_diff'] = np.max(np.abs(ifile_x['q']-ifile_u['q'])/np.abs(ifile_u['q']))
    # q recomputed from the file geometry, uniform in the file label, against the file q
    qrel_u = np.abs(ifile_q_geom(mygs,ifile_u)/ifile_u['q'][1:]-1.0)
    qrel_x = np.abs(ifile_q_geom(mygs,ifile_x)/ifile_x['q'][1:]-1.0)
    xn = (ifile_u['psi'][1:]-ifile_u['psi_bounds'][0])/(ifile_u['psi_bounds'][1]-ifile_u['psi_bounds'][0])
    info['qgeom_err_uniform_core'] = np.max(qrel_u[xn <= 0.9])
    info['qgeom_err_xpt_core'] = np.max(qrel_x[xn <= 0.9])
    info['qgeom_err_uniform_edge'] = qrel_u[-1]
    info['qgeom_err_xpt_edge'] = qrel_x[-1]
    # Explicit radial nodes are honored (and override npsi, lcfs_pad, pack_lcfs, axis_pack)
    nodes = np.r_[0.0, np.logspace(-4,-1,7), np.linspace(0.2,0.9,8), 1.0-np.logspace(-1.5,-3,6)]
    mygs.save_ifile('ifile_nodes.ifile',npsi=17,lcfs_pad=0.05,pack_lcfs=False,axis_pack=1.E-3,ntheta=65,psi_bounds=True,profile_derivs=True,
                    psi_nodes=nodes,theta_sampling='xpoint')
    ifile_n = read_ifile('ifile_nodes.ifile')
    xn = (ifile_n['psi']-ifile_n['psi_bounds'][0])/(ifile_n['psi_bounds'][1]-ifile_n['psi_bounds'][0])
    info['nodes_npsi'] = int(ifile_n['npsi'])
    info['nodes_err'] = np.max(np.abs(xn-nodes))
    info['nodes_psi_err'] = ifile_errors(mygs,'ifile_nodes.ifile')['ifile_psi_err']
    # Invalid input is rejected
    info['bad_input_rejected'] = []
    for kwargs in ({'theta_sampling': 'polar'}, {'psi_nodes': [0.1, 0.5, 0.9]}, {'psi_nodes': [0.0, 0.5, 0.5, 0.9]},
                   {'psi_nodes': [0.0, 0.5, 1.0]}, {'theta_sampling': 'xpoint', 'xpoint_width': 0.0},
                   {'theta_sampling': 'xpoint', 'xpoint_amp': -1.0}):
        try:
            mygs.save_ifile('ifile_bad.ifile',npsi=16,ntheta=33,**kwargs)
            info['bad_input_rejected'].append(False)
        except (ValueError, Exception):
            info['bad_input_rejected'].append(True)
    return info


def ITER_eq_dict_helicity(helicity):
    eq_dict = ITER_eq_dict.copy()
    eq_dict['tflux'] *= helicity
    eq_dict['dflux'] *= helicity
    eq_dict['q_0'] *= helicity
    eq_dict['q_95'] *= helicity
    # `q` and `F` follow the sign of F0; the other `fsa_*` entries do not. The copy is shallow, so
    # rebind these lists rather than scaling them in place.
    eq_dict['fsa_q'] = [val*helicity for val in eq_dict['fsa_q']]
    eq_dict['fsa_F'] = [val*helicity for val in eq_dict['fsa_F']]
    return eq_dict


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))
def test_ITER_eig(order):
    '''L/R decay times of the slowest ITER vacuum vessel eddy-current modes from the current and deprecated
    wall-mode interfaces, against reference values'''
    exp_dict = {
        'Tau_w': [6.619977E-01, 3.479492E-01, 2.554444E-01, 1.910381E-01, 1.782464E-01]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),'eig',1.0))
    assert validate_dict(results,exp_dict)
    results = mp_run(run_ITER_case,(1.0,(order,),'eig_dep',1.0))
    assert validate_dict(results,exp_dict)


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))
def test_ITER_stability(order):
    '''Linearized axisymmetric stability of the ITER equilibrium with the vessel (current and deprecated
    interfaces), then a short time-dependent run seeded with the unstable vertical mode; growth rates and the
    nonlinear flux change are compared to reference values'''
    exp_dict = {
        'gamma': [12.3620, -1.83981, -3.41613, -5.12470, -6.53393],
        'nl_change': [18.421737542246518, 27.49498837155886][order-2]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),'stab',1.0),timeout=90)
    assert validate_dict(results,exp_dict)
    results = mp_run(run_ITER_case,(1.0,(order,),'stab_dep',1.0),timeout=90)
    assert validate_dict(results,exp_dict)


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))
@pytest.mark.parametrize("helicity", (1.0,-1.0))
def test_ITER_eq(order,helicity):
    '''Free-boundary ITER H-mode equilibrium with shape, X-point, and coil regularization targets, for both
    toroidal field directions: global and flux-surface quantities, coil inductances, and the written gEQDSK and
    i-file are compared to reference values. The optional i-file records (FF', P', psi bounds) and axis packing
    are checked against the equilibrium itself'''
    from OpenFUSIONToolkit.TokaMaker.util import read_ifile
    eq_dict = ITER_eq_dict_helicity(helicity)
    eq_dict['l_i_std'] = 1.1532768229095525 # Both l_i forms are always returned
    eq_dict['l_i_iter'] = eq_dict['l_i']
    results = mp_run(run_ITER_case,(1.0,(order,),'',helicity))
    assert validate_dict(results,eq_dict)
    assert validate_eqdsk('tokamaker.eqdsk','ITER_test.eqdsk',helicity)
    assert validate_ifile('tokamaker.ifile','ITER_test.ifile',helicity)
    # Optional records are written only on request and leave the others unchanged
    ifile = read_ifile(os.path.join(test_dir,'tokamaker.ifile'))
    ifile_ext = read_ifile(os.path.join(test_dir,'tokamaker_ext.ifile'))
    assert ('ffp' not in ifile) and ('psi_bounds' not in ifile)
    assert all(np.array_equal(ifile[key],ifile_ext[key]) for key in ifile)
    # Records after the known ones are ignored (here a repeat of the last, 2-value record)
    with open(os.path.join(test_dir,'tokamaker_ext.ifile'),'rb') as fid:
        content = fid.read()
    with open(os.path.join(test_dir,'tokamaker_trail.ifile'),'wb') as fid:
        fid.write(content + content[-24:])
    ifile_trail = read_ifile(os.path.join(test_dir,'tokamaker_trail.ifile'))
    assert all(np.array_equal(ifile_ext[key],ifile_trail[key]) for key in ifile_ext)
    # i-file points lie on their flux surfaces; FF', P' and psi bounds match the equilibrium
    assert results[0]['ifile_psi_err'] < 1.E-10
    assert results[0]['ifile_ffp_err'] < 1.E-10
    assert results[0]['ifile_pp_err'] < 1.E-10
    assert results[0]['ifile_bounds_err'] < 1.E-12
    # Axis zone: first surface at axis_pack, spacing growing smoothly into the LCFS-packed points
    ifile_axis = read_ifile(os.path.join(test_dir,'tokamaker_axis.ifile'))
    x = (ifile_axis['psi']-ifile_axis['psi_bounds'][0])/(ifile_axis['psi_bounds'][1]-ifile_axis['psi_bounds'][0])
    dx = np.diff(x)
    assert abs(x[1]-1.E-4) < 1.E-10
    assert np.all(dx > 0) and np.max(dx[1:]/dx[:-1]) < 1.7
    assert results[0]['ifile_psi_err_axis'] < 1.E-10


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))
def test_ITER_ifile_theta(order):
    '''i-file poloidal sampling and explicit radial nodes on the ITER equilibrium: the default ('auto') packs this
    diverted equilibrium (same bytes as 'xpoint'); 'uniform' ignores the packing arguments; X-point-packed points lie on their flux surfaces at the reference packed angles (monotone,
    periodic end point repeated), with the same flux grid, profiles and q as uniform sampling; q recomputed from the
    file geometry (uniform in the file label) agrees with the file q, better than uniform angles at the outermost
    surface; explicit radial nodes are honored; invalid input is rejected'''
    results = mp_run(run_ITER_case,(1.0,(order,),'ifile_theta',1.0))
    assert results is not None
    info = results[0]
    assert info['default_is_xpt']
    assert info['uniform_same_bytes']
    assert info['diverted']
    assert info['xpt_psi_err'] < 1.E-10
    assert info['xpt_theta_err'] < 1.E-9
    assert info['xpt_theta_monotone'] and info['xpt_periodic']
    assert info['xpt_dtheta_ratio'] > 7.0
    assert info['xpt_same_profiles']
    assert info['xpt_q_diff'] < 1.E-12
    assert info['qgeom_err_xpt_core'] < QGEOM_CORE_TOL
    assert info['qgeom_err_xpt_edge'] < QGEOM_EDGE_TOL
    assert info['qgeom_err_xpt_edge'] < 0.1*info['qgeom_err_uniform_edge']
    assert info['nodes_npsi'] == 22
    assert info['nodes_err'] < 1.E-12
    assert info['nodes_psi_err'] < 1.E-10
    assert all(info['bad_input_rejected'])


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))
@pytest.mark.parametrize("helicity", (1.0,-1.0))
def test_ITER_eq_io(order,helicity):
    '''ITER equilibrium saved to a TokaMaker file and loaded back with replace_eq: the reloaded state matches the
    reference and is already converged (one nonlinear iteration)'''
    eq_dict = ITER_eq_dict_helicity(helicity)
    eq_dict['nl_its'] = 1
    results = mp_run(run_ITER_case,(1.0,(order,),'io',helicity))
    assert validate_dict(results,eq_dict)


def test_ITER_concurrent():
    '''Two TokaMaker objects of different element order created in sequence in one OFT environment, with the
    first deleted while the second is live; the final equilibrium must match the reference'''
    results = mp_run(run_ITER_case,(1.0,(2,3),'',1.0))
    assert validate_dict(results,ITER_eq_dict)

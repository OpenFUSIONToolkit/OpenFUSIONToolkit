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
                               create_ITER_mesh, setup_ITER_gs, set_ITER_constraints, solve_ITER_eq)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov


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
    mygs.save_ifile('tokamaker.ifile',npsi=64,ntheta=64,lcfs_pad=0.001)
    mygs.save_ifile('tokamaker_ext.ifile',npsi=64,ntheta=64,lcfs_pad=0.001,profile_derivs=True,psi_bounds=True)
    eq_info.update(ifile_errors(mygs,'tokamaker_ext.ifile'))
    mygs.save_ifile('tokamaker_axis.ifile',npsi=64,ntheta=64,lcfs_pad=0.001,psi_bounds=True,profile_derivs=True,axis_pack=1.E-4)
    eq_info['ifile_psi_err_axis'] = ifile_errors(mygs,'tokamaker_axis.ifile')['ifile_psi_err']
    mp_q.put([eq_info])
    oftpy_dump_cov()


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
    results = mp_run(run_ITER_case,(1.0,(order,),'stab',1.0))
    assert validate_dict(results,exp_dict)
    results = mp_run(run_ITER_case,(1.0,(order,),'stab_dep',1.0))
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

#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: ITER free-boundary equilibrium, wall modes, linear stability, and equilibrium I/O
import pytest
import numpy as np
from tokamaker_testing import (mp_run, validate_dict, validate_eqdsk, validate_ifile, ITER_eq_dict,
                               create_ITER_mesh, setup_ITER_gs, set_ITER_constraints, solve_ITER_eq)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov


#============================================================================
def run_ITER_case(mesh_resolution,fe_orders,test_type,helicity,mp_q):
    try:
        create_ITER_mesh(mesh_resolution)
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    mygs = None
    myOFT = OFT_env(nthreads=-1)
    for fe_order in fe_orders:
        mygs_last = mygs
        mygs = setup_ITER_gs(myOFT,fe_order,F0=helicity*5.3*6.2)
        if test_type.startswith('eig'):
            if test_type == 'eig_dep':
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
            # Brief nonlinear evolution seeded with 1% of the least stable mode
            R0 = 6.3
            psi0 = mygs.get_psi(False)
            eig_sign = eig_modes[0,(mygs.r[:,1]-R0)>0.0][abs(eig_modes[0,(mygs.r[:,1]-R0)>0.0]).argmax()]
            psi_ic = psi0-0.01*eig_modes[0,:]*(mygs.psi_bounds[1]-mygs.psi_bounds[0])/eig_sign
            mygs.set_psi(psi_ic)
            mygs.set_saddles(None)
            mygs.set_isoflux(None)
            dt = 0.1/abs(growth_rates[0])
            mygs.setup_td(dt,1.E-13,1.E-11)
            sim_time = 0.0
            for i in range(5):
                sim_time, _, _, _, _ = mygs.step_td(sim_time,dt)
            psi1 = mygs.get_psi(False)
            mp_q.put([{'gamma': growth_rates[:5], 'nl_change': np.linalg.norm(psi1-psi0)}])
            oftpy_dump_cov()
            return
        mygs.save_eqdsk('test.eqdsk',lcfs_pressure=6.E4)
        if test_type == 'io':
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
        # Flux surface averages and per-surface shape parameters (see `get_fsa`)
        fsa = EQ_obj.get_fsa(psi=np.r_[0.25, 0.5, 0.9])
        for key in ('q', 'F', 'dV/dPsi', '<|grad psi|>', '<|grad psi|^2>', '<Bp^2>',
                    '<1/B^2>', 'R_min', 'R_max', 'Z_min', 'Z_max', 'R_at_Zmin', 'R_at_Zmax'):
            eq_info['fsa_' + key] = fsa[key].tolist()
    # Deleting an instance must not disturb the one still in use
    if mygs_last is not None:
        del mygs_last
    mygs.save_eqdsk('tokamaker.eqdsk',nr=64,nz=64,lcfs_pad=0.001)
    mygs.save_ifile('tokamaker.ifile',npsi=64,ntheta=64,lcfs_pad=0.001)
    mp_q.put([eq_info])
    oftpy_dump_cov()


def ITER_eq_dict_helicity(helicity):
    eq_dict = ITER_eq_dict.copy()
    eq_dict['tflux'] *= helicity
    eq_dict['dflux'] *= helicity
    eq_dict['q_0'] *= helicity
    eq_dict['q_95'] *= helicity
    # q and F follow the sign of F0; rebind the lists since the dict copy is shallow
    eq_dict['fsa_q'] = [val*helicity for val in eq_dict['fsa_q']]
    eq_dict['fsa_F'] = [val*helicity for val in eq_dict['fsa_F']]
    return eq_dict


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_eig(order):
    '''Wall eddy-current (L/R) eigenmodes of the ITER vacuum vessel: the five slowest decay times of the
    conducting-region model, through both the current and the deprecated interface.'''
    exp_dict = {
        'Tau_w': [6.619977E-01, 3.479492E-01, 2.554444E-01, 1.910381E-01, 1.782464E-01]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),'eig',1.0))
    assert validate_dict(results,exp_dict)
    results = mp_run(run_ITER_case,(1.0,(order,),'eig_dep',1.0))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
def test_ITER_stability(order):
    '''Linear stability of the free-boundary ITER equilibrium: the growth-rate spectrum of the linearized
    time-dependent system (one vertically unstable mode and stable modes) and the change over a short nonlinear
    evolution seeded with the unstable mode, through both the current and the deprecated interface.'''
    exp_dict = {
        'gamma': [12.3620, -1.83981, -3.41613, -5.12470, -6.53393],
        'nl_change': [225.4421413167051, 338.0113029638385][order-2]
    }
    results = mp_run(run_ITER_case,(1.0,(order,),'stab',1.0))
    assert validate_dict(results,exp_dict)
    results = mp_run(run_ITER_case,(1.0,(order,),'stab_dep',1.0))
    assert validate_dict(results,exp_dict)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
@pytest.mark.parametrize("helicity", (1.0,-1.0))
def test_ITER_eq(order,helicity):
    '''Free-boundary ITER equilibrium with coil regularization, isoflux and X-point constraints: global parameters,
    coil inductances, flux-surface averages and shape, and the gEQDSK and i-file output, for both signs of the
    toroidal field.'''
    eq_dict = ITER_eq_dict_helicity(helicity)
    results = mp_run(run_ITER_case,(1.0,(order,),'',helicity))
    assert validate_dict(results,eq_dict)
    assert validate_eqdsk('tokamaker.eqdsk','ITER_test.eqdsk',helicity)
    assert validate_ifile('tokamaker.ifile','ITER_test.ifile',helicity)

@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,3))#,4))
@pytest.mark.parametrize("helicity", (1.0,-1.0))
def test_ITER_eq_io(order,helicity):
    '''Save/load round trip of the ITER equilibrium through replace_eq(): the reloaded state must reproduce the
    solution and re-converge in a single iteration.'''
    eq_dict = ITER_eq_dict_helicity(helicity)
    eq_dict['nl_its'] = 1
    results = mp_run(run_ITER_case,(1.0,(order,),'io',helicity))
    assert validate_dict(results,eq_dict)

def test_ITER_concurrent():
    '''Two TokaMaker instances of different element order on the same mesh in one environment, deleting the first; the
    second must still reproduce the ITER equilibrium.'''
    results = mp_run(run_ITER_case,(1.0,(2,3),'',1.0))
    assert validate_dict(results,ITER_eq_dict)

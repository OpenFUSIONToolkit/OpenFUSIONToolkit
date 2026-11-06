#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: ITER profiles on normalized toroidal flux (coord='phi_n')
import pytest
import numpy as np
from tokamaker_testing import mp_run, create_ITER_mesh, setup_ITER_gs, set_ITER_constraints, ITER_Hmode_profiles
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov
from OpenFUSIONToolkit.util import eC
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun


# Each case solves ITER with profiles on psi_N, maps them onto Phi_N with the converged q (dPhi/dpsi = q),
# re-solves with the toroidal-flux profiles, and compares. Differences come from profile interpolation in the
# two coordinates.
def _torflux_ITER_gs(fe_order, Ip_target=15.6E6):
    create_ITER_mesh()
    mygs = setup_ITER_gs(OFT_env(nthreads=-1),fe_order,maxits=100)
    set_ITER_constraints(mygs,Ip_target=Ip_target)
    return mygs


def _torflux_phi_map(mygs):
    '''Phi_N(psi_N) and dPhi_N/dpsi_N from the converged q profile'''
    from scipy.integrate import cumulative_trapezoid
    psi_n, q, _, _, _, _ = mygs.get_q(npsi=400,psi_pad=1.E-4)
    phi = cumulative_trapezoid(q,psi_n,initial=0.0)
    return psi_n, phi/phi[-1], q/phi[-1]


def _torflux_to_phi(x, pmap):
    '''Map psi_N nodes to Phi_N, pinning the end points to 0 and 1 exactly'''
    x_phi = np.interp(x,pmap[0],pmap[1])
    x_phi[0] = 0.0; x_phi[-1] = 1.0
    return x_phi


def _torflux_compare(mygs, ref, its=None):
    '''Relative psi, q, and Ip differences from the psi_N reference solution'''
    psi = mygs.get_psi(False)
    _, q, _, _, _, _ = mygs.get_q(npsi=100,psi_pad=1.E-3)
    stats = mygs.get_stats()
    return {
        'psi_err': np.linalg.norm(psi-ref['psi'])/np.linalg.norm(ref['psi']),
        'q_err': np.max(np.abs(q-ref['q'])/np.abs(ref['q'])),
        'Ip_err': abs(stats['Ip']-ref['Ip'])/ref['Ip'],
        'its': its
    }


def _torflux_ref(mygs, its=None):
    '''psi, q, and Ip of the current solution, for comparison with `_torflux_compare`'''
    _, q, _, _, _, _ = mygs.get_q(npsi=100,psi_pad=1.E-3)
    return {'psi': mygs.get_psi(False), 'q': q, 'Ip': mygs.get_stats()['Ip'], 'its': its}


def run_ITER_torflux_case(fe_order, test_type, mp_q):
    try:
        results = {}
        if test_type == 'bootstrap':
            from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap
            mygs = _torflux_ITER_gs(fe_order,Ip_target=13.0E6)
        else:
            mygs = _torflux_ITER_gs(fe_order)
        init = lambda: mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
        pp = create_power_flux_fun(40,4.0,1.0)
        if test_type == 'ffp':
            # FF' and P' as derivatives with respect to Phi_N, rescaled from psi_N by dPhi_N/dpsi_N
            ffp = create_power_flux_fun(40,1.5,2.0)
            mygs.set_profiles(ffp_prof=ffp,pp_prof=pp)
            init()
            _, its = mygs.solve(return_its=True)
            ref = _torflux_ref(mygs,its)
            pmap = _torflux_phi_map(mygs)
            profs = {}
            for key, prof in (('ffp', ffp), ('pp', pp)):
                y = prof['y']/np.interp(prof['x'],pmap[0],pmap[2])
                profs[key] = {'type': 'linterp', 'x': _torflux_to_phi(prof['x'],pmap), 'y': y, 'coord': 'phi_n'}
            mygs.set_profiles(ffp_prof=profs['ffp'],pp_prof=profs['pp'])
            init()
            _, its = mygs.solve(return_its=True)
            results['phi_n'] = _torflux_compare(mygs,ref,its)
            results['psi_its'] = ref['its']
            # The solver's map must satisfy dPhi_N/dpsi_N = q/Q, with Q the total toroidal flux normalization
            from scipy.integrate import cumulative_trapezoid
            x = np.linspace(0.05,0.95,19)
            phi_x, jac_x = mygs.get_torflux_map(x)
            psi_x, _ = mygs.get_torflux_map(phi_x,inverse=True)
            psi_q, q, _, _, _, _ = mygs.get_q(psi=np.linspace(0.05,0.95,400))
            Qinv = jac_x/np.interp(x,psi_q,q)
            dphi = cumulative_trapezoid(q,psi_q,initial=0.0)*np.mean(Qinv)
            results['jac_err'] = np.max(np.abs(Qinv/np.mean(Qinv)-1.0))
            results['map_err'] = np.max(np.abs(phi_x-phi_x[0]-np.interp(x,psi_q,dphi)))
            results['inv_err'] = np.max(np.abs(psi_x-x))
            # Copies and saved equilibria keep (or rebuild) the map
            prof0 = np.array(mygs.get_profiles(npsi=60)[1:])
            eq_copy = mygs.copy_eq()
            prof_c = np.array(eq_copy.get_profiles(npsi=60)[1:])
            results['copy_err'] = np.max(np.abs(prof_c-prof0))/np.max(np.abs(prof0))
            eq_copy.save_TokaMaker('torflux_eq.h5')
            init()
            mygs.replace_eq(source_file='torflux_eq.h5')
            prof_l = np.array(mygs.get_profiles(npsi=60)[1:])
            results['load_err'] = np.max(np.abs(prof_l-prof0))/np.max(np.abs(prof0))
        elif test_type == 'jphi':
            # Toroidal current density profile given on Phi_N
            x = np.linspace(0.0,1.0,65)
            jphi = {'type': 'jphi-linterp', 'x': x, 'y': create_power_flux_fun(65,2.25,2.5)['y']}
            mygs.set_profiles(ffp_prof=jphi,pp_prof=pp)
            init()
            _, its = mygs.solve(return_its=True)
            ref = _torflux_ref(mygs,its)
            pmap = _torflux_phi_map(mygs)
            # Switching to a bootstrap-split profile on Phi_N keeps F*F', also across save/load
            _, f, fp, _, _ = mygs.get_profiles(npsi=60)
            mygs.set_profiles(ffp_prof=dict(jphi,type='jphi-split-bootstrap',coord='phi_n'))
            mygs.copy_eq().save_TokaMaker('torflux_carry.h5')
            mygs.replace_eq(source_file='torflux_carry.h5')
            _, f1, fp1, _, _ = mygs.get_profiles(npsi=60)
            results['carry_err'] = np.max(np.abs(f1*fp1-f*fp))/np.max(np.abs(f*fp))
            mygs.set_profiles(ffp_prof=dict(jphi,x=_torflux_to_phi(x,pmap),coord='phi_n'),pp_prof=pp)
            init()
            _, its = mygs.solve(return_its=True)
            results['phi_n'] = _torflux_compare(mygs,ref,its)
            results['psi_its'] = ref['its']
        elif test_type == 'bootstrap':
            # Internal bootstrap solve with every input profile on Phi_N
            n = 257
            x = np.linspace(0.0,1.0,n)
            ne, Te = ITER_Hmode_profiles(n)
            jind = create_power_flux_fun(n,2.25,2.5)['y']
            jfix = 1.E5*np.exp(-((x-0.5)/0.15)**2)
            pfix = 0.2*2.0*eC*ne[0]*Te[0]*np.exp(-((x-0.3)/0.15)**2)  # 20% of the axis kinetic pressure [Pa]
            mygs.set_profiles(ffp_prof=create_power_flux_fun(40,1.5,2.0),pp_prof=pp)
            init()
            mygs.solve()
            def run(x_in, coord):
                return mygs.solve_bootstrap(
                    ffp_prof={'type': 'jphi-split-bootstrap', 'x': x_in, 'y': jind},
                    te_prof={'type': 'linterp', 'x': x_in, 'y': Te/1.E3},
                    ne_prof={'type': 'linterp', 'x': x_in, 'y': ne},
                    ti_prof={'type': 'linterp', 'x': x_in, 'y': Te/1.E3},
                    ni_prof={'type': 'linterp', 'x': x_in, 'y': ne},
                    jphi_fixed_prof={'type': 'linterp', 'x': x_in, 'y': jfix},
                    p_fixed_prof={'x': x_in, 'y': pfix},
                    Zeff=1.5, Ip_target=13.0E6, coord=coord)
            res_psi = run(x,'psi_n')
            ref = _torflux_ref(mygs)
            x_phi = _torflux_to_phi(x,_torflux_phi_map(mygs))
            res_phi = run(x_phi,'phi_n')
            results['phi_n'] = _torflux_compare(mygs,ref)
            # Compare the bootstrap current on matching surfaces, away from the edge peak
            j_ref = np.interp(res_phi['psi_n'],res_psi['psi_n'],res_psi['jphi_bs'])
            mask = (res_phi['psi_n'] > 0.05) & (res_phi['psi_n'] < 0.9)
            results['jbs_err'] = np.max(np.abs(res_phi['jphi_bs'][mask]-j_ref[mask]))/np.max(np.abs(j_ref))
            # A fixed current given on Phi_N nodes is returned unchanged on those nodes
            results['jfix_err'] = np.max(np.abs(res_phi['jphi_fixed']-jfix))/np.max(jfix)
            # Output node positions in psi_N agree with the solver's inverse map
            psi_nodes, _ = mygs.get_torflux_map(x_phi,inverse=True)
            results['node_err'] = np.max(np.abs(res_phi['psi_n']-psi_nodes))
            # solve_with_bootstrap() with Phi_N inputs, through the internal solver (the Python solver rejects them)
            res_swb = solve_with_bootstrap(mygs,ne,Te,ne,Te,1.5,13.0E6,inductive_jphi=jind,x=x_phi,coord='phi_n',
                                           jphi_fixed=jfix,p_fixed=pfix)
            psi_nodes, _ = mygs.get_torflux_map(x_phi,inverse=True)
            results['swb_node_err'] = np.max(np.abs(res_swb['psi_n']-psi_nodes))
            results['swb_jfix_err'] = np.max(np.abs(res_swb['jphi_fixed']-jfix))/np.max(jfix)
            # With the reset off, moving part of the fixed current into jphi_saw changes nothing
            res_saw = solve_with_bootstrap(mygs,ne,Te,ne,Te,1.5,13.0E6,inductive_jphi=jind,x=x_phi,coord='phi_n',
                                           jphi_fixed=0.6*jfix,jphi_saw=0.4*jfix,p_fixed=pfix)
            results['swb_jsaw_err'] = np.max(np.abs(res_saw["jphi_saw"]-0.4*jfix))/np.max(jfix)
            results['swb_jsaw_total_err'] = (np.max(np.abs(res_saw['jphi_total']-res_swb['jphi_total']))
                                             /np.max(np.abs(res_swb['jphi_total'])))
            try:
                solve_with_bootstrap(mygs,ne,Te,ne,Te,1.5,13.0E6,inductive_jphi=jind,x=x_phi,
                                     use_python_solve=True,coord='phi_n')
                results['python_rejected'] = False
            except ValueError:
                results['python_rejected'] = True
    except Exception as e:
        print(e)
        mp_q.put(None)
        return
    mp_q.put(results)
    oftpy_dump_cov()


def _torflux_check(res, psi_tol, q_tol, its_ref=None, ip_tol=1.E-4):
    '''Check psi, q, and Ip against the psi_N solution, and that convergence takes at most 5 more iterations'''
    assert res['psi_err'] < psi_tol
    assert res['q_err'] < q_tol
    # Bootstrap I_p comes from the 1-D profile quadrature (no FEM rescale)
    assert res['Ip_err'] < ip_tol
    if its_ref is not None:
        assert res['its'] <= its_ref + 5


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,))
def test_ITER_torflux_ffp(order):
    '''ITER equilibrium with FF' and P' profiles given on normalized toroidal flux matches the same profiles given
    on normalized poloidal flux; the psi_N to Phi_N map is consistent with q and survives copy and save/load'''
    results = mp_run(run_ITER_torflux_case,(order,'ffp'),timeout=120)
    assert results is not None
    _torflux_check(results['phi_n'],4.E-3,4.E-3,results['psi_its'])
    assert results['map_err'] < 1.E-4
    assert results['jac_err'] < 1.E-3
    assert results['inv_err'] < 1.E-10
    assert results['copy_err'] < 1.E-12
    assert results['load_err'] < 5.E-5


@pytest.mark.coverage
@pytest.mark.parametrize("order", (2,))
def test_ITER_torflux_jphi(order):
    '''ITER equilibrium with a toroidal current density profile on normalized toroidal flux matches the psi_N
    solution, and switching profile types on Phi_N keeps F*F' through save/load'''
    results = mp_run(run_ITER_torflux_case,(order,'jphi'),timeout=120)
    assert results is not None
    _torflux_check(results['phi_n'],1.E-3,2.E-3,results['psi_its'])
    assert results['carry_err'] < 1.E-3


@pytest.mark.slow
@pytest.mark.parametrize("order", (2,))
def test_ITER_torflux_bootstrap(order):
    '''Internal bootstrap solves, called directly and through solve_with_bootstrap(), with all profiles on
    normalized toroidal flux match the psi_N solution and return fixed and sawtooth currents on the input nodes;
    the Python bootstrap solver rejects toroidal flux coordinates'''
    results = mp_run(run_ITER_torflux_case,(order,'bootstrap'),timeout=600)
    assert results is not None
    _torflux_check(results['phi_n'],3.E-3,3.E-3,ip_tol=2.E-4)
    assert results['jbs_err'] < 1.E-2
    assert results['node_err'] < 1.E-12
    assert results['swb_node_err'] < 1.E-12
    assert results['jfix_err'] < 1.E-6
    assert results['swb_jfix_err'] < 1.E-6
    assert results['swb_jsaw_err'] < 1.E-6
    assert results['swb_jsaw_total_err'] < 1.E-3
    assert results['python_rejected']

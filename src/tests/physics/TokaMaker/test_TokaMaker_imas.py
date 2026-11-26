#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: IMAS (Data Dictionary 4.x) input and output through imas-python
import json
import pytest
import numpy as np
from tokamaker_testing import (mp_run, validate_dict, ITER_eq_dict, ITER_isoflux_pts, ITER_x_point, solve_ITER_reference,
                               set_ITER_constraints, solve_ITER_eq)
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit._interface import oftpy_dump_cov

imas = pytest.importorskip("imas", reason="imas-python not installed")
from OpenFUSIONToolkit.TokaMaker import imas_io  # noqa: E402


def _trace_averages(mygs, psi_n):
    '''Area, surface area, <|grad psi|^2/B^2>, B_min and B_max of a surface from a dense resampling of `trace_surf`'''
    pts = mygs.trace_surf(psi_n)
    pts = np.vstack((pts, pts[:1]))
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    ss = np.linspace(0.0, s[-1], 20001)
    r = np.interp(ss, s, pts[:, 0])
    z = np.interp(ss, s, pts[:, 1])
    dl = np.diff(ss)
    rm = (r[1:]+r[:-1])/2.0
    zm = (z[1:]+z[:-1])/2.0
    B = np.asarray(mygs.get_field_eval('B').eval(np.column_stack((rm, zm)))).reshape(-1, 3)
    Bm = np.linalg.norm(B, axis=1)
    Bp = np.hypot(B[:, 0], B[:, 2])
    w = dl/Bp
    return {
        'area': 0.5*abs(np.sum(r[:-1]*z[1:]-r[1:]*z[:-1])),
        'surface': 2.0*np.pi*np.sum(rm*dl),
        '<|grad psi|^2/B^2>': np.sum(w*(rm*Bp/Bm)**2)/np.sum(w),
        'B_min': Bm.min(),
        'B_max': Bm.max(),
    }


def run_output_case(fe_order, mp_q):
    mygs, EQ_obj = solve_ITER_reference(1.0, fe_order)
    if mygs is None:
        mp_q.put(None)
        return
    results = {}
    # New flux-surface averages against direct evaluation on traced surfaces
    fsa = mygs.get_fsa(psi=np.r_[0.25, 0.5, 0.9])
    for i, psi_n in enumerate((0.25, 0.5, 0.9)):
        ref = _trace_averages(mygs, psi_n)
        for key, val in ref.items():
            results['fsa_err_{0}_{1}'.format(key, i)] = abs(fsa[key][i]-val)/abs(val)
    # Write, validate, and read back
    eq_ids = imas_io.equilibrium_to_imas(mygs, 1.5)
    pf_ids = imas_io.pf_active_to_imas(mygs, 1.5)
    eq_ids.validate()
    pf_ids.validate()
    imas_io.save_imas('imas_test.nc', [eq_ids, pf_ids])
    eq_ids = imas_io.load_imas('imas_test.nc', 'equilibrium')
    pf_ids = imas_io.load_imas('imas_test.nc', 'pf_active')
    # Older DD 4 version: fields new in 4.1 are skipped, and the file is written with that version
    eq40 = imas_io.equilibrium_to_imas(mygs, 1.5, dd_version='4.0.0')
    pf40 = imas_io.pf_active_to_imas(mygs, 1.5, dd_version='4.0.0')
    eq40.validate()
    imas_io.save_imas('imas_test_dd40.nc', [eq40, pf40])
    results['dd40_version'] = imas.util.get_data_dictionary_version(imas_io.load_imas('imas_test_dd40.nc', 'equilibrium'))
    # Existing equilibrium IDS of another DD version: the pf_active companion follows it
    mygs.save_imas('imas_test_dd40_method.nc', 1.5, ids=imas.IDSFactory('4.0.0').equilibrium())
    results['dd40_method_version'] = imas.util.get_data_dictionary_version(
        imas_io.load_imas('imas_test_dd40_method.nc', 'pf_active'))
    # The vacuum field is read from the solved field, not the cached F0 property (stale after a fitF0 reconstruction)
    F0_save = EQ_obj._F0
    EQ_obj._F0 = 0.5*F0_save
    eq_stale = imas_io.equilibrium_to_imas(EQ_obj, 1.5, npsi=16, nr=16, nz=16)
    EQ_obj._F0 = F0_save
    results['b0_err'] = abs(float(eq_stale.vacuum_toroidal_field.b0[0])*float(eq_stale.vacuum_toroidal_field.r0)-F0_save)
    # Convenience method writes the same equilibrium
    mygs.save_imas('imas_test_method.nc', 1.5)
    results['method_ip_err'] = abs(float(imas_io.load_imas('imas_test_method.nc', 'equilibrium').time_slice[0].global_quantities.ip)
                                   - float(eq_ids.time_slice[0].global_quantities.ip))
    ts = eq_ids.time_slice[0]
    gq = ts.global_quantities
    p1 = ts.profiles_1d
    p2 = ts.profiles_2d[0]
    stats = EQ_obj.get_stats(li_normalization='iter')
    results['time'] = float(eq_ids.time[0])
    results['Ip'] = float(gq.ip)
    results['q_95'] = float(gq.q_95)
    results['vol'] = float(gq.volume)
    results['l_i'] = float(gq.li_3)
    results['beta_tor'] = 100.0*float(gq.beta_tor)
    results['W_MHD'] = float(gq.energy_mhd)
    results['li_iter'] = stats['l_i']
    # Profiles against get_profiles/get_fsa on the same grid (COCOS 7 -> 17: psi x 2pi, d/dpsi / 2pi)
    psi_n = np.asarray(p1.psi_norm)
    _, f, fp, p, pp = EQ_obj.get_profiles(psi=psi_n)
    q_ref = EQ_obj.get_fsa(psi=psi_n[1:])['q']
    results['prof_err_psi'] = np.max(np.abs(np.asarray(p1.psi)-2.0*np.pi*EQ_obj.psinorm_to_absolute(psi_n)))/abs(gq.psi_axis)
    results['prof_err_ffp'] = np.max(np.abs(np.asarray(p1.f_df_dpsi)*2.0*np.pi-f*fp))/np.max(np.abs(f*fp))
    results['prof_err_pp'] = np.max(np.abs(np.asarray(p1.dpressure_dpsi)*2.0*np.pi-pp))/np.max(np.abs(pp))
    results['prof_err_q'] = np.max(np.abs(np.asarray(p1.q)[1:]-q_ref)/np.abs(q_ref))
    # Exact plasma current from the IMAS j_phi (dI/dA) to the boundary surface, within the LCFS padding
    results['Ip_jphi'] = float(np.trapezoid(p1.j_phi, p1.area)) if hasattr(np, 'trapezoid') else float(np.trapz(p1.j_phi, p1.area))
    results['all_finite'] = all(np.all(np.isfinite(np.asarray(getattr(p1, key)))) for key in (
        'psi', 'phi', 'q', 'f', 'pressure', 'f_df_dpsi', 'dpressure_dpsi', 'j_phi', 'j_parallel', 'volume', 'area',
        'surface', 'gm1', 'gm2', 'gm3', 'gm4', 'gm5', 'gm6', 'gm7', 'gm8', 'gm9', 'rho_tor', 'r_inboard', 'r_outboard',
        'elongation', 'triangularity_upper', 'triangularity_lower', 'trapped_fraction', 'b_field_min', 'b_field_max'))
    # COCOS 17: B_R = -dpsi/dZ/(2 pi R) and B_Z = dpsi/dR/(2 pi R), with B from TokaMaker
    r_grid = np.asarray(p2.grid.dim1)
    z_grid = np.asarray(p2.grid.dim2)
    dpsi_dr, dpsi_dz = np.gradient(np.asarray(p2.psi), r_grid, z_grid)
    r2d = np.asarray(p2.r)
    # (RMS relative error: finite differences are inaccurate where the current density jumps; a sign or 2 pi error is O(1))
    inner = (slice(2, -2), slice(2, -2))
    bscale = np.sqrt(np.mean((np.asarray(p2.b_field_r)**2+np.asarray(p2.b_field_z)**2)[inner]))
    results['cocos_err_br'] = np.sqrt(np.mean((-dpsi_dz/(2.0*np.pi*r2d)-np.asarray(p2.b_field_r))[inner]**2))/bscale
    results['cocos_err_bz'] = np.sqrt(np.mean((dpsi_dr/(2.0*np.pi*r2d)-np.asarray(p2.b_field_z))[inner]**2))/bscale
    # Positive Ip: COCOS 17 flux decreases from the axis, so the O-point is a maximum
    results['psi_decreasing'] = bool(gq.psi_axis > gq.psi_boundary)
    results['psi_boundary_err'] = abs(float(gq.psi_boundary)-float(ts.boundary.psi))
    results['o_point_type'] = int(ts.contour_tree.node[0].critical_type.value)
    results['n_xpoints'] = len(ts.contour_tree.node)-1
    # Coil currents
    currents, _ = EQ_obj.get_coil_currents()
    results['coil_err'] = max(abs(float(coil.current.data[0])-currents[str(coil.name)]) for coil in pf_ids.coil)/max(
        abs(val) for val in currents.values())
    mp_q.put([results])
    oftpy_dump_cov()


@pytest.mark.coverage
@pytest.mark.parametrize("fe_order", (2, 3))
def test_imas_output(fe_order):
    '''Equilibrium and pf_active IDSs of the ITER reference equilibrium: written, validated, and read back from netCDF;
    globals and profiles match TokaMaker, the field satisfies COCOS 17, and the new flux-surface averages match
    direct evaluation on traced surfaces'''
    results = mp_run(run_output_case, (fe_order,), timeout=120)
    assert results is not None
    res = results[0]
    exp = {key: ITER_eq_dict[key] for key in ('Ip', 'q_95', 'vol', 'W_MHD', 'beta_tor')}
    exp['l_i'] = res['li_iter']
    exp['time'] = 1.5
    exp['Ip_jphi'] = ITER_eq_dict['Ip']
    assert validate_dict(results, exp, {'Ip_jphi': 2.E-3, 'q_95': 2.E-2})
    for key, val in res.items():
        if key.startswith('fsa_err_'):
            assert val < 1.E-4, key
    for key in ('prof_err_psi', 'prof_err_ffp', 'prof_err_pp', 'coil_err'):
        assert res[key] < 1.E-10, key
    assert res['prof_err_q'] < 1.E-10
    assert res['cocos_err_br'] < 1.E-2
    assert res['cocos_err_bz'] < 1.E-2
    assert res['all_finite']
    assert res['method_ip_err'] == 0.0
    assert res['dd40_version'] == '4.0.0'
    assert res['dd40_method_version'] == '4.0.0'
    assert res['b0_err'] < 1.E-10
    assert res['psi_boundary_err'] == 0.0
    assert res['psi_decreasing']
    assert res['o_point_type'] == 2
    assert res['n_xpoints'] >= 1


def _ITER_machine_description():
    '''wall, pf_active and pf_passive IDSs equivalent to the ITER test mesh geometry (ITER_geom.json)'''
    with open('ITER_geom.json', 'r') as fid:
        geom = json.load(fid)
    factory = imas.IDSFactory()
    wall = factory.wall()
    wall.ids_properties.homogeneous_time = 2
    wall.description_2d.resize(1)
    desc = wall.description_2d[0]
    desc.limiter.unit.resize(1)
    lim = np.asarray(geom['limiter'])
    desc.limiter.unit[0].outline.r = lim[:, 0]
    desc.limiter.unit[0].outline.z = lim[:, 1]
    desc.vessel.unit.resize(2)
    for unit, key, name in zip(desc.vessel.unit, ('inner_vv', 'outer_vv'), ('VV1', 'VV2')):
        unit.name = name
        inner, outer = np.asarray(geom[key][0]), np.asarray(geom[key][1])
        unit.annular.outline_inner.r = inner[:, 0]
        unit.annular.outline_inner.z = inner[:, 1]
        unit.annular.outline_outer.r = outer[:, 0]
        unit.annular.outline_outer.z = outer[:, 1]
        unit.annular.resistivity = 6.9E-7
    pf = factory.pf_active()
    pf.ids_properties.homogeneous_time = 2
    coils = [(key, [coil]) for key, coil in geom['coils'].items() if not key.startswith('VS')]
    coils.append(('VS', [geom['coils']['VSU'], geom['coils']['VSL']]))
    pf.coil.resize(len(coils)+1)
    for coil, (name, elements) in zip(pf.coil, coils):
        coil.name = name
        coil.element.resize(len(elements))
        for j, (element, rect) in enumerate(zip(coil.element, elements)):
            element.turns_with_sign = 1.0 if j == 0 else -1.0
            element.geometry.geometry_type = 2
            element.geometry.rectangle.r = rect['rc']
            element.geometry.rectangle.z = rect['zc']
            element.geometry.rectangle.width = rect['w']
            element.geometry.rectangle.height = rect['h']
    # A distant round conductor with a hole (annulus geometry)
    ring = pf.coil[len(coils)]
    ring.name = 'Ring'
    ring.element.resize(1)
    ring.element[0].turns_with_sign = 1.0
    ring.element[0].geometry.geometry_type = 5
    ring.element[0].geometry.annulus.r = 11.0
    ring.element[0].geometry.annulus.z = 5.0
    ring.element[0].geometry.annulus.radius_inner = 0.15
    ring.element[0].geometry.annulus.radius_outer = 0.3
    # One passive loop, given as an oblique (parallelogram) element
    pp = factory.pf_passive()
    pp.ids_properties.homogeneous_time = 2
    pp.loop.resize(1)
    pp.loop[0].name = 'TRIANGULAR_SUPPORT'
    pp.loop[0].resistivity = 1.E-6
    pp.loop[0].element.resize(1)
    obl = pp.loop[0].element[0].geometry
    obl.geometry_type = 3
    obl.oblique.r = 10.5
    obl.oblique.z = -1.0
    obl.oblique.length_alpha = 0.3
    obl.oblique.length_beta = 0.2
    obl.oblique.alpha = 0.2
    obl.oblique.beta = 0.3
    return wall, pf, pp


def run_machine_description_case(mp_q):
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import save_gs_mesh, load_gs_mesh
    wall, pf, pp = _ITER_machine_description()
    domain = imas_io.domain_from_imas(wall, pf, {'plasma': 0.15, 'vacuum': 0.3, 'boundary': 0.6, 'coil': 0.2, 'conductor': 0.3},
                                      pf_passive=pp)
    mesh_pts, mesh_lc, mesh_reg = domain.build_mesh()
    coil_dict = domain.get_coils()
    cond_dict = domain.get_conductors()
    save_gs_mesh(mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict, 'ITER_imas_mesh.h5')
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_imas_mesh.h5')
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.setup_model(mesh_pts, mesh_lc, reg=mesh_reg, cond_dict=cond_dict, coil_dict=coil_dict, order=2, F0=5.3*6.2)
    # Meshed cross-section of the ring (annulus with a hole)
    tri = mesh_pts[mesh_lc]
    tri_area = 0.5*np.abs((tri[:, 1, 0]-tri[:, 0, 0])*(tri[:, 2, 1]-tri[:, 0, 1])-(tri[:, 2, 0]-tri[:, 0, 0])*(tri[:, 1, 1]-tri[:, 0, 1]))
    ring_area = np.sum(tri_area[mesh_reg == coil_dict['RING']['reg_id']])
    results = {
        'ring_area_err': abs(ring_area-np.pi*(0.3**2-0.15**2))/(np.pi*(0.3**2-0.15**2)),
        'n_coil_sets': len(mygs.coil_sets),
        'vs_turns': sorted(sub['nturns'] for sub in mygs.coil_sets['VS']['sub_coils']),
        'conductors': sorted(name for name, cond in cond_dict.items() if 'eta' in cond),
    }
    set_ITER_constraints(mygs)
    solve_ITER_eq(mygs)
    results.update(mygs.get_stats(li_normalization='ITER'))
    # Write currents back into the machine description IDSs
    pf = imas_io.pf_active_to_imas(mygs, 0.0, ids=pf)
    currents, _ = mygs.get_coil_currents()
    results['coil_err'] = max(abs(float(coil.current.data[0])-currents[str(coil.name).upper()]) for coil in pf.coil)
    pp = imas_io.pf_passive_to_imas(mygs, 0.0, ids=pp)
    results['passive_current'] = float(pp.loop[0].current[0])
    # Writeback keeps other signals on their own times: a homogeneous IDS with another time base is
    # rejected unchanged, and a heterogeneous one gets times on the new currents only
    pf_hom = imas.IDSFactory().pf_active()
    pf_hom.ids_properties.homogeneous_time = 1
    pf_hom.time = [0.0, 1.0]
    pf_hom.coil.resize(len(pf.coil))
    for coil, src in zip(pf_hom.coil, pf.coil):
        coil.name = src.name
        coil.voltage.data = [1.0, 2.0]
    try:
        imas_io.pf_active_to_imas(mygs, 0.0, ids=pf_hom)
        results['hom_rejected'] = False
    except ValueError:
        results['hom_rejected'] = (len(pf_hom.time) == 2) and (len(pf_hom.coil[0].current.data) == 0)
    pf_old = imas.IDSFactory().pf_active()
    pf_old.ids_properties.homogeneous_time = 0
    pf_old.coil.resize(len(pf.coil))
    for coil, src in zip(pf_old.coil, pf.coil):
        coil.name = src.name
        coil.current.data = [1.0, 2.0]
        coil.current.data_error_upper = [0.1, 0.1]
        coil.current.time = [0.0, 1.0]
    imas_io.pf_active_to_imas(mygs, 0.0, ids=pf_old)
    pf_old.validate()
    results['old_errors_cleared'] = len(pf_old.coil[0].current.data_error_upper) == 0
    # A passive IDS written fresh from the mesh matches its own (exact) loop names on writeback
    pp_fresh = imas_io.pf_passive_to_imas(mygs, 0.0)
    imas_io.pf_passive_to_imas(mygs, 0.0, ids=pp_fresh, psis=[2.0*mygs.get_psi(False)])
    pp_fresh.validate()
    results['fresh_loops_updated'] = all(len(loop.current) == 1 for loop in pp_fresh.loop)
    pf_het = imas.IDSFactory().pf_active()
    pf_het.ids_properties.homogeneous_time = 0
    pf_het.coil.resize(len(pf.coil))
    for coil, src in zip(pf_het.coil, pf.coil):
        coil.name = src.name
        coil.voltage.data = [1.0, 2.0]
        coil.voltage.time = [0.0, 1.0]
    imas_io.pf_active_to_imas(mygs, 0.5, ids=pf_het)
    pf_het.validate()
    results['het_times'] = (list(pf_het.coil[0].current.time) == [0.5]) and (list(pf_het.coil[0].voltage.time) == [0.0, 1.0])
    # The updated machine description is valid and round-trips through a file
    pf.validate()
    pp.validate()
    imas_io.save_imas('ITER_imas_md.nc', [pf, pp])
    pf_read = imas_io.load_imas('ITER_imas_md.nc', 'pf_active')
    results['ring_current'] = float(pf_read.coil[len(pf_read.coil)-1].current.data[0])
    results['ring_current_ref'] = currents['RING']
    mp_q.put([results])
    oftpy_dump_cov()


@pytest.mark.coverage
def test_imas_machine_description():
    '''ITER mesh built from wall, pf_active (incl. a two-element VS coil) and pf_passive (an oblique element) IDSs;
    the reference equilibrium solved on it matches the hand-built ITER mesh, and currents are written back'''
    results = mp_run(run_machine_description_case, (), timeout=180)
    assert results is not None
    res = results[0]
    assert res['n_coil_sets'] == 14
    assert res['ring_area_err'] < 2.E-2
    assert res['ring_current'] == res['ring_current_ref']
    assert res['hom_rejected']
    assert res['old_errors_cleared']
    assert res['fresh_loops_updated']
    assert res['het_times']
    assert res['vs_turns'] == [-1.0, 1.0]
    assert res['conductors'] == ['TRIANGULAR_SUPPORT', 'VV1', 'VV2']
    exp = {key: ITER_eq_dict[key] for key in ('Ip', 'q_95', 'vol', 'kappa', 'l_i', 'W_MHD')}
    assert validate_dict(results, exp, {'q_95': 2.E-2, 'W_MHD': 2.E-2})
    assert res['coil_err'] == 0.0
    # No induced current in a static equilibrium (up to discretization error)
    assert abs(res['passive_current']) < 1.E-4*ITER_eq_dict['Ip']


def run_profiles_case(mp_q):
    mygs, EQ_obj = solve_ITER_reference(1.0, 2)
    if mygs is None:
        mp_q.put(None)
        return
    eq_ids = imas_io.equilibrium_to_imas(mygs, 0.0)
    pf_ids = imas_io.pf_active_to_imas(mygs, 0.0)
    ref = mygs.get_stats(li_normalization='ITER')
    # Re-solve from the IDS profiles, targets and boundary
    inputs = imas_io.profiles_from_imas(eq_ids)
    currents, _, _ = imas_io.coil_currents_from_imas(pf_ids, mygs)
    set_ITER_constraints(mygs)
    mygs.set_profiles(**inputs['profiles'])
    mygs.set_targets(**inputs['targets'])
    mygs.set_isoflux(np.vstack((inputs['isoflux'][::4], ITER_x_point)))
    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
    mygs.solve()
    results = mygs.get_stats(li_normalization='ITER')
    results['F0_err'] = abs(inputs['F0']-5.3*6.2)
    results['coil_keys'] = sorted(currents.keys()) == sorted(mygs.coil_sets.keys())
    results['ref'] = {key: ref[key] for key in ('Ip', 'q_95', 'vol', 'l_i', 'P_ax')}
    # Kinetic profiles from core_profiles, on rho_tor_norm
    cp = imas.IDSFactory().core_profiles()
    cp.ids_properties.homogeneous_time = 1
    cp.time = [0.0]
    cp.profiles_1d.resize(1)
    rho = np.linspace(0.0, 1.0, 21)
    cp.profiles_1d[0].grid.rho_tor_norm = rho
    cp.profiles_1d[0].electrons.temperature = 1.E4*(1.0-rho**2)+100.0
    cp.profiles_1d[0].electrons.density = 1.E20*(1.0-0.5*rho**2)
    cp.profiles_1d[0].t_i_average = 8.E3*(1.0-rho**2)+100.0
    cp.profiles_1d[0].ion.resize(1)
    cp.profiles_1d[0].ion[0].density = 0.9E20*(1.0-0.5*rho**2)
    cp.profiles_1d[0].zeff = 1.5*np.ones_like(rho)
    kinetic = imas_io.profiles_from_imas(cp)
    mygs.set_kinetic_profiles(**kinetic)
    results['kinetic_keys'] = sorted(kinetic.keys())
    results['te_axis_keV'] = float(kinetic['te_prof']['y'][0])
    results['kinetic_x'] = float(kinetic['te_prof']['x'][10])
    mp_q.put([results])
    oftpy_dump_cov()


@pytest.mark.coverage
def test_imas_profiles():
    '''Re-solving the ITER case from the profiles, targets and boundary of its own equilibrium IDS reproduces it;
    core_profiles kinetic profiles are accepted on normalized toroidal flux in keV'''
    results = mp_run(run_profiles_case, (), timeout=180)
    assert results is not None
    res = results[0]
    assert validate_dict(results, res['ref'], {'q_95': 2.E-2, 'l_i': 2.E-2})
    assert res['F0_err'] < 1.E-10
    assert res['coil_keys']
    assert res['kinetic_keys'] == ['Zeff', 'ne_prof', 'ni_prof', 'te_prof', 'ti_prof']
    assert abs(res['te_axis_keV']-10.1) < 1.E-10
    assert abs(res['kinetic_x']-0.25) < 1.E-12


def test_imas_pulse_schedule():
    '''Shape and current targets interpolated in time from a pulse_schedule IDS'''
    ps = imas.IDSFactory().pulse_schedule()
    ps.ids_properties.homogeneous_time = 1
    ps.time = [0.0, 2.0]
    pc = ps.position_control
    pc.boundary_outline.resize(ITER_isoflux_pts.shape[0])
    for node, pt in zip(pc.boundary_outline, ITER_isoflux_pts):
        node.r.reference = [pt[0], pt[0]+0.2]
        node.z.reference = [pt[1], pt[1]]
    pc.x_point.resize(1)
    pc.x_point[0].r.reference = [ITER_x_point[0, 0], ITER_x_point[0, 0]]
    pc.x_point[0].z.reference = [ITER_x_point[0, 1], ITER_x_point[0, 1]]
    pc.magnetic_axis.r.reference = [6.2, 6.4]
    ps.flux_control.ip.reference = [10.E6, 15.E6]
    targets = imas_io.targets_from_pulse_schedule(ps, time=1.0)
    assert np.allclose(targets['isoflux'][:, 0], ITER_isoflux_pts[:, 0]+0.1)
    assert np.allclose(targets['isoflux'][:, 1], ITER_isoflux_pts[:, 1])
    assert np.allclose(targets['saddles'], ITER_x_point)
    assert abs(targets['R0']-6.3) < 1.E-12
    assert abs(targets['Ip']-12.5E6) < 1.E-6
    assert 'Z0' not in targets


def run_magnetics_case(mp_q):
    from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
    mygs, EQ_obj = solve_ITER_reference(1.0, 2)
    if mygs is None:
        mp_q.put(None)
        return
    # Synthetic magnetics on the inner vessel wall, with probes at assorted poloidal angles (DD 4.1 normal convention)
    with open('ITER_geom.json', 'r') as fid:
        contour = np.asarray(json.load(fid)['inner_vv'][0])[::10]
    psi_eval = mygs.get_field_eval('PSI')
    B_eval = mygs.get_field_eval('B')
    mag = imas.IDSFactory().magnetics()
    mag.ids_properties.homogeneous_time = 1
    mag.time = [0.0]
    mag.flux_loop.resize(contour.shape[0])
    mag.b_field_pol_probe.resize(contour.shape[0])
    for i, pt in enumerate(contour):
        loop = mag.flux_loop[i]
        loop.position.resize(1)
        loop.position[0].r = pt[0]
        loop.position[0].z = pt[1]
        loop.flux.data = [2.0*np.pi*np.asarray(psi_eval.eval(pt)).reshape(-1)[0]]
        probe = mag.b_field_pol_probe[i]
        alpha = 0.7*i
        B = np.asarray(B_eval.eval(pt)).reshape(-1)
        probe.position.r = pt[0]
        probe.position.z = pt[1]
        probe.poloidal_angle = alpha
        probe.toroidal_angle = 0.0
        probe.field.data = [B[0]*np.cos(alpha)-B[2]*np.sin(alpha)]
    # Loops that must be skipped: an invalid sample on the first loop, and a saddle loop
    mag.flux_loop[0].flux.validity_timed = [-2]
    mag.flux_loop.resize(contour.shape[0]+1, keep=True)
    saddle = mag.flux_loop[contour.shape[0]]
    saddle.type.index = 2
    saddle.position.resize(2)
    saddle.flux.data = [1.0]
    # Also skipped: a typed poloidal loop with two positions, and one without Z
    mag.flux_loop.resize(contour.shape[0]+3, keep=True)
    multi = mag.flux_loop[contour.shape[0]+1]
    multi.type.index = 1
    multi.position.resize(2)
    for pos in multi.position:
        pos.r = 8.0
        pos.z = 0.0
    multi.flux.data = [1.0]
    no_z = mag.flux_loop[contour.shape[0]+2]
    no_z.position.resize(1)
    no_z.position[0].r = 8.0
    no_z.flux.data = [1.0]
    # Probes that must be skipped: a differential probe and one without a poloidal angle
    nprobe = contour.shape[0]
    mag.b_field_pol_probe.resize(nprobe+2, keep=True)
    diff = mag.b_field_pol_probe[nprobe]
    diff.type.index = 6
    diff.indices_differential = [1, 2]
    diff.field.data = [1.0]
    no_angle = mag.b_field_pol_probe[nprobe+1]
    no_angle.position.r = 8.0
    no_angle.position.z = 0.0
    no_angle.field.data = [1.0]
    mag.ip.resize(1)
    mag.ip[0].data = [EQ_obj.get_stats()['Ip']]
    pf_ids = imas_io.pf_active_to_imas(mygs, 0.0)
    myrecon = reconstruction(mygs)
    results = {}
    # Incomplete pf_active: rejected before any constraint is added
    pf_short = imas.IDSFactory().pf_active()
    pf_short.ids_properties.homogeneous_time = 1
    pf_short.time = [0.0]
    pf_short.coil.resize(1)
    pf_short.coil[0].name = pf_ids.coil[0].name
    pf_short.coil[0].current.data = [1.0]
    try:
        imas_io.load_magnetics(myrecon, mag, pf_active=pf_short, abs_err={'b_field_pol_probe': 1.E-3, 'coil': 1.0})
        results['short_pf_rejected'] = False
    except KeyError:
        results['short_pf_rejected'] = len(myrecon._flux_loops) == 0
    # A zero probe reading without an error floor has no usable uncertainty
    mag_zero = imas.IDSFactory().magnetics()
    mag_zero.ids_properties.homogeneous_time = 1
    mag_zero.time = [0.0]
    mag_zero.b_field_pol_probe.resize(1)
    mag_zero.b_field_pol_probe[0].position.r = 8.0
    mag_zero.b_field_pol_probe[0].position.z = 0.0
    mag_zero.b_field_pol_probe[0].poloidal_angle = 0.0
    mag_zero.b_field_pol_probe[0].field.data = [0.0]
    try:
        imas_io.load_magnetics(myrecon, mag_zero)
        results['zero_err_rejected'] = False
    except ValueError:
        results['zero_err_rejected'] = len(myrecon._flux_loops) == 0
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        ncons = imas_io.load_magnetics(myrecon, mag, pf_active=pf_ids, abs_err={'b_field_pol_probe': 1.E-3, 'coil': 1.0})
    messages = ' '.join(str(w.message) for w in caught)
    results['skip_warnings'] = ('flux_loop({0})'.format(contour.shape[0]+2) in messages) \
        and ('flux_loop({0})'.format(contour.shape[0]+3) in messages) and ('b_field_pol_probe({0})'.format(nprobe+1) in messages) \
        and ('b_field_pol_probe({0})'.format(nprobe+2) in messages)
    myrecon.setup_constraints()
    error_mat = myrecon.eval_error()
    eq_ids = imas_io.equilibrium_to_imas(mygs, 0.0)
    # A constraint from another code is replaced, so the structure describes this fit only
    eq_ids.time_slice[0].constraints.mse_polarization_angle.resize(1)
    eq_ids.time_slice[0].constraints.mse_polarization_angle[0].measured = 0.1
    imas_io.recon_constraints_to_imas(myrecon, eq_ids, error_mat=error_mat, freedom_degrees_n=ncons-3)
    cons = eq_ids.time_slice[0].constraints
    results['mse_cleared'] = len(cons.mse_polarization_angle) == 0
    # Re-exporting the equilibrium into the same IDS keeps the constraints; other slice times are rejected
    imas_io.equilibrium_to_imas(mygs, 0.0, ids=eq_ids, npsi=16, nr=16, nz=16)
    results['constraints_kept'] = bool(cons.ip.measured.has_value) and (len(eq_ids.time_slice[0].profiles_1d.psi) == 16)
    try:
        imas_io.equilibrium_to_imas(mygs, 1.0, ids=eq_ids, npsi=16, nr=16, nz=16)
        results['other_times_rejected'] = False
    except ValueError:
        results['other_times_rejected'] = True
    results.update({
        'ncons': ncons,
        'n_expected': 2*contour.shape[0]-1+1+len(mygs.coil_sets),
        'max_residual': float(np.max(np.abs(error_mat[:, 0]))),
        'chi2_red': float(cons.chi_squared_reduced),
        'n_flux': len(cons.flux_loop),
        'n_flux_src': len(mag.flux_loop),
        'n_probe': len(cons.b_field_pol_probe),
        'n_probe_src': len(mag.b_field_pol_probe),
        'slot0_empty': not cons.flux_loop[0].measured.has_value,
        'skipped_probes_empty': not any(cons.b_field_pol_probe[k].measured.has_value for k in (nprobe, nprobe+1)),
        'slot1_source': str(cons.flux_loop[1].source),
        'coil_sources_ok': all(str(node.source) == 'pf_active/coil({0})/current'.format(i+1) for i, node in enumerate(cons.pf_current)),
        'weight': float(cons.flux_loop[1].weight),
        'err': float(cons.flux_loop[1].measured_error_upper),
        'ip_meas': float(cons.ip.measured),
    })
    # A manual constraint is appended after the source slots; the same source twice is an error
    myrecon.add_flux_loop(np.array([8.5, 0.0]), 1.0, 0.1)
    myrecon.setup_constraints()
    imas_io.recon_constraints_to_imas(myrecon, eq_ids)
    results['manual_appended'] = (len(cons.flux_loop) == len(mag.flux_loop)+1) and (str(cons.flux_loop[1].source) == 'magnetics/flux_loop(2)/flux') \
        and (not cons.flux_loop[len(mag.flux_loop)].source.has_value)
    dup = reconstruction(mygs)
    imas_io.load_magnetics(dup, mag, abs_err={'b_field_pol_probe': 1.E-3})
    imas_io.load_magnetics(dup, mag, abs_err={'b_field_pol_probe': 1.E-3})
    try:
        imas_io.recon_constraints_to_imas(dup, eq_ids, error_mat=np.zeros((2*ncons, 4)))
        results['duplicate_rejected'] = False
    except ValueError:
        # Rejected before the previous fit is touched
        results['duplicate_rejected'] = bool(cons.ip.measured.has_value) and (len(cons.flux_loop) == len(mag.flux_loop)+1)
    # Updating an existing IDS: homogeneous slices are matched by the root time, other slice data
    # (here a 2D coordinate system) is replaced, and the constraints are kept
    target = imas.IDSFactory().equilibrium()
    target.ids_properties.homogeneous_time = 1
    target.time = [0.0]
    target.vacuum_toroidal_field.r0 = 6.2
    target.vacuum_toroidal_field.b0 = [5.3]
    target.time_slice.resize(1)
    target.time_slice[0].coordinate_system.grid.dim1 = [0.0, 1.0]
    target.time_slice[0].coordinate_system.grid.dim2 = [0.0, 1.0]
    target.time_slice[0].coordinate_system.r = np.ones((2, 2))
    target.time_slice[0].constraints.ip.measured = 1.0
    target.validate()
    imas_io.equilibrium_to_imas(mygs, 0.0, ids=target, npsi=16, nr=16, nz=16)
    target.validate()
    results['update_root_time'] = (len(target.time_slice[0].coordinate_system.r) == 0) \
        and (float(target.time_slice[0].constraints.ip.measured) == 1.0) and (len(target.time_slice[0].profiles_1d.psi) == 16)
    target.time_slice[0].time = 2.0  # A redundant local time does not stand in for the root time
    try:
        imas_io.equilibrium_to_imas(mygs, 2.0, ids=target, npsi=16, nr=16, nz=16)
        results['local_time_ignored'] = False
    except ValueError:
        results['local_time_ignored'] = list(target.time) == [0.0]
    # Rewriting the slice with an empty fit leaves nothing from the previous one
    imas_io.recon_constraints_to_imas(reconstruction(mygs), eq_ids, error_mat=np.zeros((0, 4)))
    results['rewrite_clean'] = (len(cons.flux_loop) == 0) and (not cons.ip.measured.has_value) \
        and (not cons.chi_squared_reduced.has_value) and (not cons.freedom_degrees_n.has_value) and (int(cons.constraints_n) == 0)
    mp_q.put([results])
    oftpy_dump_cov()


@pytest.mark.coverage
def test_imas_magnetics():
    '''Synthetic magnetics and pf_active IDSs of the ITER equilibrium, loaded as reconstruction constraints, are matched
    by that same equilibrium (checks flux loop and probe conventions), and the constraints are written back to IMAS'''
    results = mp_run(run_magnetics_case, (), timeout=120)
    assert results is not None
    res = results[0]
    assert res['short_pf_rejected']
    assert res['zero_err_rejected']
    assert res['ncons'] == res['n_expected']
    assert res['max_residual'] < 1.E-3
    assert res['chi2_red'] < 1.E-6
    assert res['skip_warnings']
    assert res['n_flux'] == res['n_flux_src']
    assert res['n_probe'] == res['n_probe_src']
    assert res['skipped_probes_empty']
    assert res['slot0_empty']
    assert res['slot1_source'] == 'magnetics/flux_loop(2)/flux'
    assert res['coil_sources_ok']
    assert res['rewrite_clean']
    assert res['mse_cleared']
    assert res['constraints_kept']
    assert res['other_times_rejected']
    assert res['manual_appended']
    assert res['duplicate_rejected']
    assert res['update_root_time']
    assert res['local_time_ignored']
    assert res['weight'] == 1.0
    assert res['err'] > 0.0
    assert abs(res['ip_meas']-ITER_eq_dict['Ip'])/ITER_eq_dict['Ip'] < 1.E-2


def run_torax_case(mp_q):
    from torax._src.geometry import imas as torax_imas
    mygs, EQ_obj = solve_ITER_reference(1.0, 2)
    if mygs is None:
        mp_q.put(None)
        return
    eq_ids = imas_io.equilibrium_to_imas(mygs, 0.0)
    geo = torax_imas.IMASConfig(equilibrium_object=eq_ids)._build_geometry()
    p1 = eq_ids.time_slice[0].profiles_1d
    mp_q.put([{
        'vol_err': abs(float(np.asarray(geo.volume_face)[-1])-p1.volume[-1])/p1.volume[-1],
        'phi_err': abs(float(np.asarray(geo.Phi_face)[-1])-abs(p1.phi[-1]))/abs(p1.phi[-1]),
    }])
    oftpy_dump_cov()


def test_imas_torax():
    '''TORAX's IMAS geometry loader accepts the TokaMaker equilibrium IDS'''
    pytest.importorskip("torax", reason="TORAX not installed")
    results = mp_run(run_torax_case, (), timeout=300)
    assert results is not None
    assert results[0]['vol_err'] < 1.E-3
    assert results[0]['phi_err'] < 1.E-3


def test_imas_dd3_rejected():
    '''IDSs from Data Dictionary 3.x are rejected'''
    eq3 = imas.IDSFactory("3.42.0").equilibrium()
    with pytest.raises(ValueError, match="Data Dictionary 4"):
        imas_io.check_dd4(eq3)
    with pytest.raises(ValueError, match="Data Dictionary 4"):
        imas_io.profiles_from_imas(eq3)


def test_imas_heterogeneous_time():
    '''core_profiles with per-slice times (heterogeneous mode) selects the slice nearest the requested time'''
    cp = imas.IDSFactory().core_profiles()
    cp.ids_properties.homogeneous_time = 0
    cp.profiles_1d.resize(2)
    for sl, t, te in zip(cp.profiles_1d, (0.0, 1.0), (1.E3, 2.E3)):
        sl.time = t
        sl.grid.rho_tor_norm = [0.0, 1.0]
        sl.electrons.temperature = [te, te]
    assert imas_io.profiles_from_imas(cp, time=1.0)['te_prof']['y'][0] == 2.0
    assert imas_io.profiles_from_imas(cp, time=0.1)['te_prof']['y'][0] == 1.0


def _simple_equilibrium_ids(ip=1.E6, ffp=None):
    eq = imas.IDSFactory().equilibrium()
    eq.ids_properties.homogeneous_time = 1
    eq.time = [0.0]
    eq.vacuum_toroidal_field.r0 = 2.0
    eq.vacuum_toroidal_field.b0 = [3.0]
    eq.time_slice.resize(1)
    ts = eq.time_slice[0]
    ts.global_quantities.ip = ip
    ts.profiles_1d.psi_norm = [0.0, 0.5, 1.0]
    ts.profiles_1d.f_df_dpsi = [0.0, 0.0, 0.0] if ffp is None else ffp
    ts.profiles_1d.dpressure_dpsi = [-2.0, -1.0, 0.0]
    ts.profiles_1d.pressure = [1.E5, 5.E4, 0.0]
    return eq


def test_imas_profile_inputs():
    '''Zero FF' profiles stay zero, F0 comes from vacuum_toroidal_field, and negative Ip is rejected'''
    inputs = imas_io.profiles_from_imas(_simple_equilibrium_ids())
    assert np.all(inputs['profiles']['ffp_prof']['y'] == 0.0)
    assert np.allclose(inputs['profiles']['pp_prof']['y'], [1.0, 0.5, 0.0])
    assert inputs['F0'] == 6.0
    with pytest.raises(ValueError, match="Ip > 0"):
        imas_io.profiles_from_imas(_simple_equilibrium_ids(ip=-1.E6))
    ps = imas.IDSFactory().pulse_schedule()
    ps.ids_properties.homogeneous_time = 1
    ps.time = [0.0]
    ps.flux_control.ip.reference = [-1.E6]
    with pytest.raises(ValueError, match="Ip > 0"):
        imas_io.targets_from_pulse_schedule(ps)


def test_imas_loop_current():
    '''Passive loop current from element currents accounts for turns_with_sign'''
    cur = imas_io._loop_current([[10.0], [10.0]], [1.0, 1.0], 'series')
    assert np.allclose(cur, [10.0])
    cur = imas_io._loop_current([[10.0], [-10.0]], [1.0, -1.0], 'opposed')
    assert np.allclose(cur, [10.0])
    cur = imas_io._loop_current([[20.0]], [2.0], 'two_turn')
    assert np.allclose(cur, [10.0])
    with pytest.warns(UserWarning, match='different currents per turn'):
        imas_io._loop_current([[10.0], [5.0]], [1.0, 1.0], 'inconsistent')


def test_imas_save_mixed_versions():
    '''IDSs of different Data Dictionary versions cannot be written to one entry'''
    eq = imas.IDSFactory('4.0.0').equilibrium()
    pf = imas.IDSFactory('4.1.1').pf_active()
    with pytest.raises(ValueError, match="same Data Dictionary version"):
        imas_io.save_imas('mixed.nc', [eq, pf])


def test_imas_profile_coordinates():
    '''Without psi_norm, the profile coordinate is recovered relative to the LCFS (boundary.psi_norm); multi-species ions
    give the total density and the density-weighted temperature; npsi is bounded'''
    eq = _simple_equilibrium_ids(ffp=[1.0, 0.5, 0.0])
    ts = eq.time_slice[0]
    ts.profiles_1d.psi_norm = []
    ts.profiles_1d.psi = [10.0, 5.0, 0.0]
    ts.global_quantities.psi_axis = 10.0
    ts.global_quantities.psi_boundary = 0.0
    ts.boundary.psi_norm = 0.99
    assert np.allclose(imas_io.profiles_from_imas(eq)['profiles']['ffp_prof']['x'], [0.0, 0.495, 0.99])
    cp = imas.IDSFactory().core_profiles()
    cp.ids_properties.homogeneous_time = 1
    cp.time = [0.0]
    cp.profiles_1d.resize(1)
    p1 = cp.profiles_1d[0]
    p1.grid.rho_tor_norm = [0.0, 1.0]
    p1.t_i_average = [5.E3, 5.E3]
    p1.ion.resize(2)
    for ion, (dens, temp) in zip(p1.ion, ((4.E19, 1.E4), (6.E19, 2.E4))):
        ion.density = [dens, dens]
        ion.temperature = [temp, temp]
    kinetic = imas_io.profiles_from_imas(cp)
    assert np.allclose(kinetic['ni_prof']['y'], 1.E20)
    assert np.allclose(kinetic['ti_prof']['y'], 16.0)  # (4*10 + 6*20)/10 keV
    # Thermal populations only: a fast population does not shift the weighting, and a zero-density
    # point falls back to t_i_average; electrons given only as thermal density are used
    p1.grid.rho_tor_norm = [0.0, 0.5, 1.0]
    p1.t_i_average = [1.5E3, 1.5E3, 1.5E3]
    for ion, temp, fast in zip(p1.ion, (1.E3, 2.E3), (0.0, 3.E19)):
        ion.density = [5.E19+fast, 5.E19+fast, fast]
        ion.density_fast = [fast, fast, fast]
        ion.temperature = [temp, temp, temp]
    p1.electrons.density_thermal = [1.E20, 1.E20, 1.E19]
    kinetic = imas_io.profiles_from_imas(cp)
    assert np.allclose(kinetic['ti_prof']['y'], 1.5)
    assert np.allclose(kinetic['ni_prof']['y'], [1.E20, 1.E20, 0.0])
    assert np.allclose(kinetic['ne_prof']['y'], [1.E20, 1.E20, 1.E19])
    with pytest.raises(ValueError, match="npsi"):
        imas_io.equilibrium_to_imas(None, 0.0, npsi=3)


def test_imas_passive_resistance_turns():
    '''A loop given only by its resistance gets the resistivity of a region carrying N times the loop current'''
    factory = imas.IDSFactory()
    wall = factory.wall()
    wall.ids_properties.homogeneous_time = 2
    wall.description_2d.resize(1)
    wall.description_2d[0].limiter.unit.resize(1)
    wall.description_2d[0].limiter.unit[0].outline.r = [1.0, 2.0, 2.0, 1.0]
    wall.description_2d[0].limiter.unit[0].outline.z = [-1.0, -1.0, 1.0, 1.0]
    pf = factory.pf_active()
    pf.ids_properties.homogeneous_time = 2
    pp = factory.pf_passive()
    pp.ids_properties.homogeneous_time = 2
    pp.loop.resize(1)
    loop = pp.loop[0]
    loop.name = 'winding'
    loop.resistance = 1.E-3
    loop.element.resize(1)
    loop.element[0].turns_with_sign = 2.0
    geom = loop.element[0].geometry
    geom.geometry_type = 2
    geom.rectangle.r = 3.0
    geom.rectangle.z = 0.0
    geom.rectangle.width = 0.2
    geom.rectangle.height = 0.1
    domain = imas_io.domain_from_imas(wall, pf, {'plasma': 0.1, 'vacuum': 0.2, 'boundary': 0.5, 'coil': 0.1, 'conductor': 0.1},
                                      pf_passive=pp)
    eta = domain.region_info['WINDING']['eta']
    eta_exact = 1.E-3*0.1*np.log(3.1/2.9)/(2.0*np.pi*4.0)  # R = 2 pi eta / int(dA/R), region carries 2x the loop current
    assert abs(eta-eta_exact)/eta_exact < 1.E-12
    # The same cross-section as an outline with extra vertices on one side gives the same resistivity
    geom.geometry_type = 1
    geom.outline.r = [2.9, 2.95, 3.0, 3.05, 3.1, 3.1, 2.9]
    geom.outline.z = [-0.05, -0.05, -0.05, -0.05, -0.05, 0.05, 0.05]
    domain = imas_io.domain_from_imas(wall, pf, {'plasma': 0.1, 'vacuum': 0.2, 'boundary': 0.5, 'coil': 0.1, 'conductor': 0.1},
                                      pf_passive=pp)
    assert abs(domain.region_info['WINDING']['eta']-eta_exact)/eta_exact < 1.E-12


def run_psiscale_case(mp_q):
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
    from OpenFUSIONToolkit.TokaMaker.reconstruction import reconstruction
    from OpenFUSIONToolkit.TokaMaker.util import create_isoflux, create_power_flux_fun
    gs_mesh = gs_Domain()
    gs_mesh.define_region('plasma', 0.05, 'plasma')
    gs_mesh.add_polygon(create_isoflux(40, 1.0, 0.0, 0.4, 1.4, 0.2), 'plasma')
    mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()
    myOFT = OFT_env(nthreads=-1)
    mygs = TokaMaker(myOFT)
    mygs.settings.free_boundary = False
    mygs.setup_model(mesh_pts, mesh_lc, mesh_reg, order=2, F0=1.0)
    mygs.set_profiles(ffp_prof=create_power_flux_fun(40, 1.5, 2.0), pp_prof=create_power_flux_fun(40, 2.0, 1.5))
    mygs.set_targets(Ip=1.E5, pax=1.E3)
    mygs.init_psi(1.0, 0.0, 0.4, 1.4, 0.2)
    mygs.solve()
    results = {}
    imas_io.equilibrium_to_imas(mygs, 0.0, npsi=16, nr=16, nz=16)  # Exports before the fit
    # Fixed-boundary fit of a 20% higher plasma current with no current target fits the flux scale
    Ip0 = mygs.get_stats()['Ip']
    myrecon = reconstruction(mygs)
    myrecon.set_Ip(1.2*Ip0, 0.01*Ip0)
    mygs.set_targets(Ip=myOFT.float_disable_flag, pax=myOFT.float_disable_flag)
    myrecon.settings.fitP = False
    myrecon.settings.fitF = False
    myrecon.settings.pm = False
    myrecon.reconstruct()
    results['Ip_ratio'] = mygs.get_stats()['Ip']/Ip0
    # Rejected for both the solver's flux offset and one with the axis flux near zero
    rejected = []
    for offset in (0.0, mygs.psi_bounds[1]):
        mygs.set_psi(mygs.get_psi(False)-offset, update_bounds=True)
        try:
            imas_io.equilibrium_to_imas(mygs, 0.0, npsi=16, nr=16, nz=16)
            rejected.append(False)
        except ValueError as err:
            rejected.append('psiscale' in str(err))
    results['rejected'] = rejected
    results['axis_flux'] = float(mygs.psi_bounds[1])
    mp_q.put([results])
    oftpy_dump_cov()


@pytest.mark.coverage
def test_imas_fitted_flux_scale():
    '''An equilibrium whose flux scale was fitted in a fixed-boundary reconstruction is rejected, whatever its flux offset'''
    results = mp_run(run_psiscale_case, (), timeout=120)
    assert results is not None
    res = results[0]
    assert abs(res['Ip_ratio']-1.2) < 1.E-3
    assert res['rejected'] == [True, True]
    assert abs(res['axis_flux']) < 1.E-10


def test_imas_input_validation():
    '''Non-finite magnetic samples and incomplete equilibrium inputs fail at import, naming the field'''
    mag = imas.IDSFactory().magnetics()
    mag.ids_properties.homogeneous_time = 1
    mag.time = [0.0]
    mag.flux_loop.resize(1)
    mag.flux_loop[0].position.resize(1)
    mag.flux_loop[0].position[0].r = 8.0
    mag.flux_loop[0].position[0].z = 0.0
    mag.flux_loop[0].flux.data = [np.nan]
    with pytest.raises(ValueError, match="flux_loop"):
        imas_io.load_magnetics(None, mag, abs_err={'flux_loop': 1.E-3})
    eq = _simple_equilibrium_ids()
    eq.vacuum_toroidal_field.r0 = imas.ids_defs.EMPTY_FLOAT
    with pytest.raises(ValueError, match="r0"):
        imas_io.profiles_from_imas(eq)
    eq = _simple_equilibrium_ids()
    eq.time_slice[0].profiles_1d.dpressure_dpsi = [-2.0, -1.0]
    with pytest.raises(ValueError, match="dpressure_dpsi"):
        imas_io.profiles_from_imas(eq)

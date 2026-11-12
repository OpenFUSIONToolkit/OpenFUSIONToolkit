#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# TokaMaker tests: coupled TokaMaker-TORAX pulse design workflow (pulse_design.py) on ITER
#
# Seeds are handed to TORAX as get_fsa() flux surface averages, so the run exercises the 'tokamaker' geometry
# rather than re-contoured gEQDSKs.
#
# Skipped unless the torax package is installed. Uses the ITER test mesh (see tokamaker_testing.create_ITER_mesh),
# or the mesh file named by TOKAMAKER_ITER_MESH.
from __future__ import annotations

import json
import os
import tempfile
import time
from types import SimpleNamespace
from typing import Any, Dict, Optional

import numpy as np
import pytest
from tokamaker_testing import mp_run, create_ITER_mesh, ITER_isoflux_pts, ITER_x_point

_NUMERIC = (int, float, np.integer, np.floating)


def _round_outputs(obj: Any, ndigits: int = 2) -> Any:
    '''Round numbers (as floats) for JSON regression output, leaving bool, None, and str unchanged'''
    if obj is None or isinstance(obj, bool) or isinstance(obj, str):
        return obj
    if isinstance(obj, _NUMERIC):
        return round(float(obj), ndigits)
    if isinstance(obj, dict):
        return {k: _round_outputs(v, ndigits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        seq = [_round_outputs(x, ndigits) for x in obj]
        return type(obj)(seq) if isinstance(obj, tuple) else seq
    return obj


def _default_iter_mesh_path() -> str:
    '''Mesh named by TOKAMAKER_ITER_MESH, or the ITER test mesh in the working directory (built if needed)'''
    env = os.environ.get("TOKAMAKER_ITER_MESH")
    if env:
        return os.path.abspath(env)
    create_ITER_mesh()
    return os.path.abspath("ITER_mesh.h5")


def _array_to_profile_dict(profile_array: np.ndarray, psi_grid: np.ndarray) -> Dict[float, float]:
    return {float(p): float(v) for p, v in zip(psi_grid, profile_array)}


def _parabolic_profile(edge: float, core: float, psi_grid: np.ndarray, alpha: float = 1.8) -> np.ndarray:
    return edge + (core - edge) * (1.0 - np.asarray(psi_grid) ** alpha)


def _build_min_norm_coil_reg(mygs) -> None:
    reg_terms = []
    for name in mygs.coil_sets:
        if name.startswith("CS"):
            weight = 2.0e-2 if name.startswith("CS1") else 1.0e-2
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=weight))
        elif name.startswith("PF"):
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.0e-2))
        elif name.startswith("VS"):
            reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.0e-2))
    reg_terms.append(mygs.coil_reg_term({"#VSC": 1.0}, target=0.0, weight=1.0e2))
    mygs.set_coil_reg(reg_terms=reg_terms)


def _final_tmtx_equilibrium_stats(tt) -> Dict[str, Any]:
    '''Scalars from each TokaMaker equilibrium of the last fly() pass'''
    out: Dict[str, Any] = {}
    out["tmtx_last_completed_loop"] = int(getattr(tt, "_current_loop", -1))
    equil = tt.state.get("equil", {}) or {}
    tm_times = list(tt._tm_times)
    for i in sorted(equil.keys(), key=lambda k: int(k) if isinstance(k, (int, np.integer)) else k):
        eq = equil[i]
        if eq is None:
            continue
        try:
            stats = eq.get_stats(li_normalization="iter")
        except Exception:
            continue
        pfx = f"final_equil_i{i}_"
        out[f"{pfx}time_s"] = float(tm_times[i]) if i < len(tm_times) else float("nan")
        for key, val in stats.items():
            if key == "Ip_centroid":
                out[f"{pfx}Ip_centroid_R_m"] = float(val[0])
                out[f"{pfx}Ip_centroid_Z_m"] = float(val[1])
            elif isinstance(val, (float, int, np.floating, np.integer)):
                out[f"{pfx}{key}"] = float(val)
            elif isinstance(val, np.ndarray):
                out[f"{pfx}{key}"] = val.astype(float).tolist()
            else:
                out[f"{pfx}{key}"] = val
        try:
            out[f"{pfx}diverted"] = bool(eq.diverted)
            out[f"{pfx}psi_lcfs_Wb_per_rad"] = float(eq.psi_bounds[0])
            out[f"{pfx}psi_axis_Wb_per_rad"] = float(eq.psi_bounds[1])
        except Exception:
            pass
        try:
            out[f"{pfx}vloop_V"] = float(eq.calc_loopvoltage())
        except ValueError:
            out[f"{pfx}vloop_V"] = None
    return out


def _run_tokamaker_torax(
    mesh_path: Optional[str] = None,
    nthreads: int = 2,
    t_final: float = 5.0,
    tx_dt: float = 0.5,
    Ip_flattop: float = 15.0e6,
    pax_a: float = 6.2e5,
    pax_b: float = 5.7e5,
    n_surfaces: int = 200,
    last_surface_factor: float = 0.99,
    max_loop: int = 1,
    loop0: bool = True,
) -> Dict[str, Any]:
    '''Run the minimal ITER TokaMaker_TORAX case and return its rounded scalar outputs

    `mesh_path` defaults to TOKAMAKER_ITER_MESH or the ITER test mesh. Both seed equilibria and the TORAX
    schedule use the constant current `Ip_flattop` [A]; the seeds' axis pressures `pax_a` and `pax_b` [Pa] differ
    slightly so transport has something to evolve. `n_surfaces` and `last_surface_factor` set the get_fsa()
    surfaces shared by the seeds and the TokaMaker_TORAX object. `max_loop` is the highest counted coupling pass
    and `loop0` adds the initial pass (see TokaMaker_TORAX.fly).
    '''
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.pulse_design import (
        TokaMaker_TORAX,
        seed_from_equilibrium,
        summary,
    )
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun

    mesh = os.path.abspath(mesh_path or _default_iter_mesh_path())
    if not os.path.isfile(mesh):
        raise FileNotFoundError(
            f"ITER mesh not found at {mesh!r}."
        )

    r0_geo = 6.3
    b0 = 5.2
    f0 = r0_geo * b0

    # LCFS shape and lower X-point of the ITER baseline (L-mode inverse case of ITER_baseline_ex.ipynb)
    x_point = ITER_x_point
    diverted_pts = np.vstack((ITER_isoflux_pts, x_point))

    work_dir = tempfile.mkdtemp(prefix="test_tmtx_")

    myoft = OFT_env(nthreads=nthreads)
    mygs = TokaMaker(myoft)

    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(mesh)
    mygs.settings.maxits = 500
    mygs.setup_model(mesh_pts, mesh_lc, reg=mesh_reg, cond_dict=cond_dict, coil_dict=coil_dict, order=2, F0=f0)
    mygs.set_coil_vsc({"VS": 1.0})

    coil_bounds = {key: [-50.0e6, 50.0e6] for key in mygs.coil_sets}
    mygs.set_coil_bounds(coil_bounds)

    ffp_prof = create_power_flux_fun(40, 1.5, 2.0)
    pp_prof = create_power_flux_fun(40, 4.0, 1.0)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    seed_metrics: Dict[str, float] = {}
    # Seeds are kept as get_fsa() flux surface averages rather than written out as
    # gEQDSKs, so TORAX takes the 'tokamaker' geometry and nothing round trips through
    # a rectangular grid. seeds=[...] of gEQDSK filenames selects the 'eqdsk' geometry
    # instead; the two cannot be mixed in one config.
    seeds = []

    for idx, pax in enumerate((pax_a, pax_b)):
        mygs.set_isoflux_constraints(diverted_pts)
        mygs.set_saddle_constraints(x_point)
        mygs.set_targets(Ip=Ip_flattop, pax=pax)
        _build_min_norm_coil_reg(mygs)
        mygs.init_psi()
        mygs.solve()
        seed = seed_from_equilibrium(
            mygs, n_surfaces=n_surfaces, last_surface_factor=last_surface_factor
        )
        seeds.append(seed)
        seed_metrics[f"seed{idx}_pax_Pa"] = float(pax)
        seed_metrics[f"seed{idx}_psi_lcfs_Wb_per_rad"] = float(seed["psi_lcfs"])

    prev_cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        tm_times = np.linspace(0.0, t_final, 2)
        tt = TokaMaker_TORAX(
            t_init=0.0,
            t_final=t_final,
            eqtimes=[0.0, t_final],
            seeds=seeds,
            tokamaker_obj=mygs,
            tx_dt=tx_dt,
            tm_times=tm_times,
            last_surface_factor=last_surface_factor,
            n_surfaces=n_surfaces,
            truncate_eq=False,
        )

        tt.set_TORAX_grid(grid_type="n_rho", grid=51)

        n_sample = 100
        psi_sample = np.linspace(0.0, 1.0, n_sample)
        ne_init = _parabolic_profile(0.25e20, 2.0e20, psi_sample)
        te_init = _parabolic_profile(0.10, 10.0, psi_sample)
        ne = {0.0: _array_to_profile_dict(ne_init, psi_sample)}
        te = {0.0: _array_to_profile_dict(te_init, psi_sample)}
        tt.set_ne(ne, right_bc=0.25e20)
        tt.set_Te(te, right_bc=0.1)
        tt.set_Ti(te, right_bc=0.1)
        # H-mode pedestal-top targets; the TORAX L/H formation model owns the transition timing
        # (pedestal applied once P_SOL > P_LH), ramped over transition_time.
        tt.set_pedestal(
            config_mode="INTERNAL_BOUNDARY_CONDITION",
            ped_height_Te=3.0,
            ped_height_Ti=3.0,
            ped_height_ne=0.8e20,
            ped_top=0.9,
            transition_time=0.5,
            formation_model="martin_scaling",
        )

        heat_times = {0.0: 30.0e6}
        tt.set_heating(
            generic_heat=heat_times,
            generic_heat_loc=0.25,
            nbi_current=True,
            ecrh={0.0: 20.0e6},
            ecrh_loc=0.35,
            fusion=True,
            ei_exchange=True,
        )
        tt.set_fueling(gas_puff_S_total=1e22, gas_puff_decay_length=0.05, pellet_deposition_location=0.8, pellet_width=0.1, pellet_S_total={0: 5e21})

        tt.set_Ip({0.0: Ip_flattop})
        tt.set_plasma_composition(Zeff=1.6, main_ion={"D": 0.5, "T": 0.5}, impurity="Ne")
        tt.set_evolve(density=True, Ti=True, Te=True, current=True)

        tt.fly(
            run_name="tmp",
            max_loop=max_loop,
            loop0=loop0,
            output_mode=False,
            initial_relax=True,
        )

        phys = summary(tt)
    finally:
        os.chdir(prev_cwd)

    st = tt.state
    times = np.asarray(tt._tm_times)
    out: Dict[str, Any] = {}
    out.update(seed_metrics)
    # Resolved through the same helper fly() uses, so this asserts on the geometry
    # type TORAX was actually handed, not on how the seeds were built.
    out["torax_geometry_type"] = str(
        tt._tx_geometry({0.0: seeds[0]["geometry"]})["model_name"]
    )
    out["Ip_flattop_MA"] = float(Ip_flattop / 1.0e6)
    out["t_final_s"] = float(t_final)
    out["tx_dt_s"] = float(tx_dt)
    out["n_tm_times"] = int(len(times))

    # TokaMaker-side traces, which evolve with the TORAX coupling
    out["psi_lcfs_tm_first_Wb_per_rad"] = float(st["psi_lcfs_tm"][0])
    out["psi_lcfs_tm_last_Wb_per_rad"] = float(st["psi_lcfs_tm"][-1])
    out["Ip_tm_first_MA"] = float(st["Ip_tm"][0] / 1.0e6)
    out["Ip_tm_last_MA"] = float(st["Ip_tm"][-1] / 1.0e6)
    out["q95_tm_min"] = float(np.min(st["q95_tm"][st["q95_tm"] > 0])) if np.any(st["q95_tm"] > 0) else None
    out["q95_tm_last"] = float(st["q95_tm"][-1])
    out["beta_N_tm_max"] = float(np.max(st["beta_N_tm"]))

    out.update(_final_tmtx_equilibrium_stats(tt))

    for k, v in phys.items():
        if v is not None and k not in out:
            out[k] = v

    out["_work_dir"] = work_dir
    rounded: Dict[str, Any] = {}
    for k, v in out.items():
        if k.startswith("_"):
            rounded[k] = v
        else:
            rounded[k] = _round_outputs(v, ndigits=2)
    return rounded


def run_tokamaker_torax_case(mp_q):
    '''Run the TokaMaker-TORAX case in this (child) process, putting its outputs, or None on failure, on `mp_q`'''
    import traceback
    try:
        result = _run_tokamaker_torax()
    except Exception:
        traceback.print_exc()
        result = None
    mp_q.put(result)


@pytest.mark.slow
def test_tokamaker_torax() -> None:
    '''Coupled TokaMaker-TORAX pulse simulation at constant plasma current on ITER from two seed equilibria runs to
    completion and returns its summary (smoke test, skipped without TORAX)'''
    pytest.importorskip("torax", reason="TokaMaker_TORAX requires the torax package")
    result = mp_run(run_tokamaker_torax_case, (), timeout=1800)
    assert isinstance(result, dict)
    assert "Ip_flattop_MA" in result
    assert "t_final_s" in result
    # Guards the FSA geometry path: seeds built by seed_from_equilibrium must reach
    # TORAX as the 'tokamaker' geometry, not as re-contoured gEQDSKs.
    assert result["torax_geometry_type"] == "tokamaker"


def test_save_state_json(tmp_path) -> None:
    '''save_state writes loadable JSON and skips the live TokaMaker equilibrium objects'''
    pytest.importorskip("torax", reason="TokaMaker_TORAX requires the torax package")
    from OpenFUSIONToolkit.TokaMaker.pulse_design import TokaMaker_TORAX

    # Bypass __init__: only _state is needed, so no TokaMaker or TORAX run is required.
    tt = TokaMaker_TORAX.__new__(TokaMaker_TORAX)
    tt._state = {
        'Ip': np.array([1.0e6, 2.0e6]),
        'confinement_mode': np.array([False, True]),
        'psi_lcfs_tx': {0: np.float32(0.5), 1: np.int64(2)},
        'T_e': {0: {'x': np.linspace(0.0, 1.0, 3), 'y': np.ones(3), 'type': 'linterp'}},
        'equil': {0: object(), 1: object()},  # stand-ins for TokaMaker_equilibrium objects
    }
    fname = os.path.join(str(tmp_path), "state.json")
    tt.save_state(fname)
    with open(fname) as fid:
        loaded = json.load(fid)

    assert 'equil' not in loaded
    assert loaded['Ip'] == [1.0e6, 2.0e6]
    assert loaded['confinement_mode'] == [False, True]
    assert loaded['psi_lcfs_tx'] == {'0': 0.5, '1': 2}
    assert loaded['T_e']['0']['x'] == [0.0, 0.5, 1.0]
    assert 'equil' in tt._state  # the in-memory state is left untouched


class _FakeTxScalar:
    '''Minimal stand-in for a TORAX output scalar (xarray DataArray with a time coordinate)'''

    def __init__(self, t, y):
        self.coords = {'time': SimpleNamespace(values=t)}
        self._y = y

    def to_numpy(self):
        return self._y


@pytest.mark.parametrize("pedestal_cfg, used_label", [
    ({'structure_model': {'model_name': 'set_T_top_n_top'}, 'transition_time_width': 0.5,
      'formation_model': {'model_name': 'martin_scaling'}}, r'$P_{LH}$ Martin (used)'),
    ({'structure_model': {'model_name': 'set_T_top_n_top'}, 'transition_time_width': 0.5,
      'formation_model': {'model_name': 'delabie_scaling'}}, r'$P_{LH}$ Delabie (used)'),
    # No model named: TORAX defaults the formation model to 'prescribed', which has no P_LH.
    ({'structure_model': {'model_name': 'set_T_top_n_top'}, 'formation_model': {}}, None),
    ({}, None),
])
def test_plh_used_label(monkeypatch, pedestal_cfg, used_label) -> None:
    '''The P_LH "(used)" label follows the formation model in the validated TORAX config'''
    torax = pytest.importorskip("torax", reason="TokaMaker_TORAX requires the torax package")
    import matplotlib.pyplot as plt
    from OpenFUSIONToolkit.TokaMaker import pulse_design

    tx_config = torax.ToraxConfig.from_dict({
        'profile_conditions': {'I_p': 10.0e6, 'psi_source': 'j'}, 'plasma_composition': {}, 'numerics': {},
        'geometry': {'model_name': 'circular'}, 'pedestal': pedestal_cfg,
        'sources': {}, 'solver': {}, 'transport': {},
    })

    t = np.linspace(0.0, 10.0, 11)
    scalars = {name: _FakeTxScalar(t, np.full(t.shape, val)) for name, val in [
        ('P_LH_martin', 30.0e6), ('P_LH_delabie', 20.0e6), ('P_SOL_total', 25.0e6),
        ('n_e_line_avg', 5.0e19), ('n_e_min_P_LH', 3.0e19), ('I_p', 10.0e6),
        ('a_minor', 2.0), ('R_major', 6.2), ('B_0', 5.3), ('f_GW_n_e_line_avg', 0.5),
    ]}
    tt = SimpleNamespace(
        _data_tree=SimpleNamespace(scalars=SimpleNamespace(**scalars), profiles=SimpleNamespace()),
        _tx_formation_model=tx_config.pedestal.formation_model.model_name,
        _lh_time=None, _t_init=0.0, _t_final=10.0, _current_loop=1,
    )
    figs = []
    monkeypatch.setattr(pulse_design, '_save_or_display', lambda fig, save_path=None, display=True: figs.append(fig))
    plt.switch_backend('Agg')
    pulse_design.plot_PLH_components(tt, display=False)

    # Check every P_SOL vs P_LH panel, found by its P_SOL curve rather than its grid position.
    gate_axes = []
    for ax in figs[0].axes:
        if r'$P_{SOL}$' in [line.get_label() for line in ax.get_lines()]:
            gate_axes.append(ax)
    assert gate_axes
    expected = [used_label] if used_label is not None else []
    for ax in gate_axes:
        used = [line.get_label() for line in ax.get_lines() if line.get_label().endswith('(used)')]
        assert used == expected
    plt.close(figs[0])


def main() -> None:
    t_wall0 = time.perf_counter()
    result = _run_tokamaker_torax()
    elapsed_s = time.perf_counter() - t_wall0
    slim = {k: v for k, v in result.items() if not k.startswith("_")}
    print(json.dumps(slim, indent=2, sort_keys=True))
    print(f"\nTotal wall time (script): {elapsed_s:.2f} s", flush=True)


if __name__ == "__main__":
    main()

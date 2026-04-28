"""
Minimal ITER test for pulse_planning.py: the TokaMaker + TORAX coupled pulse simulation (TokTox).

Builds two flattop-like seed equilibria (same Ip, slightly different axis pressure),
runs a 5 s simulation at constant Ip, and returns a flat dict of scalars.

Mesh file
---------
Expects ``ITER_mesh.h5`` produced from ``src/examples/TokaMaker/ITER/ITER_mesh_ex.ipynb``
(next to ``ITER_geom.json``), or set environment variable ``TOKAMAKER_ITER_MESH``
to the ``.h5`` path.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from typing import Any, Dict, Optional

import numpy as np
import pytest

_NUMERIC = (int, float, np.integer, np.floating)


def _round_outputs(obj: Any, ndigits: int = 2) -> Any:
    """Round floats (and ints as floats) for JSON-friendly regression output; preserve bool, None, str."""
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

# Repo src/python on path when run as script from elsewhere
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_TESTS_DIR, "..", ".."))
_PYTHON_SRC = os.path.join(_REPO_ROOT, "src", "python")
if _PYTHON_SRC not in sys.path:
    sys.path.insert(0, _PYTHON_SRC)


def _default_iter_mesh_path() -> str:
    env = os.environ.get("TOKAMAKER_ITER_MESH")
    if env:
        return os.path.abspath(env)
    return os.path.join(
        _REPO_ROOT,
        "src",
        "examples",
        "TokaMaker",
        "ITER",
        "ITER_mesh.h5",
    )


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


def _final_toktox_equilibrium_stats(tt) -> Dict[str, Any]:
    """Scalars from each TokaMaker equilibrium produced in the last ``fly()`` TM pass."""
    out: Dict[str, Any] = {}
    out["toktox_last_completed_loop"] = int(getattr(tt, "_current_loop", -1))
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


def _iter_baseline_shape() -> tuple[np.ndarray, np.ndarray]:
    """LCFS isoflux points and X-point from ``ITER_baseline_ex.ipynb`` (L-mode inverse case)."""
    isoflux_pts = np.array(
        [
            [8.20, 0.41],
            [8.06, 1.46],
            [7.51, 2.62],
            [6.14, 3.78],
            [4.51, 3.02],
            [4.26, 1.33],
            [4.28, 0.08],
            [4.49, -1.34],
            [7.28, -1.89],
            [8.00, -0.68],
        ]
    )
    x_point = np.array([[5.125, -3.4]])
    return isoflux_pts, x_point


def test_toktox(
    mesh_path: Optional[str] = None,
    nthreads: int = 2,
    t_final: float = 5.0,
    tx_dt: float = 0.5,
    Ip_flattop: float = 15.0e6,
    pax_a: float = 6.2e5,
    pax_b: float = 5.7e5,
    eqdsk_nr: int = 100,
    eqdsk_nz: int = 100,
    max_loop: int = 1,
) -> Dict[str, Any]:
    """
    Run the minimal ITER TokTox benchmark and return numerical outputs.

    Parameters
    ----------
    mesh_path
        Path to ``ITER_mesh.h5``. Default: ``TOKAMAKER_ITER_MESH`` or
        ``src/examples/TokaMaker/ITER/ITER_mesh.h5`` under the repo root.
    ip_flat
        Constant plasma current [A] for both seeds and the TORAX Ip schedule.
    pax_a, pax_b
        Axis pressure targets [Pa] for the two seeds (slight mismatch so transport
        has something to do while Ip stays flat).
    """
    try:
        import torax  # noqa: F401
    except ImportError:
        pytest.skip("TokTox requires the ``torax`` package.")

    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.TokaMaker import TokaMaker
    from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
    from OpenFUSIONToolkit.TokaMaker.pulse_planning import TokTox, summary
    from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun, read_eqdsk

    mesh = os.path.abspath(mesh_path or _default_iter_mesh_path())
    if not os.path.isfile(mesh):
        raise FileNotFoundError(
            f"ITER mesh not found at {mesh!r}."
        )

    r0_geo = 6.3
    b0 = 5.2
    z0 = 0.5
    f0 = r0_geo * b0

    isoflux_pts, x_point = _iter_baseline_shape()
    diverted_pts = np.vstack((isoflux_pts, x_point))

    work_dir = tempfile.mkdtemp(prefix="test_toktox_")
    eq_a = os.path.join(work_dir, "seed_t0.eqdsk")
    eq_b = os.path.join(work_dir, "seed_t5.eqdsk")

    myoft = OFT_env(nthreads=nthreads)
    mygs = TokaMaker(myoft)

    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh(mesh)
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.settings.maxits = 500
    mygs.setup(order=2, F0=f0)
    mygs.set_coil_vsc({"VS": 1.0})

    coil_bounds = {key: [-50.0e6, 50.0e6] for key in mygs.coil_sets}
    mygs.set_coil_bounds(coil_bounds)

    ffp_prof = create_power_flux_fun(40, 1.5, 2.0)
    pp_prof = create_power_flux_fun(40, 4.0, 1.0)
    mygs.set_profiles(ffp_prof=ffp_prof, pp_prof=pp_prof)

    seed_metrics: Dict[str, float] = {}

    for idx, (pax, eq_path) in enumerate(((pax_a, eq_a), (pax_b, eq_b))):
        mygs.set_isoflux_constraints(diverted_pts)
        mygs.set_saddle_constraints(x_point)
        mygs.set_targets(Ip=Ip_flattop, pax=pax)
        _build_min_norm_coil_reg(mygs)
        mygs.init_psi()
        mygs.solve()
        seed_metrics[f"seed{idx}_pax_Pa"] = float(pax)
        mygs.save_eqdsk(eq_path, cocos=2, nr=eqdsk_nr, nz=eqdsk_nz)
        g = read_eqdsk(eq_path)
        seed_metrics[f"seed{idx}_psi_lcfs_Wb_per_rad"] = float(-g["psibry"])

    prev_cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        tm_times = np.linspace(0.0, t_final, 2)
        tt = TokTox(
            t_init=0.0,
            t_final=t_final,
            eqtimes=[0.0, t_final],
            g_eqdsk_arr=["seed_t0.eqdsk", "seed_t5.eqdsk"],
            tx_dt=tx_dt,
            tm_times=tm_times,
            last_surface_factor=0.99,
            cocos=2,
            oft_env=myoft,
            oft_threads=nthreads,
            truncate_eq=False,
        )

        tt.set_tx_grid(grid_type="n_rho", grid=51)

        n_sample = 100
        psi_sample = np.linspace(0.0, 1.0, n_sample)
        ne_init = _parabolic_profile(0.25e20, 2.0e20, psi_sample)
        te_init = _parabolic_profile(0.10, 10.0, psi_sample)
        ne = {0.0: _array_to_profile_dict(ne_init, psi_sample)}
        te = {0.0: _array_to_profile_dict(te_init, psi_sample)}
        tt.set_ne(ne)
        tt.set_Te(te)
        tt.set_Ti(te)

        tt.set_right_bc(ne_right_bc=0.25e20, Te_right_bc=0.1, Ti_right_bc=0.1)
        tt.set_pedestal(set_pedestal=True, T_i_ped=3.0, T_e_ped=3.0, n_e_ped=0.8e20)

        heat_times = {0.0: 30.0e6}
        tt.set_heating(
            generic_heat=heat_times,
            generic_heat_loc=0.25,
            nbi_current=True,
            ecrh={0.0: 20.0e6},
            ecrh_loc=0.35,
        )
        tt.set_sources(fusion=True, ei_exchange=True)

        tt.set_gas_puff(S_total=1e22, decay_length=0.05) # estimations based on budney2008
        tt.set_pellet(pellet_deposition_location=0.8, pellet_width=0.1, S_total={0: 5e21})

        tt._tm = mygs
        tt.set_coil_reg(coil_bounds=coil_bounds, updownsym=False)

        tt.set_Ip({0.0: Ip_flattop})
        tt.set_plasma_composition(main_ion={"D": 0.5, "T": 0.5}, impurity="Ne")
        tt.set_Zeff(1.6)
        tt.set_evolve(density=True, Ti=True, Te=True, current=True)

        tt.fly(run_name="tmp", max_loop=max_loop, output_mode=False, initial_relax=True)

        phys = summary(tt)
    finally:
        os.chdir(prev_cwd)

    st = tt.state
    times = np.asarray(tt._tm_times)
    out: Dict[str, Any] = {}
    out.update(seed_metrics)
    out["Ip_flattop_MA"] = float(Ip_flattop / 1.0e6)
    out["t_final_s"] = float(t_final)
    out["tx_dt_s"] = float(tx_dt)
    out["n_tm_times"] = int(len(times))

    # TokaMaker-side traces (good for regression; values evolve with TORAX coupling).
    out["psi_lcfs_tm_first_Wb_per_rad"] = float(st["psi_lcfs_tm"][0])
    out["psi_lcfs_tm_last_Wb_per_rad"] = float(st["psi_lcfs_tm"][-1])
    out["Ip_tm_first_MA"] = float(st["Ip_tm"][0] / 1.0e6)
    out["Ip_tm_last_MA"] = float(st["Ip_tm"][-1] / 1.0e6)
    out["q95_tm_min"] = float(np.min(st["q95_tm"][st["q95_tm"] > 0])) if np.any(st["q95_tm"] > 0) else None
    out["q95_tm_last"] = float(st["q95_tm"][-1])
    out["beta_N_tm_max"] = float(np.max(st["beta_N_tm"]))

    out.update(_final_toktox_equilibrium_stats(tt))

    # Merge physics summary (TORAX / integrated quantities).
    for k, v in phys.items():
        if v is not None and k not in out:
            out[k] = v

    out["_work_dir"] = work_dir

    # tt.plot_scalars(display=True)
    # tt.plot_profile_evolution(display=True, one_plot=True)
    # tt.plot_lcfs_evolution(display=True, one_plot=True)

    rounded: Dict[str, Any] = {}
    for k, v in out.items():
        if k.startswith("_"):
            rounded[k] = v
        else:
            rounded[k] = _round_outputs(v, ndigits=2)
    return rounded


def main() -> None:
    t_wall0 = time.perf_counter()
    result = test_toktox()
    elapsed_s = time.perf_counter() - t_wall0
    slim = {k: v for k, v in result.items() if not k.startswith("_")}
    print(json.dumps(slim, indent=2, sort_keys=True))
    print(f"\nTotal wall time (script): {elapsed_s:.2f} s", flush=True)


if __name__ == "__main__":
    main()

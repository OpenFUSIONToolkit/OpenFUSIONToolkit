"""
ITER H-mode bootstrap test script.

By default runs BOTH bootstrap methods and prints a side-by-side comparison.
Pass --mode to restrict to a single method.

  --mode both              (default) run external then internal and compare
  --mode external          Python-side solve_with_bootstrap() from bootstrap.py
  --mode internal          GS-internal jphi-split-bootstrap (Fortran-embedded)
  --mode internal-compare  Run internal solver twice: once with --isolate-edge
                           and once with --parameterize-edge, then plot the diff

Usage
-----
  python ITER_Hmode_bootstrap_test.py [--mode both|external|internal]
                                      [--isolate-edge]
                                      [--parameterize-edge]
                                      [--scale-jBS 1.0]
                                      [--taper-edge]
                                      [--nthreads 2]
"""

import argparse
import datetime
import os
import sys
import time

import tempfile

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
mu0 = 4.0 * np.pi * 1e-7

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--mode', choices=['both', 'external', 'internal', 'internal-compare'],
                    default='both',
                    help='Which bootstrap method(s) to run (default: both)')
parser.add_argument('--isolate-edge', action='store_true',
                    help='Isolate the edge bootstrap spike from the bulk')
parser.add_argument('--parameterize-edge', action='store_true',
                    help='Fit edge spike with parametrised skew-normal (implies --isolate-edge)')
parser.add_argument('--scale-jBS', type=float, default=1.0,
                    help='Scale factor applied to the bootstrap current (default 1.0)')
parser.add_argument('--diagnose-bs', action='store_true',
                    help='Print alpha/Ip scalars, j_BS stats, and full profile tables each NL iteration')
parser.add_argument('--taper-edge', action='store_true',
                    help='Smooth the edge bootstrap spike to zero at the plasma edge')
parser.add_argument('--taper-edge-psi0', type=float, default=0.99,
                    help='psi_N (0=axis, 1=plasma edge) where taper begins (default 0.99)')
parser.add_argument('--taper-edge-shape', type=int, default=2,
                    help='Taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power (default 2)')
parser.add_argument('--nthreads', type=int, default=2,
                    help='Number of OFT threads (default 2)')
args = parser.parse_args()

run_external         = args.mode in ('both', 'external')
run_internal         = args.mode in ('both', 'internal')
run_internal_compare = args.mode == 'internal-compare'

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
oft_path = os.environ.get('OFT_INSTALL_DIR', None)
if oft_path is not None:
    sys.path.append(os.path.join(oft_path, 'python'))

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_power_flux_fun
from OpenFUSIONToolkit.TokaMaker.bootstrap import solve_with_bootstrap, Hmode_profiles

# ---------------------------------------------------------------------------
# Shared physical parameters
# ---------------------------------------------------------------------------
Ip_target  = 13.0e6   # A
P0_target  = 6.2e5    # Pa
Zeff_val   = 1.5      # effective charge (scalar)
EC         = 1.602e-19
n_sample   = 257
n_sample_plot = n_sample
psi_sample = np.linspace(0.0, 1.0, n_sample)

xphalf = 0.965
ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n_sample,
                    expin=1.6, expout=1.6, widthp=0.35, xphalf=xphalf) * 1e20
Te = Hmode_profiles(edge=1500., ped=5000., core=21000., rgrid=n_sample,
                    expin=1.3, expout=1.7, widthp=0.1, xphalf=xphalf)
ni = ne.copy()
Ti = Te.copy()
pressure   = EC * ne * Te + EC * ni * Ti

inductive_jphi = create_power_flux_fun(n_sample, 2.25, 2.5)['y']

# ---------------------------------------------------------------------------
# Helper: build a fresh, initialised TokaMaker ready for bootstrap solves
# ---------------------------------------------------------------------------
def build_gs(mygs):
    mygs.reset()
    mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('ITER_mesh.h5')
    mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
    mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    mygs.settings.maxits = 100
    mygs.setup(order=2, F0=5.3 * 6.2)

    mygs.set_coil_vsc({'VS': 1.0})
    mygs.set_coil_bounds({key: [-50.e6, 50.e6] for key in mygs.coil_sets})

    isoflux_pts = np.array([
        [ 8.20,  0.41], [ 8.06,  1.46], [ 7.51,  2.62],
        [ 6.14,  3.78], [ 4.51,  3.02], [ 4.26,  1.33],
        [ 4.28,  0.08], [ 4.49, -1.34], [ 7.28, -1.89],
        [ 8.00, -0.68],
    ])
    x_point = np.array([[5.125, -3.4]])
    mygs.set_isoflux_constraints(np.vstack((isoflux_pts, x_point)))
    mygs.set_saddle_constraints(x_point)

    reg_terms = []
    for name in mygs.coil_sets:
        w = 2.e-2 if name.startswith('CS1') else 1.e-2
        reg_terms.append(mygs.coil_reg_term({name: 1.0}, target=0.0, weight=w))
    reg_terms.append(mygs.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=1.e2))
    mygs.set_coil_reg(reg_terms=reg_terms)

    mygs.set_targets(Ip=Ip_target, pax=P0_target)
    mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
    mygs.settings.pm = False
    mygs.update_settings()
    mygs.solve()
    print("  Initial (no-bootstrap) solve complete.")
    return mygs


# ---------------------------------------------------------------------------
# Helper: print a result block
# ---------------------------------------------------------------------------
def print_result(label, elapsed, stats, coil_currents):
    print(f"\n{'='*60}")
    print(f"  {label}  ({elapsed:.3f} s)")
    print(f"{'='*60}")
    print(f"  Ip            = {stats.get('Ip', float('nan'))/1e6:.4f} MA")
    print(f"  Pax           = {stats.get('P_ax', float('nan'))/1e3:.4f} kPa")
    print(f"  Beta_pol      = {stats.get('beta_pol', float('nan')):.4f} %")
    print(f"  Beta_tor      = {stats.get('beta_tor', float('nan')):.4f} %")
    print(f"  q_0 (axis)    = {stats.get('q0', float('nan')):.4f}")
    print(f"  q_95          = {stats.get('q95', float('nan')):.4f}")
    print("  Coil Currents [MA]:")
    for key, val in coil_currents.items():
        print(f"    {key:10} {val/1e6:10.3f}")


# ---------------------------------------------------------------------------
# Helper: compare two stat dicts
# ---------------------------------------------------------------------------
def print_comparison(stats_ext, stats_int):
    keys = [('Ip', 'Ip [MA]', 1e6), ('P_ax', 'Pax [kPa]', 1e3),
            ('beta_pol', 'beta_pol [%]', 1), ('beta_tor', 'beta_tor [%]', 1),
            ('q0', 'q0', 1), ('q95', 'q95', 1)]
    print(f"\n{'='*60}")
    print("  Comparison: external vs internal")
    print(f"{'='*60}")
    print(f"  {'Quantity':<18} {'External':>12} {'Internal':>12} {'Diff %':>10}")
    print(f"  {'-'*54}")
    for key, label, scale in keys:
        v_ext = stats_ext.get(key, None)
        v_int = stats_int.get(key, None)
        if v_ext is None or v_int is None:
            continue
        pct = 100.0 * (v_int - v_ext) / v_ext if v_ext != 0 else float('nan')
        print(f"  {label:<18} {v_ext/scale:>12.4f} {v_int/scale:>12.4f} {pct:>9.3f}%")


# ---------------------------------------------------------------------------
# Shared OFT environment and TokaMaker (only one instance each allowed per process)
# ---------------------------------------------------------------------------
myOFT = OFT_env(nthreads=args.nthreads)
mygs  = TokaMaker(myOFT)

# ---------------------------------------------------------------------------
# Run external bootstrap
# ---------------------------------------------------------------------------
results_ext = None
stats_ext   = None
elapsed_ext = None

if run_external:
    print("\n" + "="*60)
    print("  EXTERNAL bootstrap  (bootstrap.py / Sauter-Redl)")
    print("="*60)
    build_gs(mygs)

    # --parameterize-edge requires --isolate-edge for the external solver.
    ext_isolate_edge = args.isolate_edge or args.parameterize_edge
    if args.parameterize_edge and not args.isolate_edge:
        print('Note: --parameterize-edge implies --isolate-edge for the external solver.')

    t0 = time.perf_counter()
    results_ext = solve_with_bootstrap(
        mygs,
        ne, Te, ni, Ti,
        np.full(n_sample, Zeff_val),
        Ip_target,
        inductive_jphi,
        scale_jBS=args.scale_jBS,
        isolate_edge_jBS=ext_isolate_edge,
        parameterize_jBS=args.parameterize_edge,
        diagnostic_plots=False,
        verbose=True,
        diagnose_bs=args.diagnose_bs,
    )
    elapsed_ext = time.perf_counter() - t0

    stats_ext = mygs.get_stats()
    coils_ext, _ = mygs.get_coil_currents()
    print_result("EXTERNAL", elapsed_ext, stats_ext, coils_ext)

    # --- Save equilibrium plot ---
    fig_eq, ax_eq = plt.subplots(1, 1, figsize=(7, 7))
    mygs.plot_machine(fig_eq, ax_eq, coil_colormap='seismic', coil_symmap=True,
                      coil_scale=1e-6, coil_clabel=r'$I_C$ [MA]')
    mygs.plot_psi(fig_eq, ax_eq, xpoint_color=None, vacuum_nlevels=4)
    mygs.plot_constraints(fig_eq, ax_eq, isoflux_color='tab:red', isoflux_marker='o')
    ax_eq.set_title('External bootstrap equilibrium')
    fig_eq.savefig(os.path.join(script_dir, 'eq_external.png'), dpi=300, bbox_inches='tight')
    plt.close(fig_eq)
    print(f"  Saved: {os.path.join(script_dir, 'eq_external.png')}")

    # --- Extract 1D profiles ---
    psi_e, F_e, Fp_e, P_e, Pp_e = mygs.get_profiles(npsi=n_sample_plot, psi_pad=1e-3)
    _, q_e, ravgs_e, _, _, _     = mygs.get_q(npsi=n_sample_plot, psi_pad=1e-3)
    jtor_e = F_e * Fp_e * ravgs_e[1] / mu0 + Pp_e * ravgs_e[0]

# ---------------------------------------------------------------------------
# Run internal bootstrap
# ---------------------------------------------------------------------------
results_int = None
stats_int   = None
elapsed_int = None

if run_internal:
    print("\n" + "="*60)
    print("  INTERNAL bootstrap  (jphi-split-bootstrap / Fortran-embedded)")
    print("="*60)
    build_gs(mygs)

    t0 = time.perf_counter()

    ne_i, Te_i, ni_i, Ti_i = ne, Te, ni, Ti
    pressure_i = pressure
    pp_vals = np.gradient(pressure_i) / np.gradient(psi_sample)
    # Enforce P' edge condition
    pp_vals[-1] = 0.0

    mygs.set_kinetic_profiles(
        te_prof={'type': 'linterp', 'x': psi_sample, 'y': Te_i / 1e3},
        ne_prof={'type': 'linterp', 'x': psi_sample, 'y': ne_i},
        ti_prof={'type': 'linterp', 'x': psi_sample, 'y': Ti_i / 1e3},
        ni_prof={'type': 'linterp', 'x': psi_sample, 'y': ni_i},
        Zeff=Zeff_val,
    )
    mygs.set_boot_ops(
        isolate_edge_jBS=args.isolate_edge,
        parameterize_jBS=args.parameterize_edge,
        scale_jBS=args.scale_jBS,
        Zeff=Zeff_val,
        diagnose_bs=args.diagnose_bs,
        taper_edge_jBS=args.taper_edge,
        taper_edge_psi0=args.taper_edge_psi0,
        taper_edge_shape=args.taper_edge_shape,
    )
    mygs.set_profiles(
        ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': inductive_jphi},
        pp_prof={'type': 'linterp', 'x': psi_sample, 'y': pp_vals / pp_vals[0]},
        foffset=5.3 * 6.2,
    )
    mygs.set_targets(Ip=Ip_target, pax=float(pressure_i[0]))
    mygs.settings.pm = True
    mygs.update_settings()
    mygs.solve()

    elapsed_int = time.perf_counter() - t0

    stats_int = mygs.get_stats()
    coils_int, _ = mygs.get_coil_currents()
    print_result("INTERNAL", elapsed_int, stats_int, coils_int)

    # --- Save equilibrium plot ---
    fig_eq, ax_eq = plt.subplots(1, 1, figsize=(7, 7))
    mygs.plot_machine(fig_eq, ax_eq, coil_colormap='seismic', coil_symmap=True,
                      coil_scale=1e-6, coil_clabel=r'$I_C$ [MA]')
    mygs.plot_psi(fig_eq, ax_eq, xpoint_color=None, vacuum_nlevels=4)
    mygs.plot_constraints(fig_eq, ax_eq, isoflux_color='tab:red', isoflux_marker='o')
    ax_eq.set_title('Internal bootstrap equilibrium')
    fig_eq.savefig(os.path.join(script_dir, 'eq_internal.png'), dpi=300, bbox_inches='tight')
    plt.close(fig_eq)
    print(f"  Saved: {os.path.join(script_dir, 'eq_internal.png')}")

    # --- Extract 1D profiles ---
    psi_i, F_i, Fp_i, P_i, Pp_i = mygs.get_profiles(npsi=n_sample_plot, psi_pad=1e-3)
    _, q_i, ravgs_i, _, _, _     = mygs.get_q(npsi=n_sample_plot, psi_pad=1e-3)
    jtor_i = F_i * Fp_i * ravgs_i[1] / mu0 + Pp_i * ravgs_i[0]


# ---------------------------------------------------------------------------
# Comparison (only when both ran)
# ---------------------------------------------------------------------------
if run_external and run_internal:
    print_comparison(stats_ext, stats_int)
    print(f"\n  Timing:  external {elapsed_ext:.3f} s  |  internal {elapsed_int:.3f} s"
          f"  (ratio {elapsed_ext/elapsed_int:.2f}x)")

    # --- Profile comparison figure (3 rows × 3 cols) ---
    def _pct_diff(a, b):
        denom = np.where(np.abs(a) > 1e-10 * np.nanmax(np.abs(a)), np.abs(a), np.nan)
        return 100.0 * (b - a) / denom

    option_parts = []
    if args.isolate_edge:     option_parts.append('isolate-edge')
    if args.parameterize_edge: option_parts.append('parameterize-edge')
    if args.taper_edge:       option_parts.append(f'taper(psi0={args.taper_edge_psi0},sh={args.taper_edge_shape})')
    if args.scale_jBS != 1.0: option_parts.append(f'scale={args.scale_jBS}')
    option_tag   = ('_' + '_'.join(option_parts)) if option_parts else '_plain'
    option_label = ' [' + ', '.join(option_parts) + ']' if option_parts else ' [plain]'

    fig, axs = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle('Bootstrap comparison: external vs internal' + option_label, fontsize=13)
    psi_lbl = r'$\hat{\psi}$'
    kw = dict(lw=1.5)

    # Row 0: F, F·F', P
    axs[0, 0].plot(psi_e, F_e,          label='External', **kw)
    axs[0, 0].plot(psi_i, F_i,          label='Internal', ls='--', **kw)
    axs[0, 0].set_ylabel('F [T·m]');  axs[0, 0].legend(fontsize=8)

    axs[0, 1].plot(psi_e, F_e * Fp_e, **kw)
    axs[0, 1].plot(psi_i, F_i * Fp_i, ls='--', **kw)
    axs[0, 1].set_ylabel("F·F' [T²]")

    axs[0, 2].plot(psi_e, P_e / 1e6, **kw)
    axs[0, 2].plot(psi_i, P_i / 1e6, ls='--', **kw)
    axs[0, 2].set_ylabel('P [MPa]')

    # Row 1: P', q, jφ
    axs[1, 0].plot(psi_e, Pp_e / 1e3, **kw)
    axs[1, 0].plot(psi_i, Pp_i / 1e3, ls='--', **kw)
    axs[1, 0].set_ylabel("P' [kPa]")

    axs[1, 1].plot(psi_e, q_e, **kw)
    axs[1, 1].plot(psi_i, q_i, ls='--', **kw)
    axs[1, 1].axhline(1.0, ls=':', c='k', alpha=0.5)
    axs[1, 1].set_ylabel('q')

    axs[1, 2].plot(psi_e, jtor_e / 1e6, **kw)
    axs[1, 2].plot(psi_i, jtor_i / 1e6, ls='--', **kw)
    axs[1, 2].set_ylabel(r'$\langle j_\phi \rangle$ [MA/m²]')

    # Row 2: percentage differences
    for col, (ye, yi, lbl) in enumerate([
            (F_e * Fp_e, F_i * Fp_i, "Δ F·F' [%]"),
            (q_e,        q_i,        'Δ q [%]'),
            (jtor_e,     jtor_i,     r'Δ $j_\phi$ [%]'),
    ]):
        axs[2, col].plot(psi_e, _pct_diff(ye, yi), c='tab:red', **kw)
        axs[2, col].axhline(0, ls=':', c='k', alpha=0.5)
        axs[2, col].set_ylabel(lbl)

    for ax in axs.flat:
        ax.grid(ls=':')
        ax.set_xlabel(psi_lbl)

    fig.tight_layout()
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    fig.text(0.99, 0.01, timestamp, ha='right', va='bottom',
             fontsize=7, color='gray', transform=fig.transFigure)
    out_path = os.path.join(script_dir, f'profiles_comparison{option_tag}_{timestamp}.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved: {out_path}")

# ---------------------------------------------------------------------------
# Single-run profile plot (external-only or internal-only)
# ---------------------------------------------------------------------------
if run_external != run_internal and not run_internal_compare:
    if run_external:
        psi_s, F_s, Fp_s, P_s, Pp_s = psi_e, F_e, Fp_e, P_e, Pp_e
        q_s, jtor_s = q_e, jtor_e
        label_s = 'External'
        tag_s   = 'external'
    else:
        psi_s, F_s, Fp_s, P_s, Pp_s = psi_i, F_i, Fp_i, P_i, Pp_i
        q_s, jtor_s = q_i, jtor_i
        label_s = 'Internal'
        tag_s   = 'internal'

    option_parts = []
    if args.isolate_edge:      option_parts.append('isolate-edge')
    if args.parameterize_edge: option_parts.append('parameterize-edge')
    if args.taper_edge:        option_parts.append(f'taper(psi0={args.taper_edge_psi0},sh={args.taper_edge_shape})')
    if args.scale_jBS != 1.0:  option_parts.append(f'scale={args.scale_jBS}')
    option_tag   = ('_' + '_'.join(option_parts)) if option_parts else '_plain'
    option_label = ' [' + ', '.join(option_parts) + ']' if option_parts else ' [plain]'

    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(f'{label_s} bootstrap profiles' + option_label, fontsize=13)
    psi_lbl = r'$\hat{\psi}$'
    kw = dict(lw=1.5)

    axs[0, 0].plot(psi_s, F_s, **kw);          axs[0, 0].set_ylabel('F [T·m]')
    axs[0, 1].plot(psi_s, F_s * Fp_s, **kw);   axs[0, 1].set_ylabel("F·F' [T²]")
    axs[0, 2].plot(psi_s, P_s / 1e6, **kw);    axs[0, 2].set_ylabel('P [MPa]')
    axs[1, 0].plot(psi_s, Pp_s / 1e3, **kw);   axs[1, 0].set_ylabel("P' [kPa]")
    axs[1, 1].plot(psi_s, q_s, **kw)
    axs[1, 1].axhline(1.0, ls=':', c='k', alpha=0.5)
    axs[1, 1].set_ylabel('q')
    axs[1, 2].plot(psi_s, jtor_s / 1e6, **kw); axs[1, 2].set_ylabel(r'$\langle j_\phi \rangle$ [MA/m²]')

    for ax in axs.flat:
        ax.grid(ls=':')
        ax.set_xlabel(psi_lbl)

    fig.tight_layout()
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    fig.text(0.99, 0.01, timestamp, ha='right', va='bottom',
             fontsize=7, color='gray', transform=fig.transFigure)
    out_path = os.path.join(script_dir, f'profiles_{tag_s}{option_tag}_{timestamp}.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved: {out_path}")

# ---------------------------------------------------------------------------
# Internal-compare: isolate-edge vs parameterize-edge (internal solver only)
# ---------------------------------------------------------------------------
if run_internal_compare:

    def _run_internal_case(label, isolate_edge, parameterize_edge):
        print(f"\n{'='*60}")
        print(f"  INTERNAL [{label}]  (jphi-split-bootstrap / Fortran-embedded)")
        print(f"{'='*60}")
        build_gs(mygs)

        t0 = time.perf_counter()
        ne_c, Te_c, ni_c, Ti_c = ne, Te, ni, Ti
        pressure_c = pressure
        pp_vals = np.gradient(pressure_c) / np.gradient(psi_sample)
        # Enforce P' edge condition
        pp_vals[-1] = 0.0

        mygs.set_kinetic_profiles(
            te_prof={'type': 'linterp', 'x': psi_sample, 'y': Te_c / 1e3},
            ne_prof={'type': 'linterp', 'x': psi_sample, 'y': ne_c},
            ti_prof={'type': 'linterp', 'x': psi_sample, 'y': Ti_c / 1e3},
            ni_prof={'type': 'linterp', 'x': psi_sample, 'y': ni_c},
            Zeff=Zeff_val,
        )
        mygs.set_boot_ops(
            isolate_edge_jBS=isolate_edge,
            parameterize_jBS=parameterize_edge,
            scale_jBS=args.scale_jBS,
            Zeff=Zeff_val,
            diagnose_bs=args.diagnose_bs,
            taper_edge_jBS=args.taper_edge,
            taper_edge_psi0=args.taper_edge_psi0,
            taper_edge_shape=args.taper_edge_shape,
        )
        mygs.set_profiles(
            ffp_prof={'type': 'jphi-split-bootstrap', 'x': psi_sample, 'y': inductive_jphi},
            pp_prof={'type': 'linterp', 'x': psi_sample, 'y': pp_vals / pp_vals[0]},
            foffset=5.3 * 6.2,
        )
        mygs.set_targets(Ip=Ip_target, pax=float(pressure_c[0]))
        mygs.settings.pm = True
        mygs.update_settings()
        mygs.solve()

        elapsed = time.perf_counter() - t0
        stats = mygs.get_stats()
        coils, _ = mygs.get_coil_currents()
        print_result(f"INTERNAL [{label}]", elapsed, stats, coils)

        psi_p, F_p, Fp_p, P_p, Pp_p = mygs.get_profiles(npsi=n_sample_plot, psi_pad=1e-3)
        _, q_p, ravgs_p, _, _, _     = mygs.get_q(npsi=n_sample_plot, psi_pad=1e-3)
        jtor_p = F_p * Fp_p * ravgs_p[1] / mu0 + Pp_p * ravgs_p[0]
        mygs.reset()
        return elapsed, stats, psi_p, F_p, Fp_p, P_p, Pp_p, q_p, jtor_p

    (elapsed_iso, stats_iso,
     psi_iso, F_iso, Fp_iso, P_iso, Pp_iso, q_iso, jtor_iso) = \
        _run_internal_case('isolate-edge', isolate_edge=True, parameterize_edge=False)

    (elapsed_par, stats_par,
     psi_par, F_par, Fp_par, P_par, Pp_par, q_par, jtor_par) = \
        _run_internal_case('parameterize-edge', isolate_edge=False, parameterize_edge=True)

    # --- Scalar comparison table ---
    cmp_keys = [('Ip', 'Ip [MA]', 1e6), ('P_ax', 'Pax [kPa]', 1e3),
                ('beta_pol', 'beta_pol [%]', 1), ('beta_tor', 'beta_tor [%]', 1)]
    print(f"\n{'='*60}")
    print("  Comparison: isolate-edge vs parameterize-edge (internal)")
    print(f"{'='*60}")
    print(f"  {'Quantity':<18} {'isolate-edge':>14} {'param-edge':>12} {'Diff %':>10}")
    print(f"  {'-'*58}")
    for key, label, scale in cmp_keys:
        v_iso = stats_iso.get(key, None)
        v_par = stats_par.get(key, None)
        if v_iso is None or v_par is None:
            continue
        pct = 100.0 * (v_par - v_iso) / v_iso if v_iso != 0 else float('nan')
        print(f"  {label:<18} {v_iso/scale:>14.4f} {v_par/scale:>12.4f} {pct:>9.3f}%")
    print(f"\n  Timing:  isolate-edge {elapsed_iso:.3f} s  |  "
          f"parameterize-edge {elapsed_par:.3f} s"
          f"  (ratio {elapsed_iso/elapsed_par:.2f}x)")

    # --- Profile comparison figure (3 rows × 3 cols) ---
    def _pct_diff_ic(a, b):
        denom = np.where(np.abs(a) > 1e-10 * np.nanmax(np.abs(a)), np.abs(a), np.nan)
        return 100.0 * (b - a) / denom

    ic_parts = []
    option_tag_ic   = ''
    option_label_ic = ''

    fig, axs = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle('Internal bootstrap: isolate-edge vs parameterize-edge' + option_label_ic, fontsize=13)
    psi_lbl = r'$\hat{\psi}$'
    kw = dict(lw=1.5)

    # Row 0: F, F·F', P
    axs[0, 0].plot(psi_iso, F_iso,         label='isolate-edge', **kw)
    axs[0, 0].plot(psi_par, F_par,         label='parameterize-edge', ls='--', **kw)
    axs[0, 0].set_ylabel('F [T·m]');  axs[0, 0].legend(fontsize=8)

    axs[0, 1].plot(psi_iso, F_iso * Fp_iso, **kw)
    axs[0, 1].plot(psi_par, F_par * Fp_par, ls='--', **kw)
    axs[0, 1].set_ylabel("F·F' [T²]")

    axs[0, 2].plot(psi_iso, P_iso / 1e6, **kw)
    axs[0, 2].plot(psi_par, P_par / 1e6, ls='--', **kw)
    axs[0, 2].set_ylabel('P [MPa]')

    # Row 1: P', q, jφ
    axs[1, 0].plot(psi_iso, Pp_iso / 1e3, **kw)
    axs[1, 0].plot(psi_par, Pp_par / 1e3, ls='--', **kw)
    axs[1, 0].set_ylabel("P' [kPa]")

    axs[1, 1].plot(psi_iso, q_iso, **kw)
    axs[1, 1].plot(psi_par, q_par, ls='--', **kw)
    axs[1, 1].axhline(1.0, ls=':', c='k', alpha=0.5)
    axs[1, 1].set_ylabel('q')

    axs[1, 2].plot(psi_iso, jtor_iso / 1e6, **kw)
    axs[1, 2].plot(psi_par, jtor_par / 1e6, ls='--', **kw)
    axs[1, 2].set_ylabel(r'$\langle j_\phi \rangle$ [MA/m²]')

    # Row 2: percentage differences (parameterize-edge − isolate-edge)
    for col, (yi, yp, lbl) in enumerate([
            (F_iso * Fp_iso, F_par * Fp_par, "Δ F·F' [%]"),
            (q_iso,          q_par,          'Δ q [%]'),
            (jtor_iso,       jtor_par,       r'Δ $j_\phi$ [%]'),
    ]):
        axs[2, col].plot(psi_iso, _pct_diff_ic(yi, yp), c='tab:red', **kw)
        axs[2, col].axhline(0, ls=':', c='k', alpha=0.5)
        axs[2, col].set_ylabel(lbl)

    for ax in axs.flat:
        ax.grid(ls=':')
        ax.set_xlabel(psi_lbl)

    fig.tight_layout()
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    fig.text(0.99, 0.01, timestamp, ha='right', va='bottom',
             fontsize=7, color='gray', transform=fig.transFigure)
    out_path = os.path.join(script_dir, f'profiles_internal_compare{option_tag_ic}_{timestamp}.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved: {out_path}")

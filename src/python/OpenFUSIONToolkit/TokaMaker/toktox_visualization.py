"""TokTox Visualization Module

All plotting, movie generation, and interactive visualization for TokTox simulations.
Functions accept a TokTox object and produce plots that can be displayed inline or saved to file.
"""

from __future__ import annotations
import os
import subprocess
import tempfile
import shutil
import platform
# from typing import TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize

# if TYPE_CHECKING:
#     from toktox import TokTox

# =========================================================================
#  Style constants
# =========================================================================
COLOR_TM = 'steelblue'
COLOR_TX = 'crimson'
LS_PRI = '-'
LS_SEC = '--'
MK_TM = '.'
MK_SZ = 3
LW = 1.6
COLORS_MULTI = [
    'darkorange', 'forestgreen', 'mediumpurple',
    'goldenrod', 'deeppink', 'teal', 'sienna',
]
VLINE_COLOR = 'black'
VLINE_LS = ':'
VLINE_LW = 1.0
GRID_ALPHA = 0.2
LEGEND_FS = 9
TITLE_FS = 13
LABEL_FS = 11
TICK_FS = 10
INFO_FS = 13
DIAG_FS = 13

MOVIE_FIG_W, MOVIE_FIG_H = 19.2, 10.8
MOVIE_DPI = 200


# =========================================================================
#  Helpers
# =========================================================================

def _in_jupyter():
    """Return True if running inside a Jupyter notebook."""
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        return shell in ('ZMQInteractiveShell',)
    except Exception:
        return False


def _save_or_display(fig, save_path=None, display=True):
    """Save figure and/or display it, then close."""
    if save_path is not None:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    if display:
        plt.show()
    else:
        plt.close(fig)


def _style(ax):
    ax.grid(True, alpha=GRID_ALPHA)
    ax.tick_params(labelsize=TICK_FS)


def _legend_if_labeled(ax, *args, **kwargs):
    """Add a legend only when labeled artists are present."""
    handles, labels = ax.get_legend_handles_labels()
    if handles and labels:
        ax.legend(*args, **kwargs)


def _vline(axes, t_now):
    for ax in (axes if hasattr(axes, '__iter__') else [axes]):
        ax.axvline(t_now, color=VLINE_COLOR, ls=VLINE_LS, lw=VLINE_LW, zorder=10)


def _tx_scalar(tt, var_name, scale=1.0):
    """Return (times, values) arrays from data_tree.scalars at full TORAX resolution."""
    dt = getattr(tt, '_data_tree', None)
    if dt is None:
        return None, None
    var = getattr(dt.scalars, var_name)
    return var.coords['time'].values, var.to_numpy() * scale


def _tx_profile_at_rho(tt, var_name, rho_val, rho_coord='rho_norm', scale=1.0):
    """Return (times, values) for a profile variable at a fixed rho, full resolution."""
    dt = getattr(tt, '_data_tree', None)
    if dt is None:
        return None, None
    var = getattr(dt.profiles, var_name)
    sliced = var.sel(**{rho_coord: rho_val}, method='nearest')
    return sliced.coords['time'].values, sliced.to_numpy() * scale


def _prof(state_dict, idx):
    """Return (x, y) arrays from a profile state dict, or (None, None)."""
    d = state_dict.get(idx)
    if d is None:
        return None, None
    return d['x'], d['y']


def _make_temp_dir():
    """Create temp directory: RAM-backed on Linux, OS default elsewhere."""
    if platform.system() == 'Linux' and os.path.isdir('/dev/shm'):
        return tempfile.mkdtemp(prefix='toktox_viz_', dir='/dev/shm')
    return tempfile.mkdtemp(prefix='toktox_viz_')


def _x_points_active(tt, i, t=None):
    """Return True when X-point targets should be applied at timestep index i."""
    diverted = getattr(tt, '_diverted_times', None)
    if diverted is None:
        return False

    div_arr = np.asarray(diverted)

    # New API: diverted window defined as (t_start, t_end).
    if (
        div_arr.ndim == 1
        and div_arr.size == 2
        and np.issubdtype(div_arr.dtype, np.number)
        and not np.issubdtype(div_arr.dtype, np.bool_)
    ):
        if t is None:
            times = getattr(tt, '_tm_times', None)
            if times is None or i >= len(times):
                return False
            t = times[i]
        return float(div_arr[0]) <= float(t) <= float(div_arr[1])

    # Backward compatibility: per-timestep diverted mask.
    if div_arr.ndim == 1 and i < div_arr.size:
        return bool(div_arr[i])

    return False


# =========================================================================
#  Profile plot (per-timestep diagnostic plot)
# =========================================================================

def profile_plot(tt, i, t, save_path=None, display=True):
    """Detailed profile comparison at a single timestep."""
    from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk

    s = tt._state
    psi_N = tt._psi_N

    tm_psi, tm_f_prof, tm_fp_prof, tm_p_prof, tm_pp_prof = tt._tm.get_profiles(npsi=len(tt._psi_N))
    tm_ffp_prof = tm_f_prof * tm_fp_prof

    fig, axes = plt.subplots(6, 3, figsize=(20, 24))
    plt.suptitle(f'loop {tt._current_loop} - t-idx {i}/{len(tt._tm_times)-1} - t = {t:.1f} s', fontsize=14)

    # Row 0: p' and p comparison
    ax = axes[0, 0]
    ax.set_title("p' and p comparison")
    ax.plot(s['pp_prof_tx'][i]['x'], s['pp_prof_tx'][i]['y'], 'b-', label="p' TX", linewidth=2)
    ax.plot(psi_N, s['pp_prof_tm'][i]['y'], 'b--', label="p' TM", linewidth=2)
    ax.set_ylabel("p' [Pa/Wb]", color='b')
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.tick_params(axis='y', labelcolor='b')
    ax.legend(fontsize=9, loc='upper left')
    ax2 = ax.twinx()
    ax2.plot(s['p_prof_tm'][i]['x'], s['p_prof_tm'][i]['y'], 'r-', label='p TM', linewidth=2)
    ax2.plot(s['ptot'][i]['x'], s['ptot'][i]['y'], 'r--', label='p TX', linewidth=2)
    ax2.set_ylabel('p [Pa]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.legend(fontsize=9, loc='upper right')
    axes[0, 1].axis('off')
    axes[0, 2].axis('off')

    # Row 1: Current densities, resistivity, FF'
    ax = axes[1, 0]
    ax.set_title('Current densities')
    ax.plot(s['j_tot'][i]['x'], s['j_tot'][i]['y'] / 1e6, 'k--', label=r'$j_{tot}$', linewidth=2)
    ax.plot(s['j_ohmic'][i]['x'], s['j_ohmic'][i]['y'] / 1e6, 'r-', label=r'$j_{ohmic}$', linewidth=1.5)
    if i in s.get('j_ni', {}):
        ax.plot(s['j_ni'][i]['x'], s['j_ni'][i]['y'] / 1e6, 'b--', label=r'$j_{NI}$', linewidth=1.5)
    if i in s.get('j_bootstrap', {}):
        ax.plot(s['j_bootstrap'][i]['x'], s['j_bootstrap'][i]['y'] / 1e6, 'g-', label=r'$j_{bootstrap}$', linewidth=1.5)
    if i in s.get('j_generic_current', {}):
        ax.plot(s['j_generic_current'][i]['x'], s['j_generic_current'][i]['y'] / 1e6,
                color='darkorange', linestyle='-', label=r'$j_{gen}$', linewidth=1.5)
    ax.set_ylabel(r'$j$ [MA/m²]')
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    ax.set_title('Resistivity')
    ax.plot(s['eta_prof'][i]['x'], s['eta_prof'][i]['y'], 'r-', label='TX', linewidth=2)
    ax.set_yscale('log')
    ax.set_ylabel(r'$\eta$ [Ohm m]')
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.legend(fontsize=9)

    ax = axes[1, 2]
    ax.set_title("FF' comparison")
    ax.plot(psi_N, s['ffp_prof_tx'][i]['y'], 'k-', label="FF' total TX", linewidth=2)
    ax.plot(s['ffp_ni_prof'][i]['x'], s['ffp_ni_prof'][i]['y'], 'b-', label="FF' NI")
    ax.plot(psi_N, s['ffp_prof_tx'][i]['y'] - s['ffp_ni_prof'][i]['y'], 'r--', label="FF' inductive")
    ax.plot(psi_N, s['ffp_prof_tm'][i]['y'], 'g--', label="FF' total TM", linewidth=1)
    ax.set_ylabel("FF'")
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.legend(fontsize=9)

    # Row 2: Psi, <1/R>, Volume
    ax = axes[2, 0]
    ax.set_title('Psi profile comparison')
    ax.plot(s['psi_tx'][i]['x'], s['psi_tx'][i]['y'], 'b-', label='Psi TX', linewidth=2)
    ax.plot(s['psi_tm'][i]['x'], s['psi_tm'][i]['y'], 'r--', label='Psi TM', linewidth=2)
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\psi$ [Wb/rad]')
    ax.legend(fontsize=9)

    ax = axes[2, 1]
    ax.set_title('<1/R> comparison')
    if i in s['R_inv_avg_tm']:
        ax.plot(s['R_inv_avg_tm'][i]['x'], s['R_inv_avg_tm'][i]['y'], 'r-', label='<1/R> TM', linewidth=2)
    if i in s['R_inv_avg_tx']:
        ax.plot(s['R_inv_avg_tx'][i]['x'], s['R_inv_avg_tx'][i]['y'], 'b--', label='<1/R> TX', linewidth=2)
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel('<1/R> [1/m]')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Row 3: q, T, n
    ax = axes[3, 0]
    ax.set_title('q profile')
    ax.plot(s['q_prof_tx'][i]['x'], s['q_prof_tx'][i]['y'], 'b--', label='TX', linewidth=1)
    ax.plot(s['q_prof_tm'][i]['x'], s['q_prof_tm'][i]['y'], 'r--', label='TM', linewidth=2)
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel('q')
    ax.legend(fontsize=9)

    ax = axes[3, 1]
    ax.set_title('T_e and T_i')
    if i in s.get('T_e', {}) and i in s.get('T_i', {}):
        ax.plot(s['T_e'][i]['x'], s['T_e'][i]['y'], 'r-', label=r'$T_e$')
        ax.plot(s['T_i'][i]['x'], s['T_i'][i]['y'], 'm--', label=r'$T_i$')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel('T [keV]')
        ax.legend(fontsize=9)
    else:
        ax.text(0.5, 0.5, 'No T profiles', ha='center', va='center')

    ax = axes[3, 2]
    ax.set_title('n_e and n_i')
    if i in s.get('n_e', {}) and i in s.get('n_i', {}):
        ax.plot(s['n_e'][i]['x'], s['n_e'][i]['y'], 'b-', label=r'$n_e$')
        ax.plot(s['n_i'][i]['x'], s['n_i'][i]['y'], 'c--', label=r'$n_i$')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$n$ [m$^{-3}$]')
        ax.legend(fontsize=9)
    else:
        ax.text(0.5, 0.5, 'No n profiles', ha='center', va='center')

    # Row 4: Chi profiles
    ax = axes[4, 0]
    ax.set_title(r'$\chi$ (NEO)')
    for key, label, color in [('chi_neo_e', r'$\chi_{NEO,e}$', 'r'), ('chi_neo_i', r'$\chi_{NEO,i}$', 'b')]:
        try:
            if i in s[key]:
                ax.plot(s[key][i]['x'], s[key][i]['y'], color=color, label=label, linewidth=1.5)
        except (KeyError, TypeError):
            pass
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\chi$ [m²/s]')
    _legend_if_labeled(ax, fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[4, 1]
    ax.set_title(r'$\chi$ (Turbulent)')
    for key, label, color in [('chi_etg_e', r'$\chi_{ETG,e}$', 'g'), ('chi_itg_e', r'$\chi_{ITG,e}$', 'm'),
                               ('chi_itg_i', r'$\chi_{ITG,i}$', 'c'), ('chi_tem_e', r'$\chi_{TEM,e}$', 'orange'),
                               ('chi_tem_i', r'$\chi_{TEM,i}$', 'purple')]:
        try:
            if i in s[key]:
                ax.plot(s[key][i]['x'], s[key][i]['y'], color=color, label=label, linewidth=1.5)
        except (KeyError, TypeError):
            pass
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\chi$ [m²/s]')
    _legend_if_labeled(ax, fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    ax = axes[4, 2]
    ax.set_title(r'$\chi$ (Turbulent - e/i)')
    for key, label, color, ls in [('chi_turb_e', r'$\chi_{turb,e}$', 'b', '-'), ('chi_turb_i', r'$\chi_{turb,i}$', 'r', '--')]:
        try:
            if i in s[key]:
                ax.plot(s[key][i]['x'], s[key][i]['y'], color=color, ls=ls, label=label, linewidth=1.5)
        except (KeyError, TypeError):
            pass
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\chi$ [m²/s]')
    _legend_if_labeled(ax, fontsize=9)
    ax.grid(True, alpha=0.3)

    # Row 5: Diffusivity profiles
    ax = axes[5, 0]
    ax.set_title(r'$D_{ITG}$')
    try:
        if i in s['D_itg_e']:
            ax.plot(s['D_itg_e'][i]['x'], s['D_itg_e'][i]['y'], 'g-', label=r'$D_{ITG,e}$', linewidth=1.5)
    except (KeyError, TypeError):
        pass
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$D$ [m²/s]')
    _legend_if_labeled(ax, fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[5, 1]
    ax.set_title(r'$D$ (NEO & TEM)')
    for key, label, color, ls in [('D_neo_e', r'$D_{NEO,e}$', 'r', '-'), ('D_tem_e', r'$D_{TEM,e}$', 'm', '--')]:
        try:
            if i in s[key]:
                ax.plot(s[key][i]['x'], s[key][i]['y'], color=color, ls=ls, label=label, linewidth=1.5)
        except (KeyError, TypeError):
            pass
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$D$ [m²/s]')
    _legend_if_labeled(ax, fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[5, 2]
    ax.set_title(r'$D_{turb}$')
    try:
        if i in s['D_turb_e']:
            ax.plot(s['D_turb_e'][i]['x'], s['D_turb_e'][i]['y'], 'c-', label=r'$D_{turb,e}$', linewidth=1.5)
    except (KeyError, TypeError):
        pass
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$D$ [m²/s]')
    _legend_if_labeled(ax, fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, hspace=0.3, wspace=0.35)
    _save_or_display(fig, save_path, display)


# =========================================================================
#  TokaMaker diagnostic plot
# =========================================================================

def tm_diagnostic_plot(tt, i, t, level_attempts, solve_succeeded, save_path=None, display=True):
    """TokaMaker input/output diagnostic plot for a single timestep."""
    # Seed FF'/p'/q references are pulled from tt._state.
    # Keep the old EQDSK reader path commented for future volume-vs-psi use.
    # from read_eqdsk_extended import read_eqdsk_extended

    s = tt._state
    psi_N = tt._psi_N

    _winning = next((a for a in level_attempts if a['succeeded']), None)
    _last = level_attempts[-1] if level_attempts else {}
    fail_msg = _last.get('error') if not solve_succeeded else None

    _level_colors = plt.cm.tab20.colors

    def _plot_levels(ax, key, seed_x=None, seed_y=None, seed_label=None):
        for attempt in level_attempts:
            color = _level_colors[attempt['level'] % len(_level_colors)]
            y_data = attempt[key]['y'].copy()
            if y_data[0] != 0:
                y_data = y_data / y_data[0]
            if attempt['succeeded']:
                ax.plot(attempt[key]['x'], y_data, color=color, linewidth=2.5, zorder=5,
                        label=f"Level {attempt['level']}: {attempt['name']} \u2713",
                        marker='o', markersize=3, markevery=max(len(attempt[key]['x']) // 10, 1), alpha=0.8)
            else:
                ax.plot(attempt[key]['x'], y_data, color=color, linewidth=1.2, linestyle='--', alpha=0.6,
                        label=f"Level {attempt['level']}: {attempt['name']} \u2717")
        if seed_x is not None:
            ax.plot(seed_x, seed_y, 'k--', linewidth=1.5, alpha=0.7, label=seed_label)

    def render_table(ax, rows, title):
        ax.axis('off')
        ax.set_title(title, fontsize=10, fontweight='bold', pad=4)
        tbl = ax.table(cellText=rows[1:], colLabels=rows[0], loc='center', cellLoc='left',
                       bbox=[0.0, 0.0, 1.0, 0.92])
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(9)
        tbl.scale(1, 1.5)
        for (row, col), cell in tbl.get_celld().items():
            if row == 0:
                cell.set_facecolor('#d0e4f7')
            elif row % 2 == 0:
                cell.set_facecolor('#f5f5f5')

    # Scalar values
    Ip_seed = abs(float(tt._Ip_seed[i]))
    pax_seed = abs(float(tt._pax_seed[i]))
    psi_lcfs_seed = float(tt._psi_lcfs_seed[i])
    Ip_tx = abs(float(s['Ip'][i]))
    pax_tx = abs(float(s['pax'][i]))
    psi_lcfs_tx = float(s['psi_lcfs_tx'][i])
    q95_tx = float(s['q95'][i])
    q0_tx = float(s['q0'][i])
    vloop_tx = float(s['vloop_tx'][i])
    beta_pol_tx = float(s['beta_pol'][i])
    beta_n_tx = float(s['beta_N_tx'][i])

    _seed_ffp_prof = s.get('ffp_prof_eqdsk', {}).get(i)
    _seed_pp_prof = s.get('pp_prof_eqdsk', {}).get(i)
    _seed_q_prof = s.get('q_prof_eqdsk', {}).get(i)

    _seed_ffp_x = _seed_ffp_norm = None
    if _seed_ffp_prof is not None:
        _seed_ffp_x = np.asarray(_seed_ffp_prof.get('x', []), dtype=float)
        _seed_ffp_norm = np.asarray(_seed_ffp_prof.get('y', []), dtype=float)
        if _seed_ffp_x.size == 0 or _seed_ffp_norm.size == 0:
            _seed_ffp_x, _seed_ffp_norm = None, None

    _seed_pp_x = _seed_pp_norm = None
    if _seed_pp_prof is not None:
        _seed_pp_x = np.asarray(_seed_pp_prof.get('x', []), dtype=float)
        _seed_pp_norm = np.asarray(_seed_pp_prof.get('y', []), dtype=float)
        if _seed_pp_x.size == 0 or _seed_pp_norm.size == 0:
            _seed_pp_x, _seed_pp_norm = None, None

    q0_seed_eq = np.nan
    q95_seed_eq = np.nan
    if _seed_q_prof is not None:
        _seed_q_x = np.asarray(_seed_q_prof.get('x', []), dtype=float)
        _seed_q_y = np.asarray(_seed_q_prof.get('y', []), dtype=float)
        if _seed_q_x.size > 0 and _seed_q_y.size > 0:
            q0_seed_eq = float(np.interp(0.0, _seed_q_x, _seed_q_y))
            q95_seed_eq = float(np.interp(0.95, _seed_q_x, _seed_q_y))

    # Fallback kept intentionally (commented):
    # _eqtimes_arr = np.array(tt._eqtimes)
    # _seed_idx = int(np.argmin(np.abs(_eqtimes_arr - t)))
    # _seed_g = read_eqdsk_extended(tt._init_files[_seed_idx])
    # _seed_ffp_norm = _seed_g['ffprim'].copy() / _seed_g['ffprim'][0]
    # _seed_pp_norm = _seed_g['pprime'].copy() / _seed_g['pprime'][0]
    # q0_seed_eq = _seed_g['q0']
    # q95_seed_eq = _seed_g['q95']
    # _seed_volume_vs_psi = _seed_g[...]  # volume(psi) if/when needed

    fig = plt.figure(figsize=(22, 12))
    gs_layout = fig.add_gridspec(3, 6, hspace=0.70, wspace=0.55)

    ax_ffp_tx = fig.add_subplot(gs_layout[0, 0:2])
    ax_pp_tx = fig.add_subplot(gs_layout[1, 0:2])
    ax_eta = fig.add_subplot(gs_layout[2, 0:2])
    ax_ffp_tm = fig.add_subplot(gs_layout[0, 2:4])
    ax_pp_tm = fig.add_subplot(gs_layout[1, 2:4])
    ax_q_tm = fig.add_subplot(gs_layout[2, 2:4])
    ax_tbl1 = fig.add_subplot(gs_layout[0, 4:6])
    ax_tbl2 = fig.add_subplot(gs_layout[1:3, 4:6])

    _plot_levels(ax_ffp_tx, 'ffp', seed_x=_seed_ffp_x, seed_y=_seed_ffp_norm, seed_label="FF' seed EQDSK (norm)")
    ax_ffp_tx.set_title("FF' tried levels (norm)", fontsize=10)
    ax_ffp_tx.set_xlabel(r'$\hat{\psi}$')
    ax_ffp_tx.set_ylabel("FF' (norm)")
    ax_ffp_tx.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2)
    ax_ffp_tx.grid(True, alpha=0.3)
    ax_ffp_tx.axhline(0, color='k', linewidth=0.5)

    _plot_levels(ax_pp_tx, 'pp', seed_x=_seed_pp_x, seed_y=_seed_pp_norm, seed_label="p' seed EQDSK (norm)")
    ax_pp_tx.set_title("p' tried levels (normalized)", fontsize=10)
    ax_pp_tx.set_xlabel(r'$\hat{\psi}$')
    ax_pp_tx.set_ylabel("p' (norm)")
    ax_pp_tx.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2)
    ax_pp_tx.grid(True, alpha=0.3)

    ax_eta.plot(s['eta_prof'][i]['x'], s['eta_prof'][i]['y'], 'r-', linewidth=2)
    ax_eta.set_yscale('log')
    ax_eta.set_title(r'Resistivity $\eta$ (input to TM)', fontsize=10)
    ax_eta.set_xlabel(r'$\hat{\psi}$')
    ax_eta.set_ylabel(r'$\eta$ [$\Omega\cdot$m]')
    ax_eta.grid(True, alpha=0.3)

    if solve_succeeded:
        input_rows = [
            ['Parameter', 'Init EQDSK', 'TORAX', 'TokaMaker'],
            ['Ip', f'{Ip_seed / 1e6:.3f} MA', f'{Ip_tx / 1e6:.3f} MA', f'{float(s["Ip_tm"][i]) / 1e6:.3f} MA'],
            ['pax', f'{pax_seed / 1e3:.2f} kPa', f'{pax_tx / 1e3:.2f} kPa', f'{float(s["pax_tm"][i]) / 1e3:.2f} kPa'],
            ['psi_lcfs', f'{psi_lcfs_seed:.4f} Wb/rad', f'{psi_lcfs_tx:.4f} Wb/rad', f'{float(s["psi_lcfs_tm"][i]):.4f} Wb/rad'],
        ]
        diag_rows = [
            ['Parameter', 'Init EQDSK', 'TORAX', 'TokaMaker'],
            ['q95', f'{q95_seed_eq:.3f}', f'{q95_tx:.3f}', f'{float(s["q95_tm"][i]):.3f}'],
            ['q0', f'{q0_seed_eq:.3f}', f'{q0_tx:.3f}', f'{float(s["q0_tm"][i]):.3f}'],
            ['v_loop', '\u2014', f'{vloop_tx:.3f} V', f'{float(s["vloop_tm"][i]):.3f} V'],
            ['beta_pol', '\u2014', f'{beta_pol_tx:.4f}', '\u2014'],
            ['beta_N', '\u2014', f'{beta_n_tx:.4f}', f'{float(s["beta_N_tm"][i]):.4f}'],
            ['l_i', '\u2014', '\u2014', f'{float(s["l_i_tm"][i]):.4f}'],
        ]

        ax_ffp_tm.plot(s['ffp_prof_tx'][i]['x'], s['ffp_prof_tx'][i]['y'], 'b--', linewidth=1.5, label="FF' TX (real)")
        ax_ffp_tm.plot(s['ffp_prof_tm'][i]['x'], s['ffp_prof_tm'][i]['y'], 'r-', linewidth=2, label="FF' TM (real)")
        ax_ffp_tm.set_title("FF' TM output vs TX input", fontsize=10)
        ax_ffp_tm.set_xlabel(r'$\hat{\psi}$')
        ax_ffp_tm.set_ylabel("FF'")
        ax_ffp_tm.legend(fontsize=8)
        ax_ffp_tm.grid(True, alpha=0.3)
        ax_ffp_tm.axhline(0, color='k', linewidth=0.5)

        ax_pp_tm.plot(s['pp_prof_tx'][i]['x'], s['pp_prof_tx'][i]['y'], 'b--', linewidth=1.5, label="p' TX (real)")
        ax_pp_tm.plot(s['pp_prof_tm'][i]['x'], s['pp_prof_tm'][i]['y'], 'r-', linewidth=2, label="p' TM (real)")
        ax_pp_tm.set_title("p' TM output vs TX input", fontsize=10)
        ax_pp_tm.set_xlabel(r'$\hat{\psi}$')
        ax_pp_tm.set_ylabel("p'")
        ax_pp_tm.legend(fontsize=8)
        ax_pp_tm.grid(True, alpha=0.3)

        if i in s['p_prof_tx'] and i in s['p_prof_tm']:
            ax_pp_tm_2 = ax_pp_tm.twinx()
            ax_pp_tm_2.plot(s['p_prof_tx'][i]['x'], s['p_prof_tx'][i]['y'], 'c--', linewidth=1.5, alpha=0.7, label='p TX')
            ax_pp_tm_2.plot(s['p_prof_tm'][i]['x'], s['p_prof_tm'][i]['y'], 'orange', linewidth=2, alpha=0.7, label='p TM')
            ax_pp_tm_2.set_ylabel('p [Pa]', color='orange')
            ax_pp_tm_2.tick_params(axis='y', labelcolor='orange')
            lines1, labels1 = ax_pp_tm.get_legend_handles_labels()
            lines2, labels2 = ax_pp_tm_2.get_legend_handles_labels()
            ax_pp_tm.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc='upper left')

        try:
            psi_geo, q_tm_vals, _, _, _, _ = tt._tm.get_q(npsi=len(tt._psi_N), psi_pad=0.02)
            ax_q_tm.plot(psi_geo, q_tm_vals, 'r--', linewidth=2, label='TokaMaker')
        except Exception:
            pass
        ax_q_tm.plot(s['q_prof_tx'][i]['x'], s['q_prof_tx'][i]['y'], 'b-', linewidth=2, label='TORAX')
        ax_q_tm.set_title('q profile (TM vs TORAX)', fontsize=10)
        ax_q_tm.set_xlabel(r'$\hat{\psi}$')
        ax_q_tm.set_ylabel('q')
        ax_q_tm.legend(fontsize=8)
        ax_q_tm.grid(True, alpha=0.3)

        render_table(ax_tbl1, input_rows, 'Scalar Inputs vs TokaMaker Outputs')
        render_table(ax_tbl2, diag_rows, 'TORAX Diagnostics vs TokaMaker')

        plt.suptitle(
            f'TM Diagnostic \u2014 loop {tt._current_loop}, t-idx {i}/{len(tt._tm_times) - 1}, t = {t:.2f} s'
            f'  |  TokaMaker: SUCCESS',
            fontsize=13, color='darkgreen',
        )
    else:
        input_rows = [
            ['Parameter', 'Init EQDSK', 'TORAX'],
            ['Ip', f'{Ip_seed / 1e6:.3f} MA', f'{Ip_tx / 1e6:.3f} MA'],
            ['pax', f'{pax_seed / 1e3:.2f} kPa', f'{pax_tx / 1e3:.2f} kPa'],
            ['psi_lcfs', f'{psi_lcfs_seed:.4f} Wb/rad', f'{psi_lcfs_tx:.4f} Wb/rad'],
        ]
        diag_rows = [
            ['Parameter', 'Init EQDSK', 'TORAX'],
            ['q95', f'{q95_seed_eq:.3f}', f'{q95_tx:.3f}'],
            ['q0', f'{q0_seed_eq:.3f}', f'{q0_tx:.3f}'],
            ['v_loop', '\u2014', f'{vloop_tx:.3f} V'],
            ['beta_pol', '\u2014', f'{beta_pol_tx:.4f}'],
            ['beta_N', '\u2014', f'{beta_n_tx:.4f}'],
        ]

        for blank_ax in [ax_ffp_tm, ax_pp_tm, ax_q_tm]:
            blank_ax.axis('off')
            blank_ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=11, color='gray',
                          transform=blank_ax.transAxes, fontweight='bold')

        render_table(ax_tbl1, input_rows, 'Scalar Inputs (Init EQDSK vs TORAX)')
        diag_rows.append(['', '', ''])
        diag_rows.append(['Failure Reason:', fail_msg if fail_msg else 'Unknown', ''])
        render_table(ax_tbl2, diag_rows, 'TORAX Diagnostics & Failure Info')

        plt.suptitle(
            f'TM Diagnostic \u2014 loop {tt._current_loop}, t-idx {i}/{len(tt._tm_times) - 1}, t = {t:.2f} s'
            f'  |  TokaMaker: FAILED',
            fontsize=13, color='darkred', fontweight='bold',
        )

    _save_or_display(fig, save_path, display)


# =========================================================================
#  TM loop summary plot
# =========================================================================

def tm_loop_summary_plot(tt, loop_level_log, save_path=None, display=True):
    """Summary figure showing per-timestep GS solve outcomes."""
    n = len(loop_level_log)
    if n == 0:
        return

    fig = plt.figure(figsize=(10, max(5, 0.4 * n + 3.5)))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.85, 0.15], hspace=0.4)

    ax_table = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])
    ax_table.axis('off')
    ax_legend.axis('off')

    col_labels = ['t-idx', 't (s)', 'Result']
    rows = []
    cell_colors = []
    for entry in loop_level_log:
        if entry['succeeded']:
            result = f"Lvl {entry['level']}"
            row_color = ['#d4edda'] * 3
        else:
            error_msg = entry['error'] if entry['error'] else 'Unknown error'
            result = error_msg[:47] + '...' if len(error_msg) > 50 else error_msg
            row_color = ['#f8d7da'] * 3
        rows.append([str(entry['i']), f"{entry['t']:.3f}", result])
        cell_colors.append(row_color)

    tbl = ax_table.table(cellText=rows, colLabels=col_labels, cellColours=cell_colors,
                         loc='center', cellLoc='center', bbox=[0.0, 0.0, 1.0, 1.0])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.4)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor('#d0e4f7')
            cell.set_text_props(fontweight='bold')

    legend_text = ("Levels: Level 0: Raw TORAX profiles | Level 1: Sign-flip clipping | "
                   "Level 2: Pedestal smoothing | Level 3: Power-flux shape")
    ax_legend.text(0.5, 0.5, legend_text, ha='center', va='center', fontsize=9,
                   transform=ax_legend.transAxes,
                   bbox=dict(boxstyle='round,pad=0.8', facecolor='#f0f0f0', edgecolor='gray', linewidth=1))

    n_ok = sum(1 for e in loop_level_log if e['succeeded'])
    n_fail = n - n_ok
    plt.suptitle(
        f'GS loop {tt._current_loop} Summary \u2014 {n_ok}/{n} timesteps succeeded, {n_fail} failed',
        fontsize=12, fontweight='bold',
        color='darkgreen' if n_fail == 0 else ('darkred' if n_ok == 0 else 'darkorange'),
    )
    _save_or_display(fig, save_path, display)


# =========================================================================
#  Profile evolution plot
# =========================================================================

def plot_profile_evolution(tt, save_path=None, display=True, one_plot=False):
    """Plot profile evolution by pulse phase by default, or as one figure when one_plot=True."""
    s = tt._state
    times = np.array(tt._tm_times)
    if len(times) == 0:
        return

    ft = getattr(tt, '_flattop', np.zeros(len(times), dtype=bool)).astype(bool)

    if np.any(ft):
        ft_indices = np.where(ft)[0]
        ft_start_i = ft_indices[0]
        ft_end_i = ft_indices[-1]
        rampup_mask = np.arange(len(times)) < ft_start_i
        flattop_mask = ft
        rampdown_mask = np.arange(len(times)) > ft_end_i
    else:
        rampup_mask = np.ones(len(times), dtype=bool)
        flattop_mask = np.zeros(len(times), dtype=bool)
        rampdown_mask = np.zeros(len(times), dtype=bool)

    phases = [
        ('Ramp-up', rampup_mask),
        ('Flattop', flattop_mask),
        ('Ramp-down', rampdown_mask),
    ]
    if one_plot:
        phases = [('Pulse', np.ones(len(times), dtype=bool))]

    cmap = cm.plasma

    for phase_name, mask in phases:
        indices = np.where(mask)[0]
        if len(indices) == 0:
            continue

        phase_times = times[indices]
        t_min, t_max = phase_times[0], phase_times[-1]
        norm = Normalize(vmin=t_min, vmax=t_max if t_max > t_min else t_min + 1e-9)

        fig, axes = plt.subplots(2, 3, figsize=(14, 10))
        if one_plot:
            fig.suptitle(f'Profile Evolution Over Time (loop {tt._current_loop})', fontsize=14)
        else:
            fig.suptitle(f'Profile Evolution - {phase_name} (loop {tt._current_loop})', fontsize=14)

        plot_specs = [
            (axes[0, 0], 'n_e', r'$n_e$ [m$^{-3}$]', 'n_e'),
            (axes[1, 0], 'n_i', r'$n_i$ [m$^{-3}$]', 'n_i'),
            (axes[0, 1], 'T_e', r'$T_e$ [keV]', 'T_e'),
            (axes[1, 1], 'T_i', r'$T_i$ [keV]', 'T_i'),
            (axes[0, 2], 'ptot', r'$p$ [Pa]', 'p'),
        ]

        for ax, key, ylabel, title in plot_specs:
            ax.set_title(title)
            ax.set_xlabel(r'$\hat{\psi}$')
            ax.set_ylabel(ylabel)
            for i_t in indices:
                if i_t in s.get(key, {}):
                    color = cmap(norm(times[i_t]))
                    ax.plot(s[key][i_t]['x'], s[key][i_t]['y'], color=color, linewidth=1.5, alpha=0.8)
            ax.set_xlim([0, 1])

        ax = axes[1, 2]
        ax.set_title('q')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$q$')
        for i_t in indices:
            if i_t in s.get('q_prof_tx', {}):
                color = cmap(norm(times[i_t]))
                ax.plot(s['q_prof_tx'][i_t]['x'], s['q_prof_tx'][i_t]['y'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_xlim([0, 1])

        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.8, aspect=30, pad=0.02)
        cbar.set_label('Time [s]', fontsize=12)

        if save_path is not None and not one_plot:
            base, ext = os.path.splitext(save_path)
            phase_tag = phase_name.lower().replace('-', '').replace(' ', '_')
            phase_save = f'{base}_{phase_tag}{ext}'
        else:
            phase_save = save_path

        _save_or_display(fig, phase_save, display)


# =========================================================================
#  Scalar time-series plot
# =========================================================================

def plot_scalars(tt, save_path=None, display=True):
    """Plot 4x3 grid of time-series scalars."""
    from scipy.interpolate import interp1d

    s = tt._state
    times = tt._tm_times

    fig, axes = plt.subplots(4, 3, figsize=(16, 12))

    # (0,0): Ip
    ax = axes[0, 0]
    ax.set_title('Ip [A]')
    ax.plot(times, s['Ip_tm'], '-o', markersize=3, label='Ip TM')
    t_Ip, y_Ip = _tx_scalar(tt, 'Ip')
    if t_Ip is not None:
        ax.plot(t_Ip, y_Ip, '-', linewidth=1, label='Ip TX')
    t_Ini, y_Ini = _tx_scalar(tt, 'I_non_inductive')
    if t_Ini is not None:
        ax.plot(t_Ini, y_Ini, '--', linewidth=1, label='Ip NI TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Ip [A]')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # (0,1): psi_lcfs and psi_axis
    ax = axes[0, 1]
    ax.set_title(r'$\psi_{lcfs}$ & $\psi_{axis}$ (TM & TX)')
    ax.plot(times, s['psi_lcfs_tm'], '-', color='tab:blue', label=r'$\psi_{lcfs}$ TM')
    ax.plot(times, s['psi_axis_tm'], '-', color='tab:orange', label=r'$\psi_{axis}$ TM')
    # TORAX psi is COCOS 11; convert to TM-native: psi_TM = -psi_COCOS11 / (2π)
    t_psi_lcfs, y_psi_lcfs = _tx_profile_at_rho(tt, 'psi', 1.0, scale=-1.0/(2.0*np.pi))
    t_psi_axis, y_psi_axis = _tx_profile_at_rho(tt, 'psi', 0.0, scale=-1.0/(2.0*np.pi))
    if t_psi_lcfs is not None:
        ax.plot(t_psi_lcfs, y_psi_lcfs, '--', color='tab:blue', linewidth=1, label=r'$\psi_{lcfs}$ TX')
    if t_psi_axis is not None:
        ax.plot(t_psi_axis, y_psi_axis, '--', color='tab:orange', linewidth=1, label=r'$\psi_{axis}$ TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'$\psi$ [Wb/rad]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (0,2): V_loop
    ax = axes[0, 2]
    ax.set_title('V_loop (TM vs TX) [V]')
    ax.plot(times, s['vloop_tm'], '-o', markersize=3, label='TokaMaker')
    t_vl, y_vl = _tx_scalar(tt, 'v_loop_lcfs')
    if t_vl is not None:
        ax.plot(t_vl, y_vl, '-', linewidth=1, label='TORAX')
    # Compute ratio only at TM timepoints
    tm_vloop = np.array(s['vloop_tm'])
    tx_vloop = np.array(s['vloop_tx'])
    ratio = tm_vloop / np.where(tx_vloop != 0, tx_vloop, np.nan)
    ax2 = ax.twinx()
    ax2.plot(times, ratio, 'g-s', markersize=3, label='TM/TX ratio')
    ax2.set_ylim(0, 2)
    ax2.set_ylabel('V_loop ratio (TM/TX)', color='g')
    ax2.tick_params(axis='y', labelcolor='g')
    ax2.legend(fontsize=8, loc='upper right')
    ax.set_xlabel('Time [s]')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, 3)

    # (1,0): Q_fusion
    ax = axes[1, 0]
    ax.set_title('Q_fusion')
    t_Q, y_Q = _tx_scalar(tt, 'Q_fusion')
    t_E, y_E = _tx_scalar(tt, 'E_fusion')
    if t_Q is not None:
        ax.plot(t_Q, y_Q, '-', linewidth=1, label='Q')
    if t_E is not None:
        ax2 = ax.twinx()
        ax2.plot(t_E, y_E, '-', color='crimson', linewidth=1, label='E_fusion')
        ax2.set_ylabel('E_fusion')
        ax2.legend(fontsize=8, loc='upper right')
    ax.set_xlabel('Time [s]')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper left')

    # (1,1): n_e
    ax = axes[1, 1]
    ax.set_title(r'$n_e$ [m$^{-3}$]')
    t_ne_avg, y_ne_avg = _tx_scalar(tt, 'n_e_line_avg')
    t_ne_core, y_ne_core = _tx_profile_at_rho(tt, 'n_e', 0.0)
    if t_ne_avg is not None:
        ax.plot(t_ne_avg, y_ne_avg, '-', linewidth=1, label='n_e line avg')
    if t_ne_core is not None:
        ax.plot(t_ne_core, y_ne_core, '-', linewidth=1, label='n_e core')
    ne_edge_y = [s['n_e'][ii]['y'][-1] if ii in s.get('n_e', {}) else np.nan for ii in range(len(times))]
    ax.plot(times, ne_edge_y, '-', linewidth=1, label='n_e edge')
    ax.set_xlabel('Time [s]')
    ax2 = ax.twinx()
    ax2.plot(times, s['f_GW'], 'm--', markersize=3, label='f_GW_line')
    ax2.plot(times, s['f_GW_vol'], 'c--', markersize=3, label='f_GW_vol')
    ax2.set_ylabel('f_GW')
    ax2.legend(fontsize=8)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (1,2): T_e
    ax = axes[1, 2]
    ax.set_title(r'$T_e$ [keV]')
    t_Te, y_Te = _tx_profile_at_rho(tt, 'T_e', 0.0)
    t_Ti, y_Ti = _tx_profile_at_rho(tt, 'T_i', 0.0)
    if t_Te is not None:
        ax.plot(t_Te, y_Te, '-', linewidth=1, label='T_e core')
    if t_Ti is not None:
        ax.plot(t_Ti, y_Ti, '-', linewidth=1, label='T_i core')
    te_edge_y = [s['T_e'][ii]['y'][-1] if ii in s.get('T_e', {}) else np.nan for ii in range(len(times))]
    ax.plot(times, te_edge_y, '-', linewidth=1, label='T_e edge')
    ax.set_xlabel('Time [s]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (2,0): Power channels
    ax = axes[2, 0]
    ax.set_title('Power channels [W]')
    for tx_name, label, fmt in [('P_ohmic_e', 'P_ohmic_e', 'r-'), ('P_radiation_e', 'P_radiation_e', 'm--'),
                                 ('P_SOL_total', 'P_SOL_total', 'c--'), ('P_alpha_total', 'P_alpha_total', 'g-.'),
                                 ('P_aux_total', 'P_aux_total', 'y-.')]:
        scale = -1.0 if tx_name == 'P_radiation_e' else 1.0
        tx_t, tx_y = _tx_scalar(tt, tx_name, scale=scale)
        if tx_t is not None:
            ax.plot(tx_t, tx_y, fmt, linewidth=1, label=label)
    ax.set_xlabel('Time [s]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (2,1): beta_N
    ax = axes[2, 1]
    ax.set_title('beta_N')
    t_bN, y_bN = _tx_scalar(tt, 'beta_N')
    if t_bN is not None:
        ax.plot(t_bN, y_bN, '-', linewidth=1, label='beta_N TX')
    ax.plot(times, s['beta_N_tm'], '--o', markersize=3, label='beta_N TM')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('beta_N')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (2,2): l_i
    ax = axes[2, 2]
    ax.set_title('l_i (li3)')
    t_li, y_li = _tx_scalar(tt, 'li3')
    if t_li is not None:
        ax.plot(t_li, y_li, '-', linewidth=1, label='l_i TX')
    ax.plot(times, s['l_i_tm'], '--o', markersize=3, label='l_i TM')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('l_i')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (3,0): q95 and q0
    ax = axes[3, 0]
    ax.set_title('Safety Factor q')
    t_q95, y_q95 = _tx_scalar(tt, 'q95')
    if t_q95 is not None:
        ax.plot(t_q95, y_q95, 'b-', linewidth=1, label='q95 TX')
    t_q0, y_q0 = _tx_profile_at_rho(tt, 'q', 0.0, rho_coord='rho_face_norm')
    if t_q0 is not None:
        ax.plot(t_q0, y_q0, 'r-', linewidth=1, label='q0 TX')
    ax.plot(times, s['q95_tm'], 'b--', label='q95 TM')
    ax.plot(times, s['q0_tm'], 'r--', label='q0 TM')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Safety Factor')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (3,1): pax
    ax = axes[3, 1]
    ax.set_title('pax [Pa]')
    t_pax, y_pax = _tx_profile_at_rho(tt, 'pressure_thermal_total', 0.0)
    if t_pax is not None:
        ax.plot(t_pax, y_pax, '-', linewidth=1, label='pax TX')
    ax.plot(times, s['pax_tm'], '--o', markersize=3, label='pax TM')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('pax [Pa]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (3,2): Flux consumption
    ax = axes[3, 2]
    ax.set_title('Flux Consumption [Wb]')
    psi_lcfs_tm_arr = np.array(s['psi_lcfs_tm'])
    ax.plot(times, -(psi_lcfs_tm_arr - psi_lcfs_tm_arr[0]) * 2 * np.pi, '-o', markersize=3, label='Flux TM')
    if t_psi_lcfs is not None:
        flux_tx = -(y_psi_lcfs - y_psi_lcfs[0]) * 2 * np.pi
        ax.plot(t_psi_lcfs, flux_tx, '-', linewidth=1, label='Flux TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Flux Consumption [Wb]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Scalars', fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    _save_or_display(fig, save_path, display)


# =========================================================================
#  Coil current plot
# =========================================================================

def _coil_net_turns(tt, cname):
    """Return net_turns for a coil set name, defaulting to 1."""
    return tt._tm.coil_sets.get(cname, {}).get('net_turns', 1.0)


def _coil_aturn_clim(tt):
    """Return (min, max) colormap limits in A-turns from _coil_bounds (A/turn * turns)."""
    coil_bounds = getattr(tt, '_coil_bounds', {})
    if not coil_bounds:
        return -1.0, 1.0
    vals = []
    for cname, (lo, hi) in coil_bounds.items():
        n = _coil_net_turns(tt, cname)
        vals.extend([lo * n, hi * n])
    return min(vals), max(vals)


def plot_coils(tt, save_path=None, display=True):
    """Plot coil current traces in MA-turns with limit bands."""
    coil_data = tt._results.get('COIL', {})
    if not coil_data:
        return
    coil_names = sorted(coil_data.keys())
    n_coils = len(coil_names)
    ncols = 3
    nrows = max(1, (n_coils + ncols - 1) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3 * nrows), squeeze=False)
    for k, cname in enumerate(coil_names):
        ax = axes[k // ncols, k % ncols]
        n_turns = _coil_net_turns(tt, cname)
        t_vals = sorted(coil_data[cname].keys())
        i_vals = [coil_data[cname][t_v] * n_turns * 1e-3 for t_v in t_vals]
        ax.plot(t_vals, i_vals, '-o', ms=2, lw=1.2)
        if hasattr(tt, '_coil_bounds') and cname in tt._coil_bounds:
            lo, hi = tt._coil_bounds[cname]
            lo_ka, hi_ka = lo * n_turns * 1e-6, hi * n_turns * 1e-6
            ax.axhline(lo_ka, color='r', ls='--', lw=0.8, label=f'limit {lo_ka:.0f} MA-t')
            ax.axhline(hi_ka, color='r', ls='--', lw=0.8, label=f'limit {hi_ka:.0f} MA-t')
        ax.set_title(cname, fontsize=9)
        ax.set_ylabel('I [MA-turns]', fontsize=8)
        ax.set_xlabel('Time [s]', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=6)
    for k in range(n_coils, nrows * ncols):
        axes[k // ncols, k % ncols].axis('off')
    plt.suptitle(f'Coil Currents (loop {tt._current_loop})', fontsize=12)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    _save_or_display(fig, save_path, display)


# =========================================================================
#  LCFS evolution plot
# =========================================================================

def plot_lcfs_evolution(tt, save_path=None, display=True, one_plot=False):
    """Plot time evolution of the last closed flux surface for each phase.

    Produces phase-split figures by default (rampup, flattop, rampdown),
    or one combined figure when one_plot=True.
    Line color encodes time, with a colorbar on the right.
    Phases not present in the simulation are skipped.
    Uses stored lcfs_geo_tm contours (extracted in _tm_update); falls back to
    trace_surf on the stored equilibrium object if the stored contour is absent.
    """
    s = tt._state
    times = np.array(tt._tm_times)
    lcfs_tm_data = s.get('lcfs_geo_tm', {})
    equil_data = s.get('equil', {})
    if not lcfs_tm_data and not equil_data:
        return

    ft = getattr(tt, '_flattop', np.zeros(len(times), dtype=bool)).astype(bool)

    # Use same phase detection as the rest of toktox
    if np.any(ft):
        ft_indices = np.where(ft)[0]
        ft_start_i = ft_indices[0]
        ft_end_i   = ft_indices[-1]
        rampup_mask   = np.arange(len(times)) < ft_start_i
        flattop_mask  = ft
        rampdown_mask = np.arange(len(times)) > ft_end_i
    else:
        rampup_mask   = np.ones(len(times), dtype=bool)
        flattop_mask  = np.zeros(len(times), dtype=bool)
        rampdown_mask = np.zeros(len(times), dtype=bool)

    phases = [
        ('Ramp-up',   rampup_mask),
        ('Flattop',   flattop_mask),
        ('Ramp-down', rampdown_mask),
    ]
    if one_plot:
        phases = [('Pulse', np.ones(len(times), dtype=bool))]

    # Limiter contours from TokaMaker (col 0 = R, col 1 = Z; swap if mirror_mode)
    lim_contours = getattr(tt._tm, 'lim_contours', None) or []
    mirror_mode  = getattr(getattr(tt._tm, 'settings', None), 'mirror_mode', False)

    cmap = cm.plasma

    for phase_name, mask in phases:
        indices = np.where(mask)[0]
        indices = [i for i in indices if lcfs_tm_data.get(i) is not None or equil_data.get(i) is not None]
        if not indices:
            continue

        phase_times = times[indices]
        t_min, t_max = phase_times[0], phase_times[-1]
        norm = Normalize(vmin=t_min, vmax=t_max if t_max > t_min else t_min + 1e-9)

        fig, ax = plt.subplots(figsize=(6, 8))
        if one_plot:
            fig.suptitle(f'LCFS Evolution (loop {tt._current_loop})', fontsize=TITLE_FS)
        else:
            fig.suptitle(f'LCFS Evolution — {phase_name} (loop {tt._current_loop})', fontsize=TITLE_FS)

        # Draw limiter
        for lc in lim_contours:
            lc = np.asarray(lc)
            if mirror_mode:
                ax.plot(lc[:, 1], lc[:, 0], color='k', linewidth=1.5, zorder=5)
            else:
                ax.plot(lc[:, 0], lc[:, 1], color='k', linewidth=1.5, zorder=5)

        for i in indices:
            # Prefer pre-extracted contour; fall back to trace_surf if missing
            lcfs = lcfs_tm_data.get(i)
            if lcfs is None:
                equil = equil_data.get(i)
                if equil is None:
                    continue
                try:
                    lcfs = equil.trace_surf(1.0)
                except Exception:
                    try:
                        lcfs = equil.trace_surf(0.99)
                    except Exception:
                        print('Failed to trace LCFS from equilibrium at time index', i)
                        continue
            if lcfs is None or len(lcfs) == 0:
                continue
            lcfs = np.asarray(lcfs)
            color = cmap(norm(times[i]))
            ax.plot(lcfs[:, 0], lcfs[:, 1], color=color, linewidth=1.0, alpha=0.85)

        ax.set_xlabel('R [m]', fontsize=LABEL_FS)
        ax.set_ylabel('Z [m]', fontsize=LABEL_FS)
        ax.set_aspect('equal')
        _style(ax)

        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label('Time [s]', fontsize=LABEL_FS)
        cbar.ax.tick_params(labelsize=TICK_FS)

        plt.tight_layout()

        if save_path is not None and not one_plot:
            base, ext = os.path.splitext(save_path)
            phase_tag = phase_name.lower().replace('-', '').replace(' ', '_')
            phase_save = f'{base}_{phase_tag}{ext}'
        else:
            phase_save = save_path

        _save_or_display(fig, phase_save, display)


# =========================================================================
#  Movie generation
# =========================================================================

def make_movie(tt, save_path=None, display=True, speed_factor=1.0, loop=None, notebook_mode=None):
    """Create pulse movie from simulation data.

    Renders equilibrium plots from stored psi snapshots (tt._tm_psi_on_nodes),
    generates composite frames in a temp directory, encodes to MP4, then cleans up.

    Parameters
    ----------
    tt : TokTox
        Fully-populated TokTox object (after fly()).
    save_path : str, optional
        Path to save the MP4 file. If None, uses a default in cwd.
    display : bool
        If True and running in Jupyter, embed the video in the notebook.
    speed_factor : float
        Playback speed relative to real time.
    loop : int, optional
        Loop to visualize. Default: last completed loop.
    notebook_mode : bool or None
        True  → embed video in notebook after saving.
        False → save to file only, do not embed.
        None  → auto-detect (embed if in Jupyter and display=True).
    """
    if loop is None:
        loop = tt._current_loop

    tmp_dir = _make_temp_dir()
    try:
        times = tt._tm_times
        n = len(times)

        psi_lcfs_tm = np.array(tt._state['psi_lcfs_tm'])
        psi_lcfs_tx = np.array(tt._state['psi_lcfs_tx'])
        flux_con_tm = -(psi_lcfs_tm - psi_lcfs_tm[0]) * 2.0 * np.pi
        flux_con_tx = -(psi_lcfs_tx - psi_lcfs_tx[0]) * 2.0 * np.pi

        # Render equil plots from stored psi snapshots
        equil_dir = os.path.join(tmp_dir, 'equil')
        os.makedirs(equil_dir, exist_ok=True)
        _render_equil_frames(tt, loop, equil_dir)

        # Render composite frames
        for idx in range(n):
            fpath = os.path.join(tmp_dir, f'frame_{idx:04d}.png')
            _render_frame(tt, loop, idx, times[idx], times,
                          flux_con_tm, flux_con_tx, fpath, tt._run_name, equil_dir)
            plt.close('all')

        # Encode video
        total_time = times[-1] - times[0] if len(times) > 1 else 1.0
        fps = max(1, speed_factor * n / total_time)

        if save_path is None:
            save_path = os.path.join(os.getcwd(), f'toktox_pulse_loop{loop:03d}.mp4')

        _encode_video_from_dir(tmp_dir, save_path, fps=fps)

        show_in_nb = notebook_mode if notebook_mode is not None else (display and _in_jupyter())
        if show_in_nb:
            _embed_video_jupyter(save_path)

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _render_equil_frames(tt, loop, equil_dir):
    """Render equilibrium plots from stored equilibrium snapshots to PNG files."""
    equil_data = tt._state.get('equil', {})

    for i in range(len(tt._tm_times)):
        out_path = os.path.join(equil_dir, f'equil_{loop:03d}.{i:03d}.png')
        equil = equil_data.get(i)
        if equil is None:
            # Create a "failed" placeholder
            fig, ax = plt.subplots(1, 1, figsize=(11, 12))
            ax.text(0.5, 0.5, 'TokaMaker failed\nto converge',
                    transform=ax.transAxes, fontsize=16,
                    ha='center', va='center', color='darkred', fontweight='bold')
            ax.axis('off')
            fig.savefig(out_path, dpi=150, bbox_inches='tight', pad_inches=0.0)
            plt.close(fig)
            continue

        fig, ax = plt.subplots(1, 1, figsize=(11, 12))
        min_bound, max_bound = _coil_aturn_clim(tt)
        cb = tt._tm.plot_machine(fig, ax, equilibrium=equil, coil_colormap='seismic', coil_symmap=False,
                                  coil_scale=1.E-6, coil_clabel=r'$I_C$ [MA-turns]')
        tt._tm.plot_constraints(fig, ax, equilibrium=equil)
        if cb is not None:
            cb.mappable.set_clim(min_bound * 1e-6, max_bound * 1e-6)
        tt._tm.plot_psi(fig, ax, equilibrium=equil, xpoint_color='r', vacuum_nlevels=4)
        x_pt = getattr(tt, '_x_point_targets', None)
        if x_pt is not None and _x_points_active(tt, i, t=tt._tm_times[i]):
            ax.plot(x_pt[:, 0], x_pt[:, 1], 'x', color='purple', markersize=10, markeredgewidth=2, label='Saddle point targets')
        sp = tt._state.get('strike_pts', {}).get(i)
        if sp is not None and len(sp) > 0:
            ax.plot(sp[:, 0], sp[:, 1], 'g^', markersize=10, markeredgewidth=2, label='Strike point targets')
        ax.set_aspect('equal')
        ax.legend(loc='upper right', fontsize=12)
        ax.tick_params(labelsize=11)
        if cb is not None:
            cb.ax.tick_params(labelsize=11)
        fig.savefig(out_path, dpi=150, bbox_inches='tight', pad_inches=0.0)
        plt.close(fig)


def _render_frame(tt, loop, idx, t_now, times, flux_con_tm, flux_con_tx, out_path, run_name, equil_dir):
    """Create a single composite movie frame."""
    from matplotlib.image import imread

    fig = plt.figure(figsize=(MOVIE_FIG_W, MOVIE_FIG_H), dpi=MOVIE_DPI)
    gs = GridSpec(6, 3, figure=fig, width_ratios=[1.2, 1.0, 1.0],
                  wspace=0.28, hspace=0.18,
                  left=0.045, right=0.94, top=0.96, bottom=0.05)

    # Column 0: info + equil
    ax_text = fig.add_subplot(gs[0:1, 0])
    ax_equil = fig.add_subplot(gs[1:5, 0])
    _draw_info(ax_text, tt, loop, run_name)
    _draw_equil_from_dir(ax_equil, tt, loop, idx, t_now, equil_dir)

    # Column 1: scalar time-series
    sax = [fig.add_subplot(gs[j, 1]) for j in range(6)]
    _draw_scalars_movie(sax, tt, times, t_now, flux_con_tm, flux_con_tx)
    for ax in sax[:-1]:
        ax.tick_params(labelbottom=False)
    sax[-1].set_xlabel('Time [s]', fontsize=LABEL_FS)

    # Column 2: profile snapshots
    pax = [fig.add_subplot(gs[j, 2]) for j in range(6)]
    _draw_profiles_movie(pax, tt, idx)
    for ax in pax[:-1]:
        ax.tick_params(labelbottom=False)
    pax[-1].set_xlabel(r'$\psi_N$', fontsize=LABEL_FS)

    fig.savefig(out_path, dpi=MOVIE_DPI)
    plt.close(fig)


def _draw_info(ax, tt, loop, run_name):
    ax.axis('off')
    t_ave = getattr(tt, '_t_ave_toggle', 'off')
    dt_val = getattr(tt, '_tx_dt', getattr(tt, '_dt', None))
    lines = [
        f'Run:   {run_name}            loop:  {loop}',
        f'n_rho: {tt._n_rho}          dt:    {dt_val} s',
        f'time range:     [{tt._t_init}, {tt._t_final}] s          times: {len(tt._tm_times)}',
        f'LSF:   {tt._last_surface_factor}',
        f't_ave: {t_ave}    window: {getattr(tt, "_t_ave_window", 0):.2f} s',
        f'causal: {getattr(tt, "_t_ave_causal", True)}    '
        f'ignore_start: {getattr(tt, "_t_ave_ignore_start", 0):.2f} s',
    ]
    ax.text(0.05, 0.95, '\n'.join(lines), transform=ax.transAxes, fontsize=INFO_FS,
            va='top', fontfamily='monospace')


def _draw_equil_from_dir(ax, tt, loop, idx, t_now, equil_dir):
    """Draw equilibrium from pre-rendered PNG in equil_dir."""
    from matplotlib.image import imread

    s = tt._state
    eq_img = os.path.join(equil_dir, f'equil_{loop:03d}.{idx:03d}.png')
    tm_ok = os.path.exists(eq_img)

    if tm_ok:
        img = imread(eq_img)
        ax.imshow(img, aspect='equal')
        ax.axis('off')
        ax.set_title(f't = {t_now:.2f} s', fontsize=TITLE_FS + 2)

        div_flags = getattr(tt, '_diverted_flags', {})
        div_flag = div_flags.get(idx)
        config_str = 'Diverted' if div_flag else ('Limited' if div_flag is not None else '?')

        diag = (
            f"i = {idx}    "
            f"Ip = {abs(s['Ip_tm'][idx]) / 1e6:.3f} MA    "
            f"pax = {s['pax_tm'][idx] / 1e3:.1f} kPa\n"
            f"R = {s['R0_mag'][idx]:.3f} m    "
            f"a = {s['a'][idx]:.3f} m    "
            f"B0 = {s['B0'][idx]:.3f} T\n"
            f"\u03ba = {s['kappa'][idx]:.3f}    "
            f"\u03b4 = {s['delta'][idx]:.3f}    "
            f"Configuration = {config_str}\n"
        )
        ax.text(0.5, -0.02, diag, transform=ax.transAxes, fontsize=DIAG_FS,
                ha='center', va='top', fontfamily='monospace')
    else:
        ax.axis('off')
        ax.text(0.5, 0.55, 'TokaMaker failed\nto converge',
                transform=ax.transAxes, fontsize=16,
                ha='center', va='center', color='darkred', fontweight='bold')
        ax.set_title(f't = {t_now:.2f} s  (FAILED)', fontsize=TITLE_FS + 2, color='darkred')


def _draw_scalars_movie(axes, tt, times, t_now, flux_con_tm, flux_con_tx):
    """Draw the 6 scalar time-series panels for a movie frame."""
    s = tt._state

    # 1: Ip
    ax = axes[0]
    ax.plot(times, np.array(s['Ip_tm']) / 1e6, color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label='Ip TM')
    ax.plot(times, np.array(s['Ip_tx']) / 1e6, color=COLOR_TX, ls=LS_PRI, lw=LW, label='Ip TX')
    ax.plot(times, np.array(s['Ip_ni_tx']) / 1e6, color=COLOR_TX, ls=LS_SEC, lw=LW, label='Ip_ni TX')
    ax.set_ylabel('Ip [MA]', fontsize=LABEL_FS)
    ax.set_title('Scalars', fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_FS, loc='upper left')
    _style(ax)

    # 2: V_loop and l_i
    ax = axes[1]
    ax.plot(times, s['vloop_tx'], color=COLOR_TX, ls=LS_PRI, lw=LW, label='Vloop TX')
    ax.set_ylabel('V_loop [V]', fontsize=LABEL_FS)
    _legend_if_labeled(ax, fontsize=LEGEND_FS, loc='upper left')
    _style(ax)
    ax2 = ax.twinx()
    if s.get('l_i_tm', None) is not None:
        ax2.plot(times, s['l_i_tm'], color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label='l_i TM')
    if s.get('l_i_tx', None) is not None:
        ax2.plot(times, s['l_i_tx'], color=COLOR_TX, ls=LS_SEC, lw=LW, label='l_i TX')
    ax2.set_ylabel('l_i', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    _legend_if_labeled(ax2, fontsize=LEGEND_FS, loc='upper right')

    # 3: Power channels and Q
    ax = axes[2]
    pkeys = [
        ('P_ohmic_e', 'Ohmic', COLORS_MULTI[0]),
        ('P_aux_total', 'Aux', COLORS_MULTI[1]),
        ('P_alpha_total', 'Fusion', COLORS_MULTI[2]),
        ('P_radiation_e', 'Radiation', COLORS_MULTI[3]),
        ('P_SOL_total', 'SOL', COLORS_MULTI[4]),
    ]
    for tx_name, label, clr in pkeys:
        scale = -1.0 if tx_name == 'P_radiation_e' else 1.0
        tx_t, tx_y = _tx_scalar(tt, tx_name, scale=scale)
        if tx_t is not None:
            ax.plot(tx_t, tx_y / 1e6, color=clr, ls=LS_PRI, lw=LW, label=label)
    ax.set_ylabel('Power [MW]', fontsize=LABEL_FS)
    ax.legend(fontsize=LEGEND_FS, loc='upper left', ncol=2)
    _style(ax)
    ax2 = ax.twinx()
    t_Q, y_Q = _tx_scalar(tt, 'Q_fusion')
    if t_Q is not None:
        ax2.plot(t_Q, y_Q, color='indigo', ls=LS_SEC, lw=LW, label='Q')
    ax2.set_ylabel('Q', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    _legend_if_labeled(ax2, fontsize=LEGEND_FS, loc='upper right')

    # 4: psi_axis and psi_lcfs + flux consumed
    ax = axes[3]
    ax.plot(times, s['psi_axis_tm'], color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label='\u03c8_axis TM')
    ax.plot(times, s['psi_axis_tx'], color=COLOR_TX, ls=LS_PRI, lw=LW, label='\u03c8_axis TX')
    ax.set_ylabel('\u03c8 [Wb/rad]', fontsize=LABEL_FS)
    ax.plot(times, s['psi_lcfs_tm'], color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label='\u03c8_lcfs TM')
    ax.plot(times, s['psi_lcfs_tx'], color=COLOR_TX, ls=LS_SEC, lw=LW, label='\u03c8_lcfs TX')
    ax.tick_params(labelsize=TICK_FS)
    _style(ax)
    _legend_if_labeled(ax, fontsize=LEGEND_FS, loc='upper left')
    ax2 = ax.twinx()
    ax2.plot(times, flux_con_tm, color='darkorange', ls='-.', lw=LW, marker=MK_TM, ms=MK_SZ, label='Flux TM')
    ax2.plot(times, flux_con_tx, color='seagreen', ls='-.', lw=LW, label='Flux TX')
    ax2.set_ylabel('Flux Consumed [Wb]', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    _legend_if_labeled(ax2, fontsize=LEGEND_FS, loc='upper right')

    coil_data = tt._results.get('COIL', {})
    cs_coils = []
    other_coils = []
    for cname, cvals in sorted(coil_data.items()):
        cname_tokens = [tok for tok in cname.upper().replace('-', '_').split('_') if tok]
        if any(tok.startswith('CS') for tok in cname_tokens):
            cs_coils.append((cname, cvals))
        else:
            other_coils.append((cname, cvals))

    # 5: CS coil currents
    ax = axes[4]
    if cs_coils:
        cs_colors = plt.cm.tab10(np.linspace(0, 1, max(len(cs_coils), 1)))
        for ci, (cname, cvals) in enumerate(cs_coils):
            ct = sorted(cvals.keys())
            n_turns = _coil_net_turns(tt, cname)
            ci_vals = [cvals[t_v] * n_turns * 1e-6 for t_v in ct]
            ax.plot(ct, ci_vals, ls=LS_PRI, lw=LW * 0.8, color=cs_colors[ci], label=cname)
        _legend_if_labeled(ax, fontsize=LEGEND_FS - 2, loc='upper left', ncol=2)
    else:
        ax.text(0.5, 0.5, 'No CS coils', transform=ax.transAxes,
                ha='center', va='center', fontsize=LABEL_FS)
    ax.set_ylabel('I_coil [MA-turns]', fontsize=LABEL_FS)
    _style(ax)

    # 6: PF and other coil currents
    ax = axes[5]
    if other_coils:
        oth_colors = plt.cm.tab20(np.linspace(0, 1, max(len(other_coils), 1)))
        for ci, (cname, cvals) in enumerate(other_coils):
            ct = sorted(cvals.keys())
            n_turns = _coil_net_turns(tt, cname)
            ci_vals = [cvals[t_v] * n_turns * 1e-6 for t_v in ct]
            ax.plot(ct, ci_vals, ls=LS_PRI, lw=LW * 0.75, color=oth_colors[ci], label=cname)
        _legend_if_labeled(ax, fontsize=LEGEND_FS - 3, loc='upper left', ncol=2)
    else:
        ax.text(0.5, 0.5, 'No PF/other coils', transform=ax.transAxes,
                ha='center', va='center', fontsize=LABEL_FS)
    ax.set_ylabel('I_coil [MA-turns]', fontsize=LABEL_FS)
    _style(ax)

    _vline(axes, t_now)


def _draw_profiles_movie(axes, tt, idx):
    """Draw the 6 radial profile panels for a movie frame."""
    s = tt._state

    # 1: q profile
    ax = axes[0]
    x, y = _prof(s['q_prof_tm'], idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label='q TM')
    x, y = _prof(s['q_prof_tx'], idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TX, ls=LS_PRI, lw=LW, label='q TX')
    ax.set_ylabel('q', fontsize=LABEL_FS)
    ax.set_title('Profiles', fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_FS, loc='best')
    _style(ax)

    # 2: n_e / T_e
    ax = axes[1]
    x, y = _prof(s.get('n_e', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TX, ls=LS_PRI, lw=LW, label='ne TX')
    ax.set_ylabel('ne [m\u207b\u00b3]', fontsize=LABEL_FS)
    ax.legend(fontsize=LEGEND_FS, loc='lower left')
    _style(ax)
    x, y = _prof(s.get('T_e', {}), idx)
    if x is not None:
        ax2 = ax.twinx()
        ax2.plot(x, y, color=COLOR_TX, ls=LS_SEC, lw=LW, label='Te TX')
        x_ti, y_ti = _prof(s.get('T_i', {}), idx)
        if x_ti is not None:
            ax2.plot(x_ti, y_ti, color='forestgreen', ls='--', lw=LW, label='Ti TX')
        ax2.set_ylabel('T [keV]', fontsize=LABEL_FS)
        ax2.tick_params(labelsize=TICK_FS)
        ax2.legend(fontsize=LEGEND_FS, loc='upper right')

    # 3: j components / FF'
    ax = axes[2]
    j_keys = [
        ('j_tot', 'j_tot', COLORS_MULTI[0]),
        ('j_ohmic', 'j_ohm', COLORS_MULTI[1]),
        ('j_ni', 'j_ni', COLORS_MULTI[2]),
        ('j_bootstrap', 'j_BS', COLORS_MULTI[3]),
        ('j_ecrh', 'j_EC', COLORS_MULTI[4]),
        ('j_generic_current', 'j_gen', COLORS_MULTI[5]),
    ]
    for skey, label, clr in j_keys:
        x, y = _prof(s.get(skey, {}), idx)
        if x is not None:
            ax.plot(x, y / 1e6, color=clr, ls=LS_PRI, lw=LW, label=label)
    ax.set_ylabel('j [MA/m\u00b2]', fontsize=LABEL_FS)
    ax.legend(fontsize=LEGEND_FS - 1, loc='upper left', ncol=2)
    _style(ax)
    ax2 = ax.twinx()
    x, y = _prof(s.get('ffp_prof_tm', {}), idx)
    if x is not None:
        ax2.plot(x, y, color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label="FF' TM")
    x, y = _prof(s.get('ffp_prof_tx', {}), idx)
    if x is not None:
        ax2.plot(x, y, color=COLOR_TX, ls=LS_SEC, lw=LW, label="FF' TX")
    ax2.set_ylabel("FF'", fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    ax2.legend(fontsize=LEGEND_FS, loc='upper right')

    # 4: p' / p
    ax = axes[3]
    x, y = _prof(s.get('pp_prof_tm', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label="p' TM")
    x, y = _prof(s.get('pp_prof_tx', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TX, ls=LS_PRI, lw=LW, label="p' TX")
    ax.set_ylabel("p'", fontsize=LABEL_FS)
    ax.legend(fontsize=LEGEND_FS, loc='center')
    _style(ax)
    ax2 = ax.twinx()
    x, y = _prof(s.get('p_prof_tm', {}), idx)
    if x is not None:
        ax2.plot(x, y, color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label='p TM')
    x, y = _prof(s.get('p_prof_tx', {}), idx)
    if x is not None:
        ax2.plot(x, y, color=COLOR_TX, ls=LS_SEC, lw=LW, label='p TX')
    ax2.set_ylabel('p [Pa]', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    ax2.legend(fontsize=LEGEND_FS, loc='upper right')

    # 5: <1/R> / psi profile
    ax = axes[4]
    x, y = _prof(s.get('R_inv_avg_tm', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label='<1/R> TM')
    x, y = _prof(s.get('R_inv_avg_tx', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TX, ls=LS_PRI, lw=LW, label='<1/R> TX')
    ax.set_ylabel('<1/R> [1/m]', fontsize=LABEL_FS)
    ax.legend(fontsize=LEGEND_FS, loc='center left')
    _style(ax)
    ax2 = ax.twinx()
    x, y = _prof(s.get('psi_tm', {}), idx)
    if x is not None:
        ax2.plot(x, y, color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label='\u03c8 TM')
    x, y = _prof(s.get('psi_tx', {}), idx)
    if x is not None:
        ax2.plot(x, y, color=COLOR_TX, ls=LS_SEC, lw=LW, label='\u03c8 TX')
    ax2.set_ylabel('\u03c8 [Wb/rad]', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    ax2.legend(fontsize=LEGEND_FS, loc='center right')

    # 6: resistivity
    ax = axes[5]
    x, y = _prof(s.get('eta_prof', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TX, ls=LS_PRI, lw=LW, label='\u03b7 TX')
    ax.set_ylabel('\u03b7 [\u03a9\u00b7m]', fontsize=LABEL_FS)
    ax.set_yscale('log')
    ax.legend(fontsize=LEGEND_FS)
    _style(ax)


def _encode_video_from_dir(frame_dir, out_path, fps=2):
    """Encode frames from a directory into an MP4 file."""
    pattern = os.path.join(frame_dir, 'frame_%04d.png')
    try:
        subprocess.run(
            ['ffmpeg', '-y', '-framerate', str(int(fps)),
             '-i', pattern,
             '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '18',
             out_path],
            check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f'Warning: ffmpeg failed to encode MP4: {e}')


def _embed_video_jupyter(video_path):
    """Embed an MP4 video in a Jupyter notebook."""
    try:
        from IPython.display import Video, display
        display(Video(video_path, embed=True, mimetype='video/mp4'))
    except ImportError:
        print(f'Video saved to: {video_path}')


# =========================================================================
#  Interactive widgets (Jupyter)
# =========================================================================

def plot_profiles_interactive(tt):
    """Interactive profile viewer with ipywidgets slider to scrub through timesteps."""
    try:
        import ipywidgets as widgets
        from IPython.display import display, clear_output
    except ImportError:
        print("ipywidgets not installed. Install with: pip install ipywidgets")
        return

    s = tt._state
    times = tt._tm_times
    out = widgets.Output()

    def _update(i):
        with out:
            clear_output(wait=True)
            t = times[i]
            fig, axes = plt.subplots(2, 3, figsize=(16, 8))
            fig.suptitle(f't = {t:.3f} s  (index {i}/{len(times) - 1})', fontsize=14)

            # n_e
            ax = axes[0, 0]
            if i in s.get('n_e', {}):
                ax.plot(s['n_e'][i]['x'], s['n_e'][i]['y'], 'b-')
            ax.set_title('n_e')
            ax.set_ylabel(r'$n_e$ [m$^{-3}$]')

            # T_e
            ax = axes[0, 1]
            if i in s.get('T_e', {}):
                ax.plot(s['T_e'][i]['x'], s['T_e'][i]['y'], 'r-')
            ax.set_title('T_e')
            ax.set_ylabel('T_e [keV]')

            # q
            ax = axes[0, 2]
            if i in s.get('q_prof_tm', {}):
                ax.plot(s['q_prof_tm'][i]['x'], s['q_prof_tm'][i]['y'], color=COLOR_TM, label='TM')
            if i in s.get('q_prof_tx', {}):
                ax.plot(s['q_prof_tx'][i]['x'], s['q_prof_tx'][i]['y'], color=COLOR_TX, label='TX')
            ax.set_title('q')
            ax.set_ylabel('q')
            ax.legend(fontsize=8)

            # j components
            ax = axes[1, 0]
            for key, label, clr in [('j_tot', 'j_tot', 'k'), ('j_ohmic', 'j_ohm', 'r'),
                                     ('j_ni', 'j_ni', 'b'), ('j_bootstrap', 'j_BS', 'g')]:
                if i in s.get(key, {}):
                    ax.plot(s[key][i]['x'], s[key][i]['y'] / 1e6, color=clr, label=label)
            ax.set_title('Current densities')
            ax.set_ylabel('j [MA/m²]')
            ax.legend(fontsize=7)

            # p / p'
            ax = axes[1, 1]
            if i in s.get('ptot', {}):
                ax.plot(s['ptot'][i]['x'], s['ptot'][i]['y'], color=COLOR_TX, label='p TX')
            if i in s.get('p_prof_tm', {}):
                ax.plot(s['p_prof_tm'][i]['x'], s['p_prof_tm'][i]['y'], color=COLOR_TM, ls='--', label='p TM')
            ax.set_title('Pressure')
            ax.set_ylabel('p [Pa]')
            ax.legend(fontsize=8)

            # T_i
            ax = axes[1, 2]
            if i in s.get('T_i', {}):
                ax.plot(s['T_i'][i]['x'], s['T_i'][i]['y'], 'm-')
            ax.set_title('T_i')
            ax.set_ylabel('T_i [keV]')

            for ax in axes.flat:
                ax.set_xlabel(r'$\hat{\psi}$')
                ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.show()

    slider = widgets.IntSlider(value=0, min=0, max=len(times) - 1, loop=1,
                                description='Time idx:', continuous_update=False)
    slider.observe(lambda change: _update(change['new']), names='value')
    display(slider, out)
    _update(0)


def plot_equil_interactive(tt, loop=None, notebook_mode=None, save_path=None):
    """Equilibrium viewer — widget slider in notebook, saved PNGs otherwise.

    Parameters
    ----------
    notebook_mode : bool or None
        True  → ipywidgets slider (Jupyter).
        False → save PNG files to save_path.
        None  → auto-detect (widget if in Jupyter, save files otherwise).
    save_path : str, optional
        Directory for PNG files when notebook_mode=False.
        Defaults to ./equil_loop<N>/ in the current directory.
    """
    if notebook_mode is None:
        notebook_mode = _in_jupyter()

    if loop is None:
        loop = tt._current_loop

    equil_data = tt._state.get('equil', {})
    times = tt._tm_times
    coil_bounds = getattr(tt, '_coil_bounds', {})

    def _draw_equil_ax(fig, ax, i):
        t = times[i]
        equil = equil_data.get(i)
        if equil is not None:
            min_bound = min(b for bounds in coil_bounds.values() for b in bounds) * 1.E-6 if coil_bounds else -1
            max_bound = max(b for bounds in coil_bounds.values() for b in bounds) * 1.E-6 if coil_bounds else 1
            cb = tt._tm.plot_machine(fig, ax, equilibrium=equil, coil_colormap='seismic', coil_symmap=False,
                                      coil_scale=1.E-6, coil_clabel=r'$I_C$ [MA-turns]')
            if cb is not None:
                cb.mappable.set_clim(min_bound, max_bound)
            tt._tm.plot_constraints(fig, ax, equilibrium=equil)
            tt._tm.plot_psi(fig, ax, equilibrium=equil, xpoint_color='r', vacuum_nlevels=4)
            x_pt = getattr(tt, '_x_point_targets', None)
            if x_pt is not None and _x_points_active(tt, i, t=t):
                ax.plot(x_pt[:, 0], x_pt[:, 1], 'rx', markersize=10, markeredgewidth=2,
                        label='Saddle point targets')
            sp = tt._state.get('strike_pts', {}).get(i)
            if sp is not None and len(sp) > 0:
                ax.plot(sp[:, 0], sp[:, 1], 'g^', markersize=10, markeredgewidth=2, label='Strike points')
            ax.set_aspect('equal')
            ax.set_title(f't = {t:.3f} s  (index {i})', fontsize=14)
        else:
            ax.text(0.5, 0.5, 'TokaMaker failed to converge', transform=ax.transAxes,
                    fontsize=14, ha='center', va='center', color='darkred')
            ax.set_title(f't = {t:.3f} s  (FAILED)', fontsize=14, color='darkred')

    if notebook_mode:
        try:
            import ipywidgets as widgets
            from IPython.display import display, clear_output
        except ImportError:
            print("ipywidgets not installed. Install with: pip install ipywidgets")
            return

        out = widgets.Output()

        def _update(i):
            with out:
                clear_output(wait=True)
                fig, ax = plt.subplots(1, 1, figsize=(8, 9))
                _draw_equil_ax(fig, ax, i)
                plt.tight_layout()
                display(fig)
                plt.close(fig)

        slider = widgets.IntSlider(value=0, min=0, max=len(times) - 1, loop=1,
                                    description='Time idx:', continuous_update=False)
        slider.observe(lambda change: _update(change['new']), names='value')
        display(slider, out)
        _update(0)
    else:
        if save_path is None:
            save_path = os.path.join(os.getcwd(), f'equil_loop{loop:03d}')
        os.makedirs(save_path, exist_ok=True)
        for i in range(len(times)):
            fig, ax = plt.subplots(1, 1, figsize=(8, 9))
            _draw_equil_ax(fig, ax, i)
            plt.tight_layout()
            out_path = os.path.join(save_path, f'equil_{i:03d}.png')
            fig.savefig(out_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
        print(f'Saved {len(times)} equilibrium plots to {save_path}')


# =========================================================================
#  Physics summary
# =========================================================================

def summary(tt):
    """Print and return key physics figures of merit from the simulation.

    Returns a dict with summary quantities.
    """
    s = tt._state
    times = np.array(tt._tm_times)

    # Flattop mask
    ft = getattr(tt, '_flattop', np.zeros(len(times), dtype=bool))
    ft_mask = ft.astype(bool)

    out = {}

    # Fusion (from data_tree at full resolution)
    t_Q, y_Q = _tx_scalar(tt, 'Q_fusion')
    if t_Q is not None:
        out['Q_max'] = float(np.nanmax(y_Q))
        out['Q_max_time'] = float(t_Q[np.nanargmax(y_Q)])
        if np.any(ft_mask):
            ft_start, ft_end = times[ft_mask][0], times[ft_mask][-1]
            ft_q_mask = (t_Q >= ft_start) & (t_Q <= ft_end)
            out['Q_flattop_avg'] = float(np.nanmean(y_Q[ft_q_mask])) if np.any(ft_q_mask) else None
        else:
            out['Q_flattop_avg'] = None

    _, y_E = _tx_scalar(tt, 'E_fusion')
    if y_E is not None:
        out['E_fusion_total_MJ'] = float(y_E[-1] / 1e6)

    # Ip
    out['Ip_max_MA'] = float(np.max(np.abs(s['Ip_tm'])) / 1e6)

    # Beta
    _, y_bN = _tx_scalar(tt, 'beta_N')
    if y_bN is not None:
        out['beta_N_max'] = float(np.nanmax(y_bN))
    out['beta_N_tm_max'] = float(np.nanmax(s['beta_N_tm']))

    # H98
    t_H98, y_H98 = _tx_scalar(tt, 'H98')
    if y_H98 is not None:
        out['H98_max'] = float(np.nanmax(y_H98))
        if np.any(ft_mask):
            ft_start, ft_end = times[ft_mask][0], times[ft_mask][-1]
            ft_h_mask = (t_H98 >= ft_start) & (t_H98 <= ft_end)
            out['H98_flattop_avg'] = float(np.nanmean(y_H98[ft_h_mask])) if np.any(ft_h_mask) else None

    # Temperatures
    _, y_Te = _tx_profile_at_rho(tt, 'T_e', 0.0)
    _, y_Ti = _tx_profile_at_rho(tt, 'T_i', 0.0)
    if y_Te is not None:
        out['T_e_core_max_keV'] = float(np.nanmax(y_Te))
    if y_Ti is not None:
        out['T_i_core_max_keV'] = float(np.nanmax(y_Ti))

    # Density
    _, y_ne = _tx_scalar(tt, 'n_e_line_avg')
    if y_ne is not None:
        out['n_e_line_avg_max'] = float(np.nanmax(y_ne))

    # Greenwald fraction
    out['f_GW_max'] = float(np.nanmax(s['f_GW']))

    # Safety factor
    out['q95_min'] = float(np.nanmin(s['q95_tm'][s['q95_tm'] > 0])) if np.any(s['q95_tm'] > 0) else None
    out['q0_min'] = float(np.nanmin(s['q0_tm'][s['q0_tm'] > 0])) if np.any(s['q0_tm'] > 0) else None

    # Flux consumption
    psi_lcfs_tm = np.array(s['psi_lcfs_tm'])
    out['flux_consumed_Wb'] = float(-(psi_lcfs_tm[-1] - psi_lcfs_tm[0]) * 2 * np.pi)

    # Power
    _, y_Pa = _tx_scalar(tt, 'P_alpha_total')
    _, y_Po = _tx_scalar(tt, 'P_ohmic_e')
    if y_Pa is not None:
        out['P_fusion_max_MW'] = float(np.nanmax(y_Pa) / 1e6 * 5)
    if y_Po is not None:
        out['P_ohmic_max_MW'] = float(np.nanmax(y_Po) / 1e6)

    # Internal inductance
    out['l_i_flattop_avg'] = float(np.nanmean(s['l_i_tm'][ft_mask])) if np.any(ft_mask) else None

    # Loop voltage
    if np.any(ft_mask):
        out['vloop_tm_flattop_avg_V'] = float(np.nanmean(np.array(s['vloop_tm'])[ft_mask]))
        out['vloop_tx_flattop_avg_V'] = float(np.nanmean(np.array(s['vloop_tx'])[ft_mask]))

    # Print formatted summary
    print(f"\n{'=' * 55}")
    print(f"  TokTox Physics Summary")
    print(f"{'=' * 55}")
    for key, val in out.items():
        if val is not None:
            if isinstance(val, float):
                print(f"  {key:30s}  {val:.4g}")
            else:
                print(f"  {key:30s}  {val}")
    print(f"{'=' * 55}\n")

    return out

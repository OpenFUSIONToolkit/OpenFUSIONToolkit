#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
r'''! TokaMaker + TORAX Coupled Pulse Design and Simulation Workflow (TokaMaker_TORAX/tmtx)

    TokaMaker (tm) = free-boundary equilibrium solver, part of OFT
    TORAX (tx) = Python Jax transport solver (developed by Google DeepMind)
        https://github.com/google-deepmind/torax

    TokaMaker_TORAX workflow couples the two codes for pulse planning, predictive kinetic equilibria, 
        and other integrated modeling applications.
    
    @authors Freddie Sheehan and John Lhota
    @date May 2026
    @ingroup doxy_oft_python
'''
import copy
import json
import logging
import os
import platform
from pathlib import Path
import pprint
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from scipy.interpolate import CubicSpline, interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, create_power_flux_fun

LCFS_WEIGHT = 100.0
N_PSI = 1000
_NBI_W_TO_MA = 1/16e6
mu_0 = 4.0 * np.pi * 1e-7
# EQDSK sampling for save_eqdsk when validating with TORAX in _run_tm:
# start at nr=nz=100, increase by 50 up to 500 until TORAX accepts eqdsk.
EQDSK_SAVE_NR_NZ_SEQUENCE = (350, 500, 650, 800)

# Default TORAX radial face count for loop 0 coarse runs (evenly spaced normalized rho).
DEFAULT_LOOP0_TX_FACE_POINTS = 51

BASE_CONFIG = {
    'plasma_composition': {},
    'profile_conditions': {
        'normalize_n_e_to_nbar': False,
        'n_e_nbar_is_fGW': False,
        'initial_psi_from_j': False,
    },
    'numerics': {
        'dt_reduction_factor': 3,
    },
    'geometry': {},
    'sources': {
        'ohmic': {},
        'impurity_radiation': {
            'model_name': 'mavrin_fit',
            'radiation_multiplier': 1.0,
        }
    },
    'mhd': {'sawtooth': {'redistribution_model': {'model_name': 'simple'},
                        'trigger_model': {'minimum_radius': 0.1,
                                            'model_name': 'simple',
                                            's_critical': 0.4}}},
    'neoclassical': {
        'bootstrap_current': {'model_name': 'redl'},
    },
    'pedestal': {
        # 'model_name': 'set_T_ped_n_ped',
    },
    # 'transport': { # recommended in TORAX documentation, TORAX hangs at t~0 when using this config.
    #     'model_name': 'combined',
    #     'transport_models': [
    #         # Base model: QLKNN applied everywhere (default ADD)
    #         {
    #             'model_name': 'qlknn',
    #             'rho_max': 1.0,
    #         },
    #         # Edge overwrite: Sets D_e and V_e in the edge, ignoring QLKNN there.
    #         # Keeps chi_i/chi_e from QLKNN (because they are disabled here).
    #         {
    #             'model_name': 'constant',
    #             'rho_min': 0.9,
    #             'D_e': 0.5,
    #             'V_e': -1.0,
    #             'merge_mode': 'overwrite',
    #             'disable_chi_i': True,
    #             'disable_chi_e': True,
    #         },
    #     ],
    # },
    'transport': { # old base config
        'model_name': 'qlknn',
        'apply_inner_patch': True,
        'rho_inner': 0.1,
        'apply_outer_patch': False,
        'D_e_outer': 0.1,
        'V_e_outer': 0.0,
        'chi_i_outer': 2.0,
        'chi_e_outer': 2.0,
        'rho_outer': 0.95,
        'chi_min': 0.05,
        'chi_max': 100,
        'D_e_min': 0.05,
        'D_e_max': 50,
        'V_e_min': -10,
        'V_e_max': 10,
        # 'smoothing_width': 0.3,
        'DV_effective': True,
        'include_ITG': True,
        'include_TEM': True,
        'include_ETG': True,
        'avoid_big_negative_s': False,
        'smooth_everywhere': True,
        'smoothing_width': 0.3,
    },
    'solver': {
        'solver_type': 'newton_raphson',
        'use_predictor_corrector': True,
        'n_corrector_steps': 10,
        'chi_pereverzev': 30,
        'D_pereverzev': 15,
        'use_pereverzev': True,
        # Loosened from TORAX defaults (n_max=30, residual_tol=1e-5,
        # residual_coarse_tol=1e-2) since the outer TX-TM coupling loop already
        # iterates to self-consistency - per-step inner tolerance does not need
        # to be machine-precision. Benchmarked on ITER 15 MA flattop + 60 s
        # rampdown: identical Te/ne/q0/li to 4 sig figs vs strict tols, ~5x
        # faster per inner solve, 19/19 raw_tx in both coupling loops. Users
        # with very stiff scenarios can re-tighten via their own solver block.
        'n_max_iterations': 10,
        'residual_tol': 1e-3,
        'residual_coarse_tol': 1e-1,
    },
    'time_step_calculator': {
        'calculator_type': 'fixed',
    },
}


# Setup output re-direct from TORAX to log file, suppressing frivolous warnings.
# Errors will still be output in terminal.
# This is the first step, needs to be given self._log_file once that is configured in self.fly().
def log_redirect_setup():
    r'''! Step 1/3 of setup to redirect noisy outputs to log file.
        Performs the initial, minimal logging setup.
        - Removes any handlers pre-configured by libraries.
        - Sets the root logger's level to capture all desired messages.
        - Adds a console handler for critical errors only.
        
    '''
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)  # Capture INFO level and above
    
    # Remove any pre-existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add a handler to show ONLY errors on the console
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.ERROR)
    formatter = logging.Formatter('CONSOLE ERROR: [%(levelname)s] %(message)s')
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

log_redirect_setup()

# Now import "noisy" packages, after running log_redirect_setup:
import torax  # noqa: E402
# Internal TORAX APIs for running the loop ourselves so we can keep the raw
# list[SimState] (run_simulation() drops the pedestal confinement_mode).
from torax._src.pedestal_model.pedestal_transition_state import ConfinementMode  # noqa: E402
from torax._src.orchestration import run_simulation as run_sim_lib  # noqa: E402
from torax._src.orchestration import run_loop as run_loop_lib  # noqa: E402
from torax._src.output_tools import output as output_lib  # noqa: E402


@contextmanager
def redirect_outputs_to_log(filename):
    r'''! Step 2/3 of setup to redirect noisy outputs to log file. 
        A context manager to temporarily redirect stdout and stderr to a file.
        @param filename Name of log file (self._log_file)
        
    '''
    if not filename:
        # If no filename is provided, do nothing.
        yield
        return

    with open(filename, 'a') as log_file:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = log_file
        sys.stderr = log_file
        try:
            # Write a separator to the log to show where stdout redirection started
            log_file.write("\n--- [Begin capturing stdout/stderr] ---\n")
            yield
        finally:
            # Write a separator to show where it ended
            log_file.write("\n--- [End capturing stdout/stderr] ---\n")
            # Crucially, restore the original streams
            sys.stdout = original_stdout
            sys.stderr = original_stderr

class MyEncoder(json.JSONEncoder):
    '''! JSON Encoder Object to store simulation results.'''
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        # Handle xarray DataArray
        if hasattr(obj, 'to_numpy'):
            return obj.to_numpy().tolist()
        return json.JSONEncoder.default(self, obj)

class TokaMaker_TORAX:
    '''! TokaMaker + TORAX Coupled Pulse Simulation Code'''


    # ─── Initialization ─────────────────────────────────────────────────────────

    def __init__(self, t_init, t_final, eqtimes, g_eqdsk_arr, tokamaker_obj, tx_dt=0.1, tm_times=None, last_surface_factor=0.99, truncate_eq=False):
        r'''! Initialize the Coupled TokaMaker + TORAX object.
                @param t_init Start time (s).
                @param t_final End time (s).
                @param eqtimes Time points of each gEQDSK file.
                @param g_eqdsk_arr Filenames of each gEQDSK file.
                @param tokamaker_obj Preconfigured TokaMaker object.
                @param tx_dt Time step (s) of TORAX simulation.
                @param tm_times Time points where TokaMaker solves equilibrium.
                @param last_surface_factor Last surface factor for Torax.
                @param truncate_eq Whether to truncate equilibrium when saving TokaMaker output to EQDSK.
                @brief Coupling fly() option loop0: if True, the first coupling pass uses a coarse
                       TORAX grid and subsampled TokaMaker times; if False, that pass is skipped (see fly()).
                
        '''
        self._tm = tokamaker_obj
        self._cocos = 2 # only cocos=2 works currently.

        self._state = {}
        self._eqtimes = eqtimes
        self._results = {}
        self._init_files = g_eqdsk_arr
        self._t_init = t_init
        self._t_final = t_final
        self._tx_dt = tx_dt # TORAX timestep
        self._last_surface_factor = last_surface_factor
        self._psi_N = np.linspace(0.0, 1.0, N_PSI) # standardized psi_N grid all values should be mapped onto
        self._truncate_eq = truncate_eq

        self._current_loop = 0

        if tm_times is None:
            self._tm_times = sorted(eqtimes)
        else:
            self._tm_times = sorted(tm_times)

        # Drop any tm_times that fall outside the [t_init, t_final] simulation window.
        _out_of_window = [t for t in self._tm_times if t < self._t_init or t > self._t_final]
        if _out_of_window:
            self._print(f"requested eq_times outside simulation time window: {_out_of_window}")
            self._tm_times = [t for t in self._tm_times if self._t_init <= t <= self._t_final]

        # Always solve an equilibrium at the start and end of the simulation window.
        # Use a small tolerance so float round-off doesn't cause duplicate near-endpoint solves.
        _endpoint_atol = 1e-9
        if not any(abs(t - self._t_init) <= _endpoint_atol for t in self._tm_times):
            self._tm_times.append(self._t_init)
        if not any(abs(t - self._t_final) <= _endpoint_atol for t in self._tm_times):
            self._tm_times.append(self._t_final)
        self._tm_times = sorted(self._tm_times)

        N = len(self._tm_times)

        # ── Geometry scalars (seed values from EQDSK, updated by TM each loop) ──
        self._state['R0_mag']   = np.zeros(N)
        self._state['Z']        = np.zeros(N)
        self._state['a']        = np.zeros(N)
        self._state['kappa']    = np.zeros(N)
        self._state['delta']    = np.zeros(N)
        self._state['B0']       = np.zeros(N)

        # ── Global plasma scalars ────────────────────────────────────────────────
        self._state['Ip']       = np.zeros(N)   # EQDSK seed; updated by TX/TM each loop
        self._state['Ip_tm']    = np.zeros(N)   # TM-solved Ip
        self._state['Ip_tx']    = np.zeros(N)   # TX-solved Ip
        self._state['Ip_ni_tx'] = np.zeros(N)   # TX non-inductive Ip
        self._state['pax']      = np.zeros(N)   # EQDSK seed axis pressure
        self._state['pax_tm']   = np.zeros(N)   # TM-solved axis pressure
        self._state['beta_N_tm']= np.zeros(N)
        self._state['beta_N_tx']= np.zeros(N)
        self._state['beta_pol'] = np.zeros(N)
        self._state['l_i_tm']   = np.zeros(N)
        self._state['vloop_tm'] = np.zeros(N)
        self._state['vloop_tx'] = np.zeros(N)
        self._state['f_GW']     = np.zeros(N)
        self._state['f_GW_vol'] = np.zeros(N)

        # H-mode flag at each tm_time (True=H-mode, False=L-mode); filled from the
        # TORAX pedestal confinement_mode when the L/H state machine is active.
        self._state['confinement_mode'] = np.zeros(N, dtype=bool)

        # ── Safety factor scalars ────────────────────────────────────────────────
        self._state['q0']    = np.zeros(N)
        self._state['q95']   = np.zeros(N)
        self._state['q0_tm'] = np.zeros(N)
        self._state['q95_tm']= np.zeros(N)

        # ── Poloidal flux scalars ────────────────────────────────────────────────
        self._state['psi_axis_tm']     = np.zeros(N)
        self._state['psi_lcfs_tm']     = np.zeros(N)
        self._state['psi_axis_tx']     = np.zeros(N)
        self._state['psi_lcfs_tx']     = np.zeros(N)
        self._state['psi_grid_prev_tm']= np.zeros(N)  # psi on nodes from previous TM solve, for warm-starting

        # ── Safety factor profiles {i: {'x':..., 'y':..., 'type':...}} ─────────
        self._state['q_prof_eqdsk'] = {}
        self._state['q_prof_tm']    = {}
        self._state['q_prof_tx']    = {}

        # ── Poloidal flux profiles {i: {'x':..., 'y':..., 'type':...}} ─────────
        self._state['psi_tm'] = {}
        self._state['psi_tx'] = {}

        # ── LCFS geometry and GS source profiles {i: ...} ───────────────────────
        self._state['lcfs_geo']       = {}  # (N,2) LCFS contour at each timepoint
        self._state['ffp_prof']       = {}  # FF' profile sent to GS solve (updated each loop)
        self._state['pp_prof']        = {}  # p' profile sent to GS solve (updated each loop)
        self._state['ffp_ni_prof']    = {}  # FF' from non-inductive current only
        self._state['eta_prof']       = {}  # resistivity profile
        self._state['ffp_prof_eqdsk'] = {}  # seed FF' from EQDSK (normalized)
        self._state['pp_prof_eqdsk']  = {}  # seed p' from EQDSK (normalized)
        self._state['p_prof_eqdsk']   = {}  # seed pressure from EQDSK (normalized)
        self._state['ffp_prof_tx']    = {}  # FF' from TORAX
        self._state['pp_prof_tx']     = {}  # p' from TORAX
        self._state['ffp_prof_tm']    = {}  # FF' from TM solve (FF, not FF')
        self._state['pp_prof_tm']     = {}  # p' from TM solve
        self._state['p_prof_tm']      = {}  # pressure from TM solve
        self._state['p_prof_tx']      = {}  # pressure from TORAX
        self._state['f_prof_tm']      = {}  # F=RBphi from TM solve

        # ── Kinetic profiles from TORAX {i: {'x':..., 'y':..., 'type':...}} ────
        self._state['T_e']  = {}
        self._state['T_i']  = {}
        self._state['n_e']  = {}
        self._state['n_i']  = {}
        self._state['ptot'] = {}

        # ── Current density profiles from TORAX {i: ...} ────────────────────────
        self._state['j_tot']             = {}
        self._state['j_ohmic']           = {}
        self._state['j_ni']              = {}
        self._state['j_bootstrap']       = {}
        self._state['j_ecrh']            = {}
        self._state['j_external']        = {}
        self._state['j_generic_current'] = {}

        # ── Thermal conductivity profiles from TORAX {i: ...} ───────────────────
        self._state['chi_neo_e']  = {}
        self._state['chi_neo_i']  = {}
        self._state['chi_etg_e']  = {}
        self._state['chi_itg_e']  = {}
        self._state['chi_itg_i']  = {}
        self._state['chi_tem_e']  = {}
        self._state['chi_tem_i']  = {}
        self._state['chi_turb_e'] = {}
        self._state['chi_turb_i'] = {}

        # ── Particle diffusivity profiles from TORAX {i: ...} ───────────────────
        self._state['D_itg_e']  = {}
        self._state['D_neo_e']  = {}
        self._state['D_tem_e']  = {}
        self._state['D_turb_e'] = {}

        # ── Flux-surface geometry profiles {i: ...} ──────────────────────────────
        self._state['R_avg_tm']     = {}
        self._state['R_inv_avg_tm'] = {}
        self._state['R_inv_avg_tx'] = {}

        # ── TokaMaker equilibrium objects and strike points {i: ...} ────────────
        self._state['equil']      = {}  # {i: TokaMaker equilibrium object}
        self._state['strike_pts'] = {}  # {i: (N,2) [R,Z] strike points, or empty}
        self._state['lcfs_geo_tm']= {}  # {i: (N,2) LCFS contour traced from TM solve}
        self._state['x_pts_tm']   = {}  # {i: (n_xpts,2) X-point locations from TM solve, or None}

        # ── Warm-start psi: persists across loops ────────────────────────────────
        self._psi_warm_start = {}  # {i: psi_array}

        # ── Steady-state mode ────────────────────────────────────────────
        self._steady_state_mode = False
        self._steady_state_tx_seed = None
        self._steady_state_tm_psi_seed = None

        self._results['q'] = {}
        self._results['n_e'] = {}
        self._results['T_e'] = {}
        self._results['T_i'] = {}


        R = []
        Z = []  
        a = []
        kappa = []
        delta = []
        B0 = []
        pax = []
        Ip = []
        lcfs = []
        ffp_prof = []
        pp_prof = []
        q_prof = []
        psi_axis = []
        psi_lcfs = []
        pres_prof = []

        def interp_tm_lcfs(lcfs, n=N_PSI):
            r'''! Interpolates LCFS geometry points between input eqtimes to fill in all tm_times.'''
            x = lcfs[:, 0]
            y = lcfs[:, 1]
            u = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
            x[-1] = x[0]
            y[-1] = y[0]
            u = np.r_[0, u]
            cs_x = CubicSpline(u, x, bc_type='periodic')
            cs_y = CubicSpline(u, y, bc_type='periodic')
            u_new = np.linspace(u.min(), u.max(), n)
            x_new_cs = cs_x(u_new)
            y_new_cs = cs_y(u_new)
            return np.column_stack([x_new_cs, y_new_cs])

        for i, t in enumerate(self._eqtimes):
            g = read_eqdsk(g_eqdsk_arr[i])
            zmax = np.max(g['rzout'][:,1])
            zmin = np.min(g['rzout'][:,1])
            rmax = np.max(g['rzout'][:,0])
            rmin = np.min(g['rzout'][:,0])
            minor_radius = (rmax - rmin) / 2.0
            rgeo = (rmax + rmin) / 2.0
            highest_pt_idx = np.argmax(g['rzout'][:,1])
            lowest_pt_idx = np.argmin(g['rzout'][:,1])
            rupper = g['rzout'][highest_pt_idx][0]
            rlower = g['rzout'][lowest_pt_idx][0]
            delta_upper = (rgeo - rupper) / minor_radius
            delta_lower = (rgeo - rlower) / minor_radius

            R.append(g['rcentr'])
            Z.append(g['zaxis'])  # magnetic axis Z, not grid midpoint (zmid)
            a.append(minor_radius)
            kappa.append((zmax - zmin) / (2.0 * minor_radius))
            delta.append((delta_upper + delta_lower) / 2.0)
            
            B0.append(g['bcentr'])
            pax.append(g['pres'][0])
            Ip.append(abs(g['ip']))

            # Eqdsks are written via TokaMaker save_eqdsk(cocos=2), which negates
            # TM's native psi on write. Flip sign back to store TM-native in _state
            # (psi_axis > psi_lcfs, in Wb/rad). Sign is preserved (can be ±).
            psi_axis.append(-g['psimag'])
            psi_lcfs.append(-g['psibry'])

            lcfs.append(interp_tm_lcfs(g['rzout']))

            psi_eqdsk = np.linspace(0.0, 1.0, g['nr'])
            # cocos=2 file psi = -psi_TM, so d/dpsi_file = -d/dpsi_TM.
            # Flip pprime and ffprime to get TM-native derivatives.
            ffp = -np.interp(self._psi_N, psi_eqdsk, g['ffprim'])
            pp  = -np.interp(self._psi_N, psi_eqdsk, g['pprime'])
            q = np.interp(self._psi_N, psi_eqdsk, g['qpsi'])

            ffp_prof.append(ffp)
            pp_prof.append(pp)
            q_prof.append(q)

            pres = np.interp(self._psi_N, psi_eqdsk, g['pres'])
            pres_prof.append(pres)

        self.lcfs = lcfs

        def interp_prof(profs, time):
            r'''! Linearly interpolates profile between input eqtimes.'''
            if time <= self._eqtimes[0]:
                return profs[0]
            for i in range(1, len(self._eqtimes)):
                if time == self._eqtimes[i]:
                    return profs[i]
                elif time > self._eqtimes[i-1] and time <= self._eqtimes[i]:
                    dt = self._eqtimes[i] - self._eqtimes[i-1]
                    alpha = (time - self._eqtimes[i-1]) / dt
                    return (1.0 - alpha) * np.array(profs[i-1]) + alpha * np.array(profs[i])
            return profs[-1]

        for i, t in enumerate(self._tm_times):
            # Default Scalars
            self._state['R0_mag'][i] = np.interp(t, self._eqtimes, R)
            self._state['Z'][i] = np.interp(t, self._eqtimes, Z)
            self._state['a'][i] = np.interp(t, self._eqtimes, a)
            self._state['kappa'][i] = np.interp(t, self._eqtimes, kappa)
            self._state['delta'][i] = np.interp(t, self._eqtimes, delta)
            self._state['B0'][i] = np.interp(t, self._eqtimes, B0)
            self._state['pax'][i] = np.interp(t, self._eqtimes, pax)
            self._state['Ip'][i] = np.interp(t, self._eqtimes, Ip)
            self._state['psi_axis_tm'][i] = np.interp(t, self._eqtimes, psi_axis)
            self._state['psi_lcfs_tm'][i] = np.interp(t, self._eqtimes, psi_lcfs)

            # Default Profiles
            self._state['lcfs_geo'][i] = interp_prof(lcfs, t)
            self._state['ffp_prof'][i] = {'x': self._psi_N.copy(), 'y': interp_prof(ffp_prof, t), 'type': 'linterp'}
            self._state['pp_prof'][i] = {'x': self._psi_N.copy(), 'y': interp_prof(pp_prof, t), 'type': 'linterp'}
            self._state['ffp_ni_prof'][i] = {'x': [], 'y': [], 'type': 'linterp'}

            self._state['ffp_prof_eqdsk'][i] = self._state['ffp_prof'][i].copy()
            self._state['pp_prof_eqdsk'][i] = self._state['pp_prof'][i].copy()
            self._state['p_prof_eqdsk'][i] = {'x': self._psi_N.copy(), 'y': interp_prof(pres_prof, t), 'type': 'linterp'}
            self._state['pp_prof_eqdsk'][i]['y'] /= self._state['pp_prof_eqdsk'][i]['y'][0]
            self._state['p_prof_eqdsk'][i]['y'] /= self._state['p_prof_eqdsk'][i]['y'][0]
            self._state['ffp_prof_eqdsk'][i]['y'] /= self._state['ffp_prof_eqdsk'][i]['y'][0]
            self._state['q_prof_eqdsk'][i] = {'x': self._psi_N.copy(), 'y': interp_prof(q_prof, t), 'type': 'linterp'}

            self._state['eta_prof'][i]= {
                'x': self._psi_N.copy(),
                'y': np.zeros(N_PSI),
                'type': 'linterp',
            }
            
        # Save seed values from initial equilibria
        self._psi_axis_seed = self._state['psi_axis_tm'].copy()
        self._psi_lcfs_seed = self._state['psi_lcfs_tm'].copy()
        self._Ip_seed       = self._state['Ip'].copy()
        self._pax_seed      = self._state['pax'].copy()
        self._state['pax_tm'] = self._state['pax'].copy()

        self._psi_init = None
        self._n_e_init = None
        self._T_e_init = None
        self._T_i_init = None
        self._relax_profiles_snapshot = None  # buffer while appending to relax debug history
        # Main TORAX snapshots at t_init after each completed coupling loop (for relax diagnostics).
        self._relax_mainrun_profile_history = []

        self._Ip = None
        self._Zeff = None

        self._nbi_heating = None
        self._generic_heat = None
        self._generic_heat_loc = None
        self._generic_heat_width = None
        self._use_nbi_current = False
        self._ecrh_heating = None
        self._ecrh_loc = None
        self._ecrh_width = None
        self._nbi_loc = None
        self._ohmic_power = None
        self._nbar = None
        self._n_e = None
        self._T_i = None
        self._T_e = None

        self._set_pedestal = None
        self._T_i_ped = None
        self._T_e_ped = None
        self._n_e_ped = None
        self._ped_top = None

        # Optional full replacement for the TORAX pedestal section.
        self._pedestal_config = None

        # ── Pedestal config mode (set_pedestal) ──────────────────────────────
        # 'ADAPTIVE_SOURCE' | 'ADAPTIVE_TRANSPORT' | 'internal_manual'.
        self._ped_mode = 'off'
        # Enable the L/H formation model under ADAPTIVE_SOURCE (new-API pedestals only).
        self._ped_formation_model = False
        # internal_manual pedestal knobs (heights in keV / m^-3; widths/times in rho_norm / s).
        self._ped_height_Te = None
        self._ped_height_Ti = None
        self._ped_height_ne = None
        self._ped_width = 0.15
        self._ped_transition_time = 0.5
        self._ped_timing = 'detect'      # 'detect' (loop-1 ADAPTIVE_SOURCE) or 'input'
        self._lh_time = None             # L->H transition time (s); set by user or detection
        self._hl_time = None             # H->L back-transition time (s); None if none
        self._ped_formation_model_name = 'delabie_scaling'  # L/H formation model for detect loop
        # IBC penalty stiffness (TORAX numerics.adaptive_{T,n}_source_prefactor); None = TORAX default.
        # IBCs are soft penalties: higher = evolved tracks the imposed pedestal more exactly,
        # lower = TORAX transport smooths the imposed shape more.
        self._ped_T_source_prefactor = None
        self._ped_n_source_prefactor = None
        # Per-tm-time inner-edge (value, slope) matched to the evolved core after the first IBC
        # loop, {field: {tm_time: (value, slope)}}; empty until _smooth_ibc_inner_edge runs.
        self._ped_smoothed = {}
        # Shape params used to build the current loop's IBC (for the debug figure); None until built.
        self._ped_ibc_snapshot = None
        # Core-band L->H transition IBC (internal_manual): when a core_height is set, a band over
        # the core (rho 0->ped_rho) ramps the evolved L-mode core into a smooth H-mode-like core
        # shape during [lh, lh+transition], then retreats to the pedestal-only band. Off if None.
        self._core_height_Te = None
        self._core_height_Ti = None
        self._core_height_ne = None
        self._core_exp_a = 2.0    # core-shape exponents: f = ped + (core-ped)*(1-(rho/ped_rho)^a)^b
        self._core_exp_b = 2.0    # a>=2,b>=2 -> flat at core AND flat into the pedestal top
        # Evolved L-mode core near t=lh, {field: {rho: value}} on the core grid; captured per loop.
        self._ped_lmode_profile = {}
        # Relaxed on-axis value read back after the transition retreats, {field: value}; used as the
        # core target on the next loop (loop 1 uses the user core_height_*).
        self._ped_core_relaxed = {}

        self._Te_right_bc = None
        self._Ti_right_bc = None
        self._ne_right_bc = None

        self._gp_s = None
        self._gp_dl = None

        self._pellet_deposition_location = None
        self._pellet_width = None
        self._pellet_s_total = None

        # Transport / numerics overrides — None means "use value from
        # loaded config (or base config if no loaded config)".  Only set
        # to a real value when an explicit set_*() call is made AFTER
        # load_TORAX_config().
        self._normalize_to_nbar = None
        self._evolve_density = None
        self._evolve_current = None
        self._evolve_Ti = None
        self._evolve_Te = None
        self._chi_min = None
        self._chi_max = None
        self._De_min = None
        self._De_max = None
        self._Ve_min = None
        self._Ve_max = None

        self._main_ion = None
        self._impurity = None
        self._enable_fusion = True
        self._enable_ei_exchange = True

        self._loaded_config = None   # set by load_TORAX_config()
        self._tx_grid_type = None
        self._tx_grid = None
        self._tx_grid_stash = []  # stack of (type, grid) for loop-1 coarse TORAX grid save/restore

        self._diverted_times = None
        self._x_point_targets = None
        self._x_point_weight = 100.0
        self._strike_point_targets = None
        self._trim_lcfs = True

        self._eqdsk_skip = []

        # Psi snapshots for movie generation (populated during _run_tm)
        self._tm_psi_on_nodes = {}  # {loop: {i: psi_array}}

        # Temp/output directory state (set in fly())
        self._out_dir = None
        self._eqdsk_dir = None
        self._eqdsk_dir_is_temp = False
        self._save_outputs = False
        self._debug_mode = False
        self._output_mode = False
        self._output_file_tag = None
        self._run_timestamp = None
        self._diagnostics = False
        self._logging_configured = False
        self._log_file = None

    # ─── Static Utilities ───────────────────────────────────────────────────────

    @staticmethod
    def _numpy_to_plain_python(obj):
        r'''! Recursively convert numpy scalars/arrays to plain Python types.
        
                Used before pformat-saving config dicts so the saved .py files are
                loadable without numpy (no array([...]) references).
                
        '''
        if isinstance(obj, dict):
            return {TokaMaker_TORAX._numpy_to_plain_python(k): TokaMaker_TORAX._numpy_to_plain_python(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            converted = [TokaMaker_TORAX._numpy_to_plain_python(v) for v in obj]
            return type(obj)(converted)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    @staticmethod
    def _tx_config_merge(base, override):
        r'''! Recursively merge override into base TORAX config (in-place).
        
                For every key in override:
                  - If both values are dicts, recurse.
                  - Otherwise the override value wins.
                Keys in base that are absent from override are kept as-is.
        
                @param base     Dict to merge into (modified in-place).
                @param override Dict whose keys take precedence.
                @return base (for convenience).
                
        '''
        for key, val in override.items():
            if key in base and isinstance(base[key], dict) and isinstance(val, dict):
                TokaMaker_TORAX._tx_config_merge(base[key], val)
            else:
                base[key] = copy.deepcopy(val)
        return base

    @staticmethod
    def _flatten_time_dependent(config):
        r'''! Recursively flatten time-dependent config values to their initial value only.
        
                A dict whose keys are ALL numeric (int/float) is treated as time-dependent:
                only the entry with the smallest key is retained.
                Structural dicts (with any string keys) are recursed into.
                Tuples of form (times_array, values_array) are flattened to the first entry.
                
        '''
        for key in list(config.keys()):
            val = config[key]
            if isinstance(val, dict):
                if val and all(isinstance(k, (int, float, np.integer, np.floating)) for k in val.keys()):
                    # Time-dependent dict: keep only the first (smallest t) entry.
                    # Do NOT recurse — the nested value may be a rho-profile dict.
                    first_key = min(val.keys())
                    config[key] = {first_key: copy.deepcopy(val[first_key])}
                else:
                    # Structural dict with string keys: recurse.
                    TokaMaker_TORAX._flatten_time_dependent(val)
            elif isinstance(val, (tuple, list)) and len(val) == 2:
                t_arr, v_arr = val
                try:
                    if (hasattr(t_arr, '__len__') and hasattr(v_arr, '__len__')
                            and not isinstance(t_arr, str) and not isinstance(v_arr, str)
                            and len(t_arr) > 1
                            and max(t_arr) > 1.0):  # rho grids have max==1.0; skip them
                        first_t = t_arr[0]
                        first_v = v_arr[0]
                        if isinstance(t_arr, np.ndarray):
                            config[key] = (np.array([first_t]), np.array([first_v]))
                        elif isinstance(t_arr, tuple):
                            config[key] = ((first_t,), (first_v,))
                        else:
                            config[key] = ([first_t], [first_v])
                except TypeError:
                    pass
            elif isinstance(val, (tuple, list)) and len(val) == 3:
                t_arr, rho_arr, v_arr = val
                try:
                    if (hasattr(t_arr, '__len__') and hasattr(rho_arr, '__len__') and hasattr(v_arr, '__len__')
                            and not isinstance(t_arr, str) and not isinstance(v_arr, str)
                            and len(t_arr) > 1
                            and max(t_arr) > 1.0):
                        first_t = t_arr[0]
                        first_v = v_arr[0]
                        if isinstance(t_arr, np.ndarray):
                            config[key] = (np.array([first_t]), rho_arr, np.array([first_v]))
                        elif isinstance(t_arr, tuple):
                            config[key] = ((first_t,), rho_arr, (first_v,))
                        else:
                            config[key] = ([first_t], rho_arr, [first_v])
                except TypeError:
                    pass

    @staticmethod
    def _relax_flat_profile_to_rho_y(profile_val):
        r'''! (rho_norm, y) from a profile_conditions value after merge + flatten.
        
                Handles:
                  - nested {time: inner} (TORAX time-sliced profiles; flatten keeps one key);
                  - (rho, y) static radial tuples;
                  - (times, rho, [profiles]) including a single time index.
                
        '''
        if profile_val is None:
            return None, None

        if isinstance(profile_val, dict) and profile_val:
            keys = list(profile_val.keys())
            if keys and all(isinstance(k, (int, float, np.integer, np.floating)) for k in keys):
                sample_vals = list(profile_val.values())
                all_scalar = all(
                    isinstance(v, (int, float, np.integer, np.floating))
                    for v in sample_vals
                )
                # TORAX set_ne / set_Te format: {t: {psi_N: value, ...}} after
                # flatten has one time key; inner dict is *only* numeric keys in [0, 1].
                if all_scalar:
                    rho_try = np.array(sorted(keys), dtype=float)
                    kmin, kmax = float(np.min(rho_try)), float(np.max(rho_try))
                    if kmax <= 1.0 + 1e-5 and kmin >= -1e-5:
                        y_try = np.array([float(profile_val[k]) for k in sorted(keys)], dtype=float)
                        if rho_try.size >= 1:
                            return rho_try, y_try
                    # e.g. Ip: {0.: ..., 5.: ...}: not a radial table — peel one time
                inner = profile_val[min(keys)]
                return TokaMaker_TORAX._relax_flat_profile_to_rho_y(inner)

        if isinstance(profile_val, (tuple, list)) and len(profile_val) == 2:
            rho, y = profile_val
            try:
                rho = np.asarray(rho, dtype=float).ravel()
                y = np.asarray(y, dtype=float).ravel()
                if rho.size == 0 or y.size == 0:
                    return None, None
                if rho.shape == y.shape:
                    return rho, y
            except (TypeError, ValueError):
                pass
            return None, None

        if isinstance(profile_val, (tuple, list)) and len(profile_val) == 3:
            _t, rho, vs = profile_val
            try:
                rho = np.asarray(rho, dtype=float).ravel()
                if isinstance(vs, np.ndarray):
                    if vs.ndim == 2 and vs.shape[0] >= 1:
                        y = np.asarray(vs[0], dtype=float).ravel()
                    else:
                        y = np.asarray(vs, dtype=float).ravel()
                elif isinstance(vs, (list, tuple)) and len(vs) > 0:
                    y = np.asarray(vs[0], dtype=float).ravel()
                else:
                    y = np.asarray(vs, dtype=float).ravel()
                if rho.size == 0 or y.size == 0:
                    return None, None
                if rho.shape == y.shape:
                    return rho, y
            except (TypeError, ValueError):
                pass

        return None, None


    # ─── Setup & Configuration ──────────────────────────────────────────────────

    def load_TORAX_config(self, config):
        r'''! Load a TORAX config dict.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html
                The loaded config is deep-merged on top of BASE_CONFIG when the
                simulation config is built.  Any key present in the loaded config
                will override the corresponding BASE_CONFIG key; keys only in
                BASE_CONFIG are kept as-is.  Geometry is always overwritten (eqdsk-based).
        
                Explicit set_*() calls made AFTER load_TORAX_config() will
                override both the base and the loaded config.
        
                @param config Dictionary (TORAX config format).
                
        '''
        self._loaded_config = copy.deepcopy(config)

    def set_TORAX_grid(self, grid_type, grid):
        r'''! Set TORAX grid type and grid points.
                TORAX grid documentation: https://torax.readthedocs.io/en/latest/configuration.html#geometry
                @param grid_type Grid type ('n_rho' or 'face_centers').
                @param grid Grid points (integer or np.array).
                
        '''
        self._tx_grid_type = grid_type
        self._tx_grid = grid
        if grid_type not in ['n_rho', 'face_centers']:
            raise ValueError(f'Invalid grid type: {type}. Must be "n_rho" or "face_centers".')

    def _push_tx_grid(self): 
        r'''! Save current TORAX grid onto stack for later _pop_tx_grid. Used for reducing grid in loop0.'''
        g = self._tx_grid
        if isinstance(g, np.ndarray):
            g = g.copy()
        elif isinstance(g, list):
            g = copy.deepcopy(g)
        self._tx_grid_stash.append((self._tx_grid_type, g))

    def _pop_tx_grid(self):
        r'''! Restore TORAX grid from _push_tx_grid stack. Used for reducing grid in loop0.'''
        self._tx_grid_type, self._tx_grid = self._tx_grid_stash.pop()

    def _coupling_iteration_is_first(self):
        r'''! True on the first TokaMaker_TORAX loop (TORAX uses the seed gEQDSK map).'''
        if self._current_loop == 0:
            return True
        if not getattr(self, '_fly_loop0', True) and self._current_loop == 1:
            return True
        return False

    @staticmethod
    def _loop0_tm_solve_indices(n_tm, stride):
        r'''! TokaMaker timestep indices for coupling loop 0 subsampling: endpoints plus every stride-th.'''
        if n_tm <= 0:
            return []
        if stride <= 1:
            return list(range(n_tm))
        return sorted({0, n_tm - 1}.union(range(0, n_tm, stride)))

    @contextmanager
    def _loop0_coarse_tx_relax_scope(self, stage):
        r'''! Apply loop-0 coarse TORAX grid during initial relax only; restore after.'''
        active = (
            stage == 'initial'
            and getattr(self, '_loop0_relax_coarse_tx', True)
            and getattr(self, '_loop0_coarse_tx', True)
        )
        if active:
            self._push_tx_grid()
            npt = int(getattr(self, '_loop0_tx_face_points', DEFAULT_LOOP0_TX_FACE_POINTS))
            self.set_TORAX_grid('face_centers', np.linspace(0.0, 1.0, npt))
        try:
            yield
        finally:
            if active:
                self._pop_tx_grid()

    @contextmanager
    def _loop0_coarse_tx_main_scope(self):
        r'''! Apply loop-0 coarse TORAX grid for main _run_tx; restore after.'''
        active = self._current_loop == 0 and getattr(self, '_loop0_coarse_tx', True)
        if active:
            self._push_tx_grid()
            npt = int(getattr(self, '_loop0_tx_face_points', DEFAULT_LOOP0_TX_FACE_POINTS))
            self.set_TORAX_grid('face_centers', np.linspace(0.0, 1.0, npt))
        try:
            yield
        finally:
            if active:
                self._pop_tx_grid()

    def set_TokaMaker_coil_reg(self, coil_bounds=None, updownsym=False,
                     default_weight=1.0E-1, disable_coils=None,
                     disable_weight=1.0E4, symmetry_weight=1.0E3,
                     disable_virtual_vsc=True, vsc_weight=1.0E4):
        r'''! Set coil regularization using the dict-based TokaMaker reg_terms API.
        
                Coil bounds are hard current limits in Amperes per turn (A/turn). The total
                current in a coil region is I_coil [A/turn] * n_turns, where n_turns comes
                from the mesh file coil definition.
        
                During the pulse, coil current targets are set automatically: the initial
                equilibrium currents seed i=0, and each subsequent timestep uses the previous
                timestep's solved currents as loose targets.
        
                @param coil_bounds Dict of {coil_name: [min, max]} hard current bounds [A/turn]. Default ±5 MA/turn.
                @param updownsym Enforce up-down symmetry for coil pairs (U/L naming convention).
                @param default_weight Regularization weight for normal coils (default 0.1).
                @param disable_coils List of coil name prefixes to disable (e.g. ['DV1', 'DV2']).
                @param disable_weight Regularization weight for disabled coils (default 1e4).
                @param symmetry_weight Regularization weight for symmetry constraints (default 1e3).
                @param disable_virtual_vsc Disable the virtual VSC coil (default True).
                @param vsc_weight Regularization weight for disabled VSC (default 1e4).
                
        '''
        if coil_bounds is None:
            coil_bounds = {key: [-5.0E6, 5.0E6] for key in self._tm.coil_sets}
        self._tm.set_coil_bounds(coil_bounds)
        self._coil_bounds = coil_bounds  # store for re-application after solve

        if disable_coils is None:
            disable_coils = []

        self._coil_reg_config = {
            'coil_bounds': coil_bounds,
            'updownsym': updownsym, 'default_weight': default_weight,
            'disable_coils': disable_coils, 'disable_weight': disable_weight,
            'symmetry_weight': symmetry_weight,
            'disable_virtual_vsc': disable_virtual_vsc, 'vsc_weight': vsc_weight,
        }
        self._apply_tm_coil_reg(targets=None)

    def _apply_tm_coil_reg(self, targets=None):
        r'''! Internal: build and apply reg_terms from stored config plus per-timestep targets.
                @param targets Dict of {coil_name: current [A/turn]} to use as soft targets, or None for zeros.
                
        '''
        cfg = self._coil_reg_config
        updownsym    = cfg['updownsym']
        default_weight = cfg['default_weight']
        disable_coils  = cfg['disable_coils']
        disable_weight = cfg['disable_weight']
        symmetry_weight = cfg['symmetry_weight']
        disable_virtual_vsc = cfg['disable_virtual_vsc']
        vsc_weight   = cfg['vsc_weight']

        reg_terms = []
        processed_for_symmetry = set()

        for name in self._tm.coil_sets:
            if name in processed_for_symmetry:
                continue

            t_target = targets[name] if (targets and name in targets) else 0.0

            if updownsym and 'U' in name:
                # Enforce up-down symmetry: I_upper - I_lower = 0
                lower_name = name.replace('U', 'L')
                if lower_name in self._tm.coil_sets:
                    reg_terms.append(self._tm.coil_reg_term(
                        {name: 1.0, lower_name: -1.0}, target=0.0, weight=symmetry_weight))
                    processed_for_symmetry.add(name)
                    processed_for_symmetry.add(lower_name)
                    continue

            weight = disable_weight if any(name.startswith(p) for p in disable_coils) else default_weight
            reg_terms.append(self._tm.coil_reg_term({name: 1.0}, target=t_target, weight=weight))

        # Virtual VSC coil
        vsc_w = vsc_weight if disable_virtual_vsc else default_weight
        reg_terms.append(self._tm.coil_reg_term({'#VSC': 1.0}, target=0.0, weight=vsc_w))

        self._tm.set_coil_reg(reg_terms=reg_terms)


    # ─── Property Setters ───────────────────────────────────────────────────────

    def set_Ip(self, Ip):
        r'''! Set plasma current (Amps), used for both codes and not evolved by either.
                @param ip Plasma current.
                
        '''
        self._Ip = Ip

    def set_ne(self, n_e, right_bc=None):
        r'''! Set electron density profiles and optional right boundary condition.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#profile-conditions
                @param n_e Electron density (m^-3).
                @param right_bc Right boundary value at rho=1 (m^-3), scalar or time-varying map.
        '''
        self._n_e = n_e
        if right_bc is not None:
            self._ne_right_bc = right_bc

    def set_Te(self, T_e, right_bc=None):
        r'''! Set electron temperature profiles and optional right boundary condition.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#profile-conditions
                @param T_e Electron temperature (keV).
                @param right_bc Right boundary value at rho=1 (keV), scalar or time-varying map.
        '''
        self._T_e = T_e
        if right_bc is not None:
            self._Te_right_bc = right_bc

    def set_Ti(self, T_i, right_bc=None):
        r'''! Set ion temperature profiles and optional right boundary condition.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#profile-conditions
                @param T_i Ion temperature (keV).
                @param right_bc Right boundary value at rho=1 (keV), scalar or time-varying map.
        '''
        self._T_i = T_i
        if right_bc is not None:
            self._Ti_right_bc = right_bc

    def set_plasma_composition(self, Zeff=None, main_ion=None, impurity=None):
        r'''! Set plasma effective charge and composition.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#profile-conditions
                @param Zeff Effective charge target.
                @param main_ion Main ion species dict, e.g. {'D': 0.5, 'T': 0.5}.
                @param impurity Impurity species string, e.g. 'Ne', 'Ar', 'W'.
        '''
        if Zeff is not None:
            self._Zeff = Zeff
        if main_ion is not None:
            self._main_ion = main_ion
        if impurity is not None:
            self._impurity = impurity

    def set_heating(self, generic_heat=None, generic_heat_loc=None, generic_heat_width=0.25, generic_heat_absorption_fraction=1.0, generic_heat_electron_heat_fraction=0.4, nbi_current=False, ecrh=None, ecrh_loc=None, ecrh_width=0.1, ohmic=None, fusion=True, ei_exchange=True):
        r'''! Set TORAX heating and source toggles.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#sources
                @param generic_heat Generic heating ({time: power_W}).
                @param generic_heat_loc Generic heating deposition location (normalized rho).
                @param generic_heat_width Generic heating deposition width (normalized rho).
                @param generic_heat_absorption_fraction Fraction of generic heating that is absorbed (default here 1.0, TORAX default is 0.0).
                @param generic_heat_electron_heat_fraction Fraction of generic heating that is converted to electron heat (default 0.4, from arc v3a config, TORAX default is 0.6).
                @param nbi_current Enable NBI current drive estimate from generic heating.
                @param ecrh ECRH heating ({time: power_W}).
                @param ecrh_loc ECRH deposition location (normalized rho).
                @param ecrh_width ECRH deposition width (normalized rho).
                @param ohmic Optional explicit ohmic power ({time: power_W}).
                @param fusion Enable fusion alpha heating source.
                @param ei_exchange Enable electron-ion energy exchange source.
        '''
        if generic_heat is not None and generic_heat_loc is not None:
            self._generic_heat = generic_heat
            self._generic_heat_loc = generic_heat_loc
            self._generic_heat_width = generic_heat_width
            self._generic_heat_absorption_fraction = generic_heat_absorption_fraction
            self._generic_heat_electron_heat_fraction = generic_heat_electron_heat_fraction
        if ecrh is not None and ecrh_loc is not None:
            self._ecrh_heating = ecrh
            self._ecrh_loc = ecrh_loc
            self._ecrh_width = ecrh_width
        if ohmic is not None:
            self._ohmic_power = ohmic
        
        self._use_nbi_current = nbi_current
        self._enable_fusion = fusion
        self._enable_ei_exchange = ei_exchange

    def set_pedestal(self, config_mode='off', config=None,
                     ped_height_Te=None, ped_height_Ti=None, ped_height_ne=None,
                     ped_width=0.15, transition_time=0.5, timing='detect',
                     lh_time=None, hl_time=None, formation_model='delabie_scaling',
                     T_source_prefactor=None, n_source_prefactor=None,
                     core_height_Te=None, core_height_Ti=None, core_height_ne=None,
                     core_exp_a=2.0, core_exp_b=2.0,
                     set_pedestal=None, T_i_ped=None, T_e_ped=None, n_e_ped=None, ped_top=0.90):
        r'''! Configure the TORAX pedestal model.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#pedestal

                Three modes via config_mode:
                - 'off'                   : no pedestal (default).
                - 'ADAPTIVE_SOURCE'    : TORAX set_T_ped_n_ped pedestal in ADAPTIVE_SOURCE mode,
                                            with use_formation_model_with_adaptive_source enabled by
                                            default (gives the L/H state machine). This is the legacy
                                            behavior; pass set_pedestal/T_*_ped/n_e_ped/ped_top as before.
                - 'ADAPTIVE_TRANSPORT' : same pedestal model in ADAPTIVE_TRANSPORT mode.
                - 'internal_manual'       : pedestal imposed by us as a TORAX internal boundary
                                            condition (IBC), not a TORAX pedestal model. The pedestal
                                            top sits at rho = 1 - ped_width and an mtanh shape is imposed
                                            out to rho=1, ramped over transition_time across each L<->H
                                            edge. Edge (rho=1) values come from set_ne/set_Te/set_Ti.
                                            If a core_height_* is given, the IBC additionally spans the
                                            CORE (rho 0->ped_rho) during the L->H transition only,
                                            ramping the TORAX-evolved L-mode core into a smooth
                                            H-mode-like core shape, then retreating to the pedestal band.

                @param config_mode 'off', 'ADAPTIVE_SOURCE', 'ADAPTIVE_TRANSPORT', or 'internal_manual'.
                @param config Optional full pedestal-section dict (tx_adaptive modes only); overrides defaults.
                @param ped_height_Te Electron-temperature pedestal-top value (keV), internal_manual mode.
                @param ped_height_Ti Ion-temperature pedestal-top value (keV), internal_manual mode (defaults to Te).
                @param ped_height_ne Electron-density pedestal-top value (m^-3), internal_manual mode.
                @param ped_width Pedestal width in rho_norm; inner edge is at 1 - ped_width.
                @param transition_time Ramp duration (s) over which the pedestal rises from L-mode.
                @param timing 'detect' (find L/H times from a loop-1 ADAPTIVE_SOURCE run) or 'input'.
                @param lh_time L->H transition time (s); required when timing='input'.
                @param hl_time H->L back-transition time (s); None means the pedestal stays on.
                @param formation_model L/H formation-model name used by the internal_manual 'detect' loop
                                       (e.g. 'delabie_scaling', 'martin_scaling').
                @param T_source_prefactor Stiffness of the internal_manual T_e/T_i IBC (TORAX
                                       numerics.adaptive_T_source_prefactor). The IBC is a soft
                                       penalty, not a hard BC: higher = the evolved T tracks the
                                       imposed pedestal more exactly; lower = TORAX transport smooths
                                       the imposed shape more. None = TORAX default (2e10).
                @param n_source_prefactor Same for the n_e IBC (numerics.adaptive_n_source_prefactor;
                                       TORAX default 2e8). None = TORAX default.
                @param core_height_Te On-axis (rho=0) H-mode T_e target (keV) for the full-domain L->H
                                       transition IBC; None disables the full-domain ramp for T_e.
                                       Loop 1 uses this value; loops 2+ use the relaxed evolved core.
                @param core_height_Ti On-axis H-mode T_i target (keV); None disables it for T_i.
                @param core_height_ne On-axis H-mode n_e target (m^-3); None disables it for n_e.
                @param core_exp_a/core_exp_b Core-shape exponents for the transition core profile
                                       f = ped + (core-ped)*(1-(rho/ped_rho)^a)^b. a>=2,b>=2 (default
                                       2,2) give a flat core and a flat hand-off into the pedestal top.
                @param set_pedestal Legacy tx_adaptive toggle (True=set_T_ped_n_ped, False=no_pedestal).
                @param T_i_ped Legacy ion-temperature pedestal (tx_adaptive).
                @param T_e_ped Legacy electron-temperature pedestal (tx_adaptive).
                @param n_e_ped Legacy electron-density pedestal (tx_adaptive).
                @param ped_top Legacy pedestal-top location rho_norm_ped_top (tx_adaptive).
        '''
        if config_mode not in ('off', 'ADAPTIVE_SOURCE', 'ADAPTIVE_TRANSPORT', 'internal_manual'):
            raise ValueError("config_mode must be 'off', 'ADAPTIVE_SOURCE', 'ADAPTIVE_TRANSPORT', or 'internal_manual'.")
        self._ped_mode = config_mode

        if config_mode == 'off':
            # No pedestal: clear both the tx_adaptive and internal_manual state so a prior
            # set_pedestal() call doesn't leak a pedestal section into the config.
            self._pedestal_config = None
            self._set_pedestal = None
            self._ped_formation_model = False
            return

        if config_mode == 'internal_manual':
            if timing not in ('detect', 'input'):
                raise ValueError("timing must be 'detect' or 'input'.")
            if timing == 'input' and lh_time is None:
                raise ValueError("timing='input' requires lh_time.")
            self._pedestal_config = None
            self._set_pedestal = None
            self._ped_height_Te = ped_height_Te
            self._ped_height_Ti = ped_height_Ti if ped_height_Ti is not None else ped_height_Te
            self._ped_height_ne = ped_height_ne
            self._ped_width = ped_width
            self._ped_transition_time = transition_time
            self._ped_timing = timing
            self._lh_time = lh_time
            self._hl_time = hl_time
            self._ped_formation_model_name = formation_model
            self._ped_T_source_prefactor = T_source_prefactor
            self._ped_n_source_prefactor = n_source_prefactor
            self._ped_smoothed = {}
            # Core-band L->H transition IBC (optional, per field via core_height_*).
            self._core_height_Te = core_height_Te
            self._core_height_Ti = core_height_Ti if core_height_Ti is not None else core_height_Te
            self._core_height_ne = core_height_ne
            self._core_exp_a = core_exp_a   # core-shape exponents (>=2 = flat core & flat into ped)
            self._core_exp_b = core_exp_b
            self._ped_lmode_profile = {}
            self._ped_core_relaxed = {}
            return

        # ── ADAPTIVE_SOURCE / ADAPTIVE_TRANSPORT ──
        # A full pedestal dict overrides everything else (former load_pedestal_config).
        self._pedestal_config = copy.deepcopy(config) if config is not None else None
        # Legacy default: pedestal off unless explicitly enabled.
        self._set_pedestal = False if set_pedestal is None else set_pedestal
        # Enable the L/H formation model for ADAPTIVE_SOURCE via the new API only; the
        # legacy set_pedestal=... toggle keeps TORAX's default (off) so old runs are unchanged.
        self._ped_formation_model = set_pedestal is None
        if T_i_ped is not None:
            self._T_i_ped = T_i_ped
        if T_e_ped is not None:
            self._T_e_ped = T_e_ped
        if n_e_ped is not None:
            self._n_e_ped = n_e_ped
        self._ped_top = ped_top

    def load_pedestal_config(self, pedestal_config):
        r'''! Deprecated: use set_pedestal(config=...). Fully replaces the TORAX pedestal section.
                @param pedestal_config Dictionary in TORAX pedestal config format, or None to clear.
        '''
        if pedestal_config is None:
            self._pedestal_config = None
            return
        if not isinstance(pedestal_config, dict):
            raise TypeError('pedestal_config must be a dictionary or None.')
        self._pedestal_config = copy.deepcopy(pedestal_config)

    def set_evolve(self, density=True, Ti=True, Te=True, current=True):
        r'''! Set variables as either prescribed (False) or evolved (True) for TORAX.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#numerics
                @param density Evolve density.
                @param Ti Evolve ion temperature.
                @param Te Evolve electron temperature.
                @param current Evolve current.
                
        '''
        self._evolve_density = density
        self._evolve_current = current
        self._evolve_Ti = Ti
        self._evolve_Te = Te

    def set_fueling(self, gas_puff_S_total=None, gas_puff_decay_length=None, pellet_deposition_location=None, pellet_width=None, pellet_S_total=None,
                    generic_particle_location=None, generic_particle_width=None, generic_particle_S_total=None):
        r'''! Set gas puff and pellet fueling particle sources for TORAX.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#sources
                @param gas_puff_S_total Gas puff particle source (particles/s).
                @param gas_puff_decay_length Gas puff decay length from edge (normalized rho).
                @param pellet_deposition_location Pellet deposition location (normalized rho).
                @param pellet_width Pellet deposition width (normalized rho).
                @param pellet_S_total Pellet particle source (particles/s).
                @param generic_particle_location Generic particle source location (normalized rho).
                @param generic_particle_width Generic particle source width (normalized rho).
                @param generic_particle_width Generic particle source amount (particles/s).
        '''
        self._gp_s = gas_puff_S_total
        self._gp_dl = gas_puff_decay_length
        self._pellet_deposition_location = pellet_deposition_location
        self._pellet_width = pellet_width
        self._pellet_s_total = pellet_S_total

        if [generic_particle_location, generic_particle_width, generic_particle_S_total].count(None) in [1,2]:
            raise ValueError("Must specify all three generic particle parameters or none of them.")
        
        self._generic_particle_location = generic_particle_location
        self._generic_particle_width = generic_particle_width
        self._generic_particle_s_total = generic_particle_S_total

    def set_transport_coefs(self, chi_min=None, chi_max=None, De_min=None, De_max=None, Ve_min=None, Ve_max=None):
        r'''! Set transport coefficient bounds for TORAX.
                TORAX input config documentation: https://torax.readthedocs.io/en/latest/configuration.html#transport
                @param chi_min Minimum ion thermal diffusivity (m^2/s).
                @param chi_max Maximum ion thermal diffusivity (m^2/s).
                @param De_min Minimum electron diffusion coefficient (m^2/s).
                @param De_max Maximum electron diffusion coefficient (m^2/s).
                @param Ve_min Minimum electron thermal velocity (m/s).
                @param Ve_max Maximum electron thermal velocity (m/s).
        '''
        if chi_min is not None:
            self._chi_min = chi_min
        if chi_max is not None:
            self._chi_max = chi_max
        if De_min is not None:
            self._De_min = De_min
        if De_max is not None:
            self._De_max = De_max
        if Ve_min is not None:
            self._Ve_min = Ve_min
        if Ve_max is not None:
            self._Ve_max = Ve_max

    def set_x_points(self, diverted_times=None, x_point_targets=None, x_point_weight=100.0, strike_point_targets=None,
                     trim_lcfs=True, trim_lcfs_perc_limit=0.80):
        r'''! Configure diverted window, X-point targets, and optional strike points.

                @param diverted_times Tuple (t_start, t_end) defining the diverted plasma window.
                @param x_point_targets X-point target locations, shape (n_xpoints, 2) with [R, Z] pairs.
                @param x_point_weight Weight for saddle-point constraints.
                @param strike_point_targets Strike point locations, shape (n_points, 2) with [R, Z] pairs,
                                            or None to disable.
                @param trim_lcfs When True (default) and X-points are active, LCFS isoflux targets near
                                 the X-point(s) are removed. When False, all defined isoflux (LCFS) points
                                 are used together with the defined X-point(s).

        '''
        if diverted_times is not None and len(diverted_times) != 2:
            raise ValueError('diverted_times must be a (t_start, t_end) tuple.')

        self._diverted_times = diverted_times
        self._x_point_targets = None if x_point_targets is None else np.atleast_2d(x_point_targets)
        self._x_point_weight = x_point_weight
        self._strike_point_targets = None if strike_point_targets is None else np.atleast_2d(strike_point_targets)
        self._trim_lcfs = trim_lcfs
        self._trim_lcfs_perc_limit = trim_lcfs_perc_limit


    # ─── TORAX (TX) Methods ───────────────────────────────────────────────────

    # ── Time-averaging helpers ──────────────────────────────────────────────

    def _get_time_window(self, time, tx_times):
        r'''! Return the (t_start, t_end) averaging window for a given timepoint.
        
                Respects self._t_ave_toggle, self._t_ave_window, self._t_ave_causal,
                and self._t_ave_ignore_start.  If averaging is disabled for this
                timepoint the returned window collapses to (time, time).
        
                @param time     The target timepoint (seconds).
                @param tx_times Sorted numpy array of all available TORAX times.
                @return (t_start, t_end) — inclusive bounds for the averaging window.
                
        '''
        # Check if averaging is active for this timepoint
        if self._t_ave_toggle == 'off':
            return (time, time)
        if self._t_ave_toggle == 'flattop':
            # Inline flattop check: only average when inside detected flat-top
            if not hasattr(self, '_flattop') or not np.any(self._flattop):
                return (time, time)
            ft_times = np.array(self._tm_times)[self._flattop]
            if not (ft_times[0] <= time <= ft_times[-1]):
                return (time, time)
        # 'pulse' → always average

        half = self._t_ave_window / 2.0
        if self._t_ave_causal:
            t_start = time - self._t_ave_window
            t_end   = time
        else:
            t_start = time - half
            t_end   = time + half

        # Clamp to available data range
        t_start = max(t_start, float(tx_times[0]))
        t_end   = min(t_end,   float(tx_times[-1]))

        # Enforce ignore-start: averaging window must not dip below this threshold
        t_earliest = float(tx_times[0]) + self._t_ave_ignore_start
        t_start = max(t_start, t_earliest)

        # If the window collapsed (e.g. very early in the pulse), just use the
        # single requested timepoint so we still return something sensible.
        if t_start > t_end:
            return (time, time)

        return (t_start, t_end)

    # ── Profile extraction ──────────────────────────────────────────────────

    def _interp_tx_profile_onto_psi(self, data_tree, var_name, time, profile_type='linterp'):
        r'''! Interpolate a single TORAX profile snapshot onto self._psi_N.
        
                No averaging, just one timeslice.
                Returns a plain numpy array on self._psi_N.
                
        '''
        var = getattr(data_tree.profiles, var_name)
        var_data = var.sel(time=time, method='nearest').to_numpy()

        # Detect rho coordinate
        if 'rho_cell_norm' in var.coords:
            grid = 'rho_cell_norm'
        elif 'rho_face_norm' in var.coords:
            grid = 'rho_face_norm'
        elif 'rho_norm' in var.coords:
            grid = 'rho_norm'
        else:
            raise ValueError(f"Variable {var_name} does not have a recognized rho coordinate")

        # Psi mapping
        psi_norm_face = data_tree.profiles.psi_norm.sel(time=time, method='nearest').to_numpy()
        psi_rho_norm = data_tree.profiles.psi.sel(time=time, method='nearest').to_numpy()
        psi_norm_rho_norm = (psi_rho_norm - psi_rho_norm[0]) / (psi_rho_norm[-1] - psi_rho_norm[0])
        psi_norm_rho_norm[1] = (psi_norm_face[0] + psi_norm_face[1]) / 2.0

        if grid == 'rho_cell_norm':
            psi_on_grid = psi_norm_rho_norm[1:-1]
        elif grid == 'rho_face_norm':
            psi_on_grid = psi_norm_face
        elif grid == 'rho_norm':
            psi_on_grid = psi_norm_rho_norm

        psi_on_grid_real = psi_on_grid * self._last_surface_factor

        left_fill  = float(var_data[0])
        right_fill = 0.0 if profile_type == 'jphi-linterp' else float(var_data[-1])

        return interp1d(psi_on_grid_real, var_data, kind='linear',
                        fill_value=(left_fill, right_fill),
                        bounds_error=False)(self._psi_N)

    def _extract_tx_profile(self, data_tree, var_name, time, load_into_state='state',
                            normalize=False, profile_type='linterp'):
        r'''! Extract a TORAX profile onto self._psi_N with optional time-averaging.
        
                Replaces the former _pull_tx_onto_psi.  When time-averaging is active
                the profile is interpolated at every TORAX timestep inside the window
                and the results are averaged pointwise on the psi_N grid.
        
                @param data_tree     TORAX output data tree.
                @param var_name      Name of variable (e.g., 'T_i', 'j_ohmic', 'FFprime').
                @param time          Target time value.
                @param load_into_state  'state' → return dict; else return plain array.
                @param normalize     If True, normalize profile by the core value.
                @param profile_type  'linterp' or 'jphi-linterp'.
                
        '''
        tx_times = data_tree.profiles.psi.coords['time'].values
        t_start, t_end = self._get_time_window(time, tx_times)

        if t_start == t_end:
            # No averaging — single snapshot
            data_on_psi = self._interp_tx_profile_onto_psi(data_tree, var_name, time, profile_type)
        else:
            # Collect all TORAX timesteps inside the window
            mask = (tx_times >= t_start) & (tx_times <= t_end)
            win_times = tx_times[mask]
            if len(win_times) == 0:
                data_on_psi = self._interp_tx_profile_onto_psi(data_tree, var_name, time, profile_type)
            else:
                stack = np.stack([
                    self._interp_tx_profile_onto_psi(data_tree, var_name, wt, profile_type)
                    for wt in win_times
                ])
                data_on_psi = np.mean(stack, axis=0)

        # Normalize if requested
        if normalize:
            var = getattr(data_tree.profiles, var_name)
            if 'rho_cell_norm' in var.coords:
                psi_norm_face = data_tree.profiles.psi_norm.sel(time=time, method='nearest').to_numpy()
                psi_rho_norm = data_tree.profiles.psi.sel(time=time, method='nearest').to_numpy()
                psi_norm_rho_norm = (psi_rho_norm - psi_rho_norm[0]) / (psi_rho_norm[-1] - psi_rho_norm[0])
                psi_norm_rho_norm[1] = (psi_norm_face[0] + psi_norm_face[1]) / 2.0
                psi_on_grid_real = psi_norm_rho_norm[1:-1] * self._last_surface_factor
                core_idx = np.argmin(np.abs(self._psi_N - psi_on_grid_real[0]))
                data_on_psi /= data_on_psi[core_idx]
            else:
                data_on_psi /= data_on_psi[0]

        if load_into_state == 'state':
            return {'x': self._psi_N.copy(), 'y': data_on_psi.copy(), 'type': profile_type}
        else:
            return data_on_psi

    # ── Scalar extraction ───────────────────────────────────────────────────

    def _extract_tx_scalar(self, data_tree, var_name, time, source='scalars', scale=1.0):
        r'''! Extract a scalar value from TORAX with optional time-averaging.
        
                @param data_tree  TORAX output data tree.
                @param var_name   Attribute name on data_tree.scalars (or .profiles).
                @param time       Target time (seconds).
                @param source     'scalars' or 'profiles' — which subtree to read from.
                @param scale      Multiplicative factor applied after extraction.
                @return float — the (optionally averaged) scalar value.
                
        '''
        container = getattr(data_tree, source)
        var = getattr(container, var_name)
        tx_times = var.coords['time'].values
        t_start, t_end = self._get_time_window(time, tx_times)

        if t_start == t_end:
            val = float(var.sel(time=time, method='nearest'))
        else:
            mask = (tx_times >= t_start) & (tx_times <= t_end)
            win_times = tx_times[mask]
            if len(win_times) == 0:
                val = float(var.sel(time=time, method='nearest'))
            else:
                val = float(var.sel(time=win_times).mean())

        return val * scale

    def _extract_tx_scalar_at_rho(self, data_tree, var_name, time, rho_val, rho_coord='rho_norm', scale=1.0):
        r'''! Extract a profile value at a specific rho location as a scalar, with time-averaging.
        
                @param data_tree   TORAX output data tree.
                @param var_name    Attribute name on data_tree.profiles.
                @param time        Target time (seconds).
                @param rho_val     Radial coordinate value to select.
                @param rho_coord   Name of the rho coordinate ('rho_norm', 'rho_face_norm', etc.).
                @param scale       Multiplicative factor applied after extraction.
                @return float — the (optionally averaged) scalar value.
                
        '''
        var = getattr(data_tree.profiles, var_name)
        tx_times = var.coords['time'].values
        t_start, t_end = self._get_time_window(time, tx_times)

        if t_start == t_end:
            val = float(var.sel(time=time, **{rho_coord: rho_val}, method='nearest'))
        else:
            mask = (tx_times >= t_start) & (tx_times <= t_end)
            win_times = tx_times[mask]
            if len(win_times) == 0:
                val = float(var.sel(time=time, **{rho_coord: rho_val}, method='nearest'))
            else:
                val = float(var.sel(time=win_times, **{rho_coord: rho_val}, method='nearest').mean())

        return val * scale

    def _extract_tx_scalar_timeseries(self, data_tree, var_name, source='scalars', scale=1.0):
        r'''! Extract a full time-series scalar from TORAX with time-averaging applied at each point.
        
                Returns dict {'x': times_list, 'y': values_array} suitable for self._results.
        
                @param data_tree  TORAX output data tree.
                @param var_name   Attribute name on data_tree.scalars (or .profiles).
                @param source     'scalars' or 'profiles'.
                @param scale      Multiplicative factor applied after extraction.
                
        '''
        container = getattr(data_tree, source)
        var = getattr(container, var_name)
        tx_times = var.coords['time'].values
        raw = var.to_numpy()

        t_start_arr = np.empty_like(tx_times)
        t_end_arr   = np.empty_like(tx_times)
        for idx, t in enumerate(tx_times):
            t_start_arr[idx], t_end_arr[idx] = self._get_time_window(t, tx_times)

        averaged = np.empty(len(tx_times))
        for idx, t in enumerate(tx_times):
            if t_start_arr[idx] == t_end_arr[idx]:
                averaged[idx] = raw[idx]
            else:
                mask = (tx_times >= t_start_arr[idx]) & (tx_times <= t_end_arr[idx])
                averaged[idx] = np.mean(raw[mask])

        return {
            'x': list(tx_times),
            'y': averaged * scale,
        }

    def _extract_tx_scalar_at_rho_timeseries(self, data_tree, var_name, rho_val,
                                               rho_coord='rho_norm', scale=1.0):
        r'''! Extract a profile-at-fixed-rho time-series with time-averaging.
        
                Returns dict {'x': times_list, 'y': values_array}.
                
        '''
        var = getattr(data_tree.profiles, var_name)
        tx_times = var.coords['time'].values
        raw = var.sel(**{rho_coord: rho_val}, method='nearest').to_numpy()

        averaged = np.empty(len(tx_times))
        for idx, t in enumerate(tx_times):
            ts, te = self._get_time_window(t, tx_times)
            if ts == te:
                averaged[idx] = raw[idx]
            else:
                mask = (tx_times >= ts) & (tx_times <= te)
                averaged[idx] = np.mean(raw[mask])

        return {
            'x': list(tx_times),
            'y': averaged * scale,
        }


    def _apply_tx_set_overrides(self, myconfig):
        r'''! Apply user set_*() overrides to a TORAX config dict (in place).
        
                Only applied when the corresponding attribute is not None (i.e. the user
                called the setter explicitly; None means fall through to the loaded/base
                config). Used by both _get_tx_config (coupling loop 0+) and _run_tx_relax
                so each relax simulation sees the same plasma
                conditions as the main sim.
        
                Does NOT touch geometry, numerics.{t_initial, t_final, fixed_dt},
                or profile_conditions.{psi, initial_psi_mode, initial_psi_from_j} —
                those are loop-specific and set by the calling method.
        
                @param myconfig Config dict (modified in place).
                
        '''
        myconfig.setdefault('profile_conditions', {})
        myconfig.setdefault('numerics', {})

        if self._Ip is not None:
            myconfig['profile_conditions']['Ip'] = self._Ip

        if self._n_e is not None:
            myconfig['profile_conditions']['n_e'] = self._n_e

        if self._T_e is not None:
            myconfig['profile_conditions']['T_e'] = self._T_e

        if self._T_i is not None:
            myconfig['profile_conditions']['T_i'] = self._T_i

        if self._Zeff is not None:
            myconfig.setdefault('plasma_composition', {})
            myconfig['plasma_composition']['Z_eff'] = self._Zeff

        if self._main_ion is not None:
            myconfig.setdefault('plasma_composition', {})
            myconfig['plasma_composition']['main_ion'] = self._main_ion

        if self._impurity is not None:
            myconfig.setdefault('plasma_composition', {})
            myconfig['plasma_composition']['impurity'] = self._impurity

        if self._enable_fusion:
            myconfig.setdefault('sources', {})
            myconfig['sources'].setdefault('fusion', {})

        if self._enable_ei_exchange:
            myconfig.setdefault('sources', {})
            myconfig['sources'].setdefault('ei_exchange', {})

        if self._ecrh_loc is not None:
            myconfig.setdefault('sources', {})
            myconfig['sources'].setdefault('ecrh', {})
            myconfig['sources']['ecrh']['P_total'] = self._ecrh_heating
            myconfig['sources']['ecrh']['gaussian_location'] = self._ecrh_loc
            myconfig['sources']['ecrh']['gaussian_width'] = self._ecrh_width

        if self._generic_heat is not None:
            nbi_times, nbi_pow = zip(*self._generic_heat.items())
            myconfig.setdefault('sources', {})
            myconfig['sources'].setdefault('generic_heat', {})
            myconfig['sources']['generic_heat']['P_total'] = (nbi_times, nbi_pow)
            myconfig['sources']['generic_heat']['gaussian_location'] = self._generic_heat_loc
            myconfig['sources']['generic_heat']['gaussian_width'] = self._generic_heat_width
            myconfig['sources']['generic_heat']['absorption_fraction'] = self._generic_heat_absorption_fraction 
            myconfig['sources']['generic_heat']['electron_heat_fraction'] = self._generic_heat_electron_heat_fraction 

            if self._use_nbi_current:
                myconfig['sources'].setdefault('generic_current', {})
                myconfig['sources']['generic_current']['use_absolute_current'] = True
                myconfig['sources']['generic_current']['I_generic'] = (nbi_times, _NBI_W_TO_MA * np.array(nbi_pow))
                myconfig['sources']['generic_current']['gaussian_location'] = self._generic_heat_loc

        if self._generic_particle_location is not None:
            myconfig['sources']['generic_particle'] = {}
            myconfig['sources']['generic_particle']['deposition_location'] = self._generic_particle_location
            myconfig['sources']['generic_particle']['particle_width'] = self._generic_particle_width
            myconfig['sources']['generic_particle']['S_total'] = self._generic_particle_s_total
            
        if self._ped_mode == 'off' and self._pedestal_config is None:
            # No pedestal: leave myconfig['pedestal'] unset (TORAX default no_pedestal).
            pass
        elif self._pedestal_config is not None:
            # Full pedestal dict replacement (set_pedestal(config=...)).
            myconfig['pedestal'] = copy.deepcopy(self._pedestal_config)
        elif self._ped_mode == 'internal_manual':
            if self._manual_ibc_active():
                # Timing is known: TORAX pedestal off, impose the pedestal as an IBC.
                myconfig['pedestal'] = {'model_name': 'no_pedestal'}
                myconfig['profile_conditions']['internal_boundary_conditions'] = self._build_manual_ibc()
                # IBC penalty stiffness (soft target, not a hard BC): higher tracks the imposed
                # pedestal more exactly, lower lets TORAX transport smooth the imposed shape.
                if self._ped_T_source_prefactor is not None:
                    myconfig['numerics']['adaptive_T_source_prefactor'] = self._ped_T_source_prefactor
                if self._ped_n_source_prefactor is not None:
                    myconfig['numerics']['adaptive_n_source_prefactor'] = self._ped_n_source_prefactor
            else:
                # Detection loop: ADAPTIVE_SOURCE + formation model gives the L/H state
                # machine; set_pedestal=False so the (later) IBC solely owns the edge.
                myconfig['pedestal'] = {
                    'model_name': 'set_T_ped_n_ped',
                    'set_pedestal': False,
                    'mode': 'ADAPTIVE_SOURCE',
                    'use_formation_model_with_adaptive_source': True,
                    'transition_time_width': self._ped_transition_time,
                    'formation_model': {'model_name': self._ped_formation_model_name},
                }
        else:
            ped_enabled = (
                self._set_pedestal is not None
                and not (
                    isinstance(self._set_pedestal, (bool, np.bool_))
                    and (not bool(self._set_pedestal))
                )
            )

            if ped_enabled:
                # Build in required key order:
                # model_name -> set_pedestal -> mode -> rho_norm_ped_top -> T_i/T_e/n_e.
                ped_cfg = {}
                ped_cfg['model_name'] = 'set_T_ped_n_ped'
                ped_cfg['set_pedestal'] = self._set_pedestal
                ped_cfg['mode'] = self._ped_mode
                if self._ped_mode == 'ADAPTIVE_SOURCE' and self._ped_formation_model:
                    ped_cfg['use_formation_model_with_adaptive_source'] = True
                if self._ped_top is not None:
                    ped_cfg['rho_norm_ped_top'] = self._ped_top
                if self._T_i_ped is not None:
                    ped_cfg['T_i_ped'] = self._T_i_ped
                if self._T_e_ped is not None:
                    ped_cfg['T_e_ped'] = self._T_e_ped
                if self._n_e_ped is not None:
                    ped_cfg['n_e_ped_is_fGW'] = False
                    ped_cfg['n_e_ped'] = self._n_e_ped
                myconfig['pedestal'] = ped_cfg
            else:
                # Replace pedestal section entirely to avoid invalid extra keys
                # for the no_pedestal model.
                myconfig['pedestal'] = {'model_name': 'no_pedestal'}

        if self._nbar is not None:
            myconfig['profile_conditions']['nbar'] = self._nbar
        if self._normalize_to_nbar is not None:
            myconfig['profile_conditions']['normalize_n_e_to_nbar'] = self._normalize_to_nbar

        if self._ne_right_bc is not None:
            myconfig['profile_conditions']['n_e_right_bc_is_fGW'] = False
            myconfig['profile_conditions']['n_e_right_bc'] = self._ne_right_bc

        if self._Te_right_bc is not None:
            myconfig['profile_conditions']['T_e_right_bc'] = self._Te_right_bc
        if self._Ti_right_bc is not None:
            myconfig['profile_conditions']['T_i_right_bc'] = self._Ti_right_bc

        if self._evolve_density is not None:
            myconfig['numerics']['evolve_density'] = self._evolve_density
        if self._evolve_current is not None:
            myconfig['numerics']['evolve_current'] = self._evolve_current
        if self._evolve_Ti is not None:
            myconfig['numerics']['evolve_ion_heat'] = self._evolve_Ti
        if self._evolve_Te is not None:
            myconfig['numerics']['evolve_electron_heat'] = self._evolve_Te

        if self._gp_s is not None and self._gp_dl is not None:
            myconfig.setdefault('sources', {})
            myconfig['sources']['gas_puff'] = {
                'S_total': self._gp_s,
                'puff_decay_length': self._gp_dl,
            }

        if (
            self._pellet_deposition_location is not None
            and self._pellet_width is not None
            and self._pellet_s_total is not None
        ):
            myconfig.setdefault('sources', {})
            myconfig['sources']['pellet'] = {
                'S_total': self._pellet_s_total,
                'pellet_width': self._pellet_width,
                'pellet_deposition_location': self._pellet_deposition_location,
            }

        myconfig.setdefault('transport', {})
        if self._chi_min is not None:
            myconfig['transport']['chi_min'] = self._chi_min
        if self._chi_max is not None:
            myconfig['transport']['chi_max'] = self._chi_max
        if self._De_min is not None:
            myconfig['transport']['D_e_min'] = self._De_min
        if self._De_max is not None:
            myconfig['transport']['D_e_max'] = self._De_max
        if self._Ve_min is not None:
            myconfig['transport']['V_e_min'] = self._Ve_min
        if self._Ve_max is not None:
            myconfig['transport']['V_e_max'] = self._Ve_max

    def _get_tx_config(self):
        r'''! Generate config object for Torax simulation.
        
                Build order
                -----------
                1. Deep-copy BASE_CONFIG.
                2. Deep-merge the loaded config on top (if load_TORAX_config() was called).
                   Every key in the loaded config overwrites the matching base key;
                   keys only in BASE_CONFIG are kept as-is.
                3. Override geometry (always set to use TokaMaker equilibria).
                   For loop 1+ with injected psi and relax=False, the seed EQDSK is
                   forced at t_initial (psi from initial TORAX relax on seed); if relax=True,
                   loop N relax keeps psi on the TM i=0 surface so TM EQDSK is used.
                4. Override t_initial / t_final / fixed_dt from __init__ params.
                5. Use psi profile from initial TORAX relax (if available) from profile_conditions.
                6. Apply any explicit set_*() overrides (only when the attribute is
                   not None, i.e. the user called the setter after load_TORAX_config).
                7. If relax_kinetics was True for the last relax, override n_e, T_e, T_i
                   with profiles taken from the end of that relax simulation.
                8. If fly(..., steady_state_mode=True) and this is not the first loop,
                   override psi and kinetic profiles with profiles saved from the previous main TORAX run at
                   t_final (see _capture_steady_state_tx_seed).
        
                @return Torax config object.
                
        '''

        # ── 1. Start from base config ──────────────────────────────────────
        myconfig = copy.deepcopy(BASE_CONFIG)

        # ── 2. Deep-merge loaded config ────────────────────────────────────
        if self._loaded_config is not None:
            self._tx_config_merge(myconfig, self._loaded_config)

        # ── 3. Geometry (always set by TokaMaker_TORAX) ─────────────────────────────
        myconfig['geometry'] = {
            'geometry_type': 'eqdsk',
            'geometry_directory': os.getcwd(),
            'last_surface_factor': self._last_surface_factor,
            'n_surfaces': 50,
            'Ip_from_parameters': True, # True tells TX to pull from config, not from eqdsk, in case eqdsks fail TX retains correct Ip targets
        }
        if self._coupling_iteration_is_first():
            eq_safe = []
            t_safe = []
            t_skipped = []
            for i, t in enumerate(self._eqtimes):
                eq = self._init_files[i]
                if self._test_eqdsk_tx_config(eq):
                    eq_safe.append(eq)
                    t_safe.append(t)
                else:
                    if not self._skip_bad_init_eqdsks:
                        raise ValueError(f'Bad initial gEQDSK at t={t}: {eq}')
                    t_skipped.append(t)
            self._log(f'\tTX: {len(eq_safe)}/{len(self._eqtimes)} initial EQDSKs valid'
                      + (f', skipped {len(t_skipped)}' if t_skipped else ''))
            myconfig['geometry']['geometry_configs'] = {
                t: {'geometry_file': eq_safe[i], 'cocos': self._cocos} for i, t in enumerate(t_safe)
            }
        else:
            # For times where TM succeeded last loop, use the TM-solved EQDSK.
            # For times where TM failed, omit from the map (TORAX interpolates from neighbors),
            # except t=0 which always gets a seed fallback (see below).
            full_eqdsk_map = {}
            n_tm = 0
            if self._steady_state_mode:
                i_last = len(self._tm_times) - 1
                eq_last = os.path.join(self._eqdsk_dir, f'{self._current_loop - 1:03d}.{i_last:03d}.eqdsk')
                tm_last_ok = eq_last not in self._eqdsk_skip and os.path.isfile(eq_last)
                if tm_last_ok:
                    for t in self._tm_times:
                        full_eqdsk_map[t] = eq_last
                    n_tm = len(self._tm_times)
                    self._log(
                        f'Loop {self._current_loop}: steady_state_mode geometry — all TX times use '
                        f'final TM EQDSK from loop {self._current_loop - 1}: {os.path.basename(eq_last)}.'
                    )
                else:
                    self._log(
                        f'Loop {self._current_loop}: steady_state_mode: final TM EQDSK missing or skipped '
                        f'({os.path.basename(eq_last)}); falling back to per-time TM map.'
                    )

            if not self._steady_state_mode or n_tm == 0:
                # Build a geometry entry for EVERY tm_time so a failed TM solve can
                # never leave a gap (a gap makes TORAX silently interpolate geometry
                # across the failure, which can corrupt a/R/B_0 over the whole time
                # series — e.g. a failed t=0 once flattened the entire ramp). Fallback
                # priority per failed time: (1) the most-recent EARLIER tm_time that
                # solved this loop (carry-forward last-good equilibrium), else
                # (2) the nearest seed EQDSK (by eqtime), else (3) the next later
                # solved tm_time. Each fallback is TORAX-validated before use.
                full_eqdsk_map = {}
                n_tm = 0
                solved = {}  # i -> eqdsk path for TM solves that succeeded
                for i, t in enumerate(self._tm_times):
                    eqdsk = os.path.join(self._eqdsk_dir, f'{self._current_loop - 1:03d}.{i:03d}.eqdsk')
                    # TM stage already saved and validated each EQDSK with TORAX (_run_tm).
                    if eqdsk not in self._eqdsk_skip and os.path.isfile(eqdsk):
                        solved[i] = eqdsk

                def _nearest_seed(t):
                    j = int(np.argmin(np.abs(np.asarray(self._eqtimes, float) - t)))
                    return self._init_files[j]

                last_good = None  # most-recent solved eqdsk (carry-forward)
                n_fallback = 0
                for i, t in enumerate(self._tm_times):
                    if i in solved:
                        full_eqdsk_map[t] = solved[i]
                        last_good = solved[i]
                        n_tm += 1
                        continue
                    # Failed solve at this time — pick a valid fallback in priority order.
                    candidates = []
                    if last_good is not None:
                        candidates.append(last_good)                 # (1) last-good
                    candidates.append(_nearest_seed(t))              # (2) nearest seed
                    nxt = next((solved[k] for k in sorted(solved) if k > i), None)
                    if nxt is not None:
                        candidates.append(nxt)                       # (3) next solved
                    chosen = next((c for c in candidates
                                   if self._test_eqdsk_tx_config(c)), None)
                    if chosen is not None:
                        full_eqdsk_map[t] = chosen
                        n_fallback += 1
                    else:
                        self._log(f'Warning: Loop {self._current_loop}: TM failed at '
                                  f't={t:.2f} and no valid fallback EQDSK found; '
                                  f'leaving gap (TORAX will interpolate).')
                if n_fallback:
                    self._log(f'Loop {self._current_loop}: {n_fallback} failed TM '
                              f'timestep(s) filled by last-good/seed fallback.')

            # Injected psi from initial relax used the seed EQDSK.  If we used
            # the TM EQDSK at t_init without a loop N re-relax, the metric changes
            # and j becomes inconsistent.  When relax=False, force seed at t_init.
            # When relax=True, loop N relax aligns psi with TM i=0 — keep TM.
            # steady_state_mode uses the final TM EQDSK everywhere and seeds psi from the
            # previous TORAX t_final; do not replace t_init with the seed file.
            if (self._psi_init is not None and not self._relax
                    and not self._steady_state_mode):
                seed_eqdsk_tinit = self._init_files[0]
                if self._test_eqdsk_tx_config(seed_eqdsk_tinit):
                    full_eqdsk_map[self._t_init] = seed_eqdsk_tinit
                    self._log(
                        f'Loop {self._current_loop}: seed EQDSK at t_init={self._t_init} s for TORAX '
                        f'(psi from initial relax on seed; inter-loop relax disabled).'
                    )

            if n_tm == 0:
                self._log(f'Warning: Loop {self._current_loop}: no valid TM EQDSKs from loop {self._current_loop-1}, using all seed EQDSKs.')
            else:
                self._log(f'Loop {self._current_loop}: using {n_tm}/{len(self._tm_times)} TM-solved EQDSKs, {len(self._tm_times)-n_tm} seed fallbacks.')
            
            myconfig['geometry']['geometry_configs'] = {
                t: {'geometry_file': eqdsk_f, 'cocos': self._cocos} for t, eqdsk_f in full_eqdsk_map.items()
            }

        if self._tx_grid_type == 'n_rho':
            myconfig['geometry']['n_rho'] = self._tx_grid
        elif self._tx_grid_type == 'face_centers':
            myconfig['geometry']['face_centers'] = self._tx_grid

        # ── 4. Override t_initial / t_final / fixed_dt from __init__ ───────
        myconfig.setdefault('numerics', {})
        myconfig['numerics']['t_initial'] = self._t_init
        myconfig['numerics']['t_final'] = self._t_final
        myconfig['numerics']['fixed_dt'] = self._tx_dt

        # ── 5. Psi profile from last TORAX relax (initial / inter-loop) ───────────────────
        myconfig.setdefault('profile_conditions', {})
        if self._psi_init is not None:
            myconfig['profile_conditions']['psi'] = self._psi_init
            myconfig['profile_conditions']['initial_psi_mode'] = 'profile_conditions'
            myconfig['profile_conditions']['initial_psi_from_j'] = False
        else:
            myconfig['profile_conditions']['initial_psi_mode'] = 'geometry'  # if initial relax was skipped

        # ── 6. Explicit set_*() overrides ──────────────────────────────────
        self._apply_tx_set_overrides(myconfig)

        # ── 7. Kinetic profiles from last relax (when relax_kinetics was True) ──
        # Applied after set_*() so relaxed n_e / T_e / T_i replace the initial
        # slice of user schedules, matching how psi from relax seeds the flux.
        if self._n_e_init is not None:
            myconfig['profile_conditions']['n_e'] = self._n_e_init
        if self._T_e_init is not None:
            myconfig['profile_conditions']['T_e'] = self._T_e_init
        if self._T_i_init is not None:
            myconfig['profile_conditions']['T_i'] = self._T_i_init

        # ── 8. Steady-state coupling: previous loop TORAX t_final → IC for this loop ──
        # If an inter-loop relax already ran this loop, its output should define the
        # main-run ICs (psi always, and kinetics when relax_kinetics=True). In that case
        # do not overwrite with the raw steady-state seed here.
        if (self._steady_state_mode
                and not self._coupling_iteration_is_first()
                and not (self._relax and self._current_loop >= 1)
                and self._steady_state_tx_seed):
            seed = self._steady_state_tx_seed
            pc = myconfig['profile_conditions']
            for k in ('psi', 'n_e', 'T_e', 'T_i'):
                if k in seed and seed[k] is not None:
                    pc[k] = copy.deepcopy(seed[k])
            pc['initial_psi_mode'] = 'profile_conditions'
            pc['initial_psi_from_j'] = False
            self._log(
                f'Loop {self._current_loop}: steady_state_mode: TORAX initial profiles from '
                f'previous loop t_final={self._t_final:g} s.'
            )

        if self._output_mode in ('normal', 'debug') and self._out_dir is not None:
            cfg_name = f'tx_config_loop{self._current_loop:03d}.py'
            if self._output_file_tag is not None:
                cfg_name = f'{self._output_file_tag}_{cfg_name}'
            config_filename = os.path.join(self._out_dir, cfg_name)
            with open(config_filename, 'w') as f:
                f.write('# Torax configuration\n')
                f.write(f'# Loop {self._current_loop}\n\n')
                f.write('tx_config = ')
                f.write(pprint.pformat(self._numpy_to_plain_python(myconfig), width=100, sort_dicts=False))

        tx_config = torax.ToraxConfig.from_dict(myconfig)
        return tx_config

    def _test_eqdsk_tx_config(self, eqdsk, *, quiet=False):
        r'''! Return whether TORAX accepts eqdsk as ToraxConfig geometry.
        
                @param quiet If True, do not print or append to the coupling log on failure
                       (used for intermediate-resolution retries in _run_tm). If False,
                       failure is only reported when output_mode is 'debug'.
                
        '''
        myconfig = copy.deepcopy(BASE_CONFIG)
        if self._loaded_config is not None:
            self._tx_config_merge(myconfig, self._loaded_config)
        myconfig['geometry'] = {
            'geometry_type': 'eqdsk',
            'geometry_directory': os.getcwd(),
            'last_surface_factor': self._last_surface_factor,
            'Ip_from_parameters': False,
            'geometry_file': eqdsk,
            'cocos': self._cocos,
        }
        try:
            _ = torax.ToraxConfig.from_dict(myconfig)
            return True
        except Exception as e:
            if not quiet and self._output_mode == 'debug':
                self._log(f"TEST EQDSK FAILED: {eqdsk} — {repr(e)}")
                self._print(f'    EQDSK rejected by TORAX: {os.path.basename(eqdsk)}')
            return False

    def _capture_relax_tx_profiles_from_datatree(self, data_tree, time_val=None):
        r'''! Store psi, n_e, T_e, T_i at time_val for debug relax figures / history.
        
                Updates _relax_profiles_snapshot temporarily so the caller can append a copy
                to _relax_mainrun_profile_history. Not used to seed inter-loop relax runs
                (those profiles come from the user config each time).
                
        '''
        if time_val is None:
            time_val = self._t_init
        snap = {}
        _sel = dict(time=time_val, method='nearest')
        for _name in ('psi', 'n_e', 'T_e', 'T_i'):
            _xr = getattr(data_tree.profiles, _name).sel(**_sel)
            _rho = _xr.coords['rho_norm'].to_numpy()
            _arr = _xr.to_numpy()
            snap[_name] = ([self._t_init], _rho.tolist(), [_arr.tolist()])
        self._relax_profiles_snapshot = snap

    def _capture_steady_state_tx_seed(self, data_tree):
        r'''! Store psi, n_e, T_e, T_i at t_final for steady_state_mode next coupling loop.
        
                Profile tuples use t_init as the time key (TORAX profile_conditions convention),
                matching _psi_init / relax snapshots.
                
        '''
        t_fin = float(self._t_final)
        _sel = dict(time=t_fin, method='nearest')
        seed = {}
        for _name in ('psi', 'n_e', 'T_e', 'T_i'):
            _xr = getattr(data_tree.profiles, _name).sel(**_sel)
            _rho = _xr.coords['rho_norm'].to_numpy()
            _arr = _xr.to_numpy()
            seed[_name] = ([self._t_init], _rho.tolist(), [_arr.tolist()])
        self._steady_state_tx_seed = seed

    def _manual_ibc_active(self):
        r'''! True when the internal_manual pedestal IBC should be imposed this loop, i.e. the
                L/H timing is known: timing='input' (always), or timing='detect' once the
                detection loop (loop 1) has run and found an L->H time.'''
        if self._ped_mode != 'internal_manual':
            return False
        if self._ped_timing == 'input':
            return True
        return self._current_loop > 1 and self._lh_time is not None

    def _build_manual_ibc(self):
        r'''! Build the internal_boundary_conditions block for T_e, T_i, n_e from the MANUAL
                pedestal knobs. The pedestal foot (rho=1) is the edge BC from set_ne/set_Te/set_Ti
                (so n_e_right_bc and the IBC agree there); the inner-edge gradient uses the
                smoothed slope from _smooth_ibc_inner_edge when available, else a
                straight line from ped_top down to foot.

                If a core_height is set for a field, a CORE band (rho 0->just inside ped_rho) is
                merged in that ramps the evolved L-mode core into a smooth H-mode-like core shape
                over the L->H transition, then retreats. The two bands are disjoint (the pedestal
                band owns [ped_rho,1] at all times); they meet at ped_rho with the same value.

                @return {field: SparseTimeVaryingArray} for the fields with a height set.
        '''
        if self._lh_time is None:
            return {}   # no L->H transition: no pedestal to impose
        windows = [(self._lh_time, self._hl_time)]
        ped_rho = 1.0 - self._ped_width
        lh = self._lh_time
        transition = self._ped_transition_time
        # Shape params scaled to ped_width (match the notebook's 0.93/0.025/0.015 at width 0.15):
        # the steep tanh edge is a narrow band near the separatrix, so the inner edge sits on the
        # flat core line (value ~ ped_top) and the foot saturates to the edge BC.
        w = self._ped_width
        knee = 1.0 - 0.47 * w
        hwid = 0.17 * w
        blend = 0.10 * w
        # Core transition-band grid: ends ONE cell short of ped_rho so the band location
        # (0, core_hi) is strictly disjoint from the pedestal band (ped_rho, 1) -- TORAX rejects
        # bands that even touch at an endpoint. The tiny (core_hi, ped_rho) gap is bridged by
        # interpolation; both sides equal ped_height there, so it's seamless.
        core_dx = ped_rho / PED_N_SAMPLE
        rhos_core = np.linspace(0.0, ped_rho - core_dx, PED_N_SAMPLE)
        ibc = {}
        # Snapshot of the shape params actually used to build this loop's IBC, for the
        # debug figure (IBC_debug): per field the ped_top/slope and whether
        # they came from the previous loop's smoothing or the user-height fallback.
        # IBC penalty stiffness actually in effect (None -> TORAX default), for the figure label.
        T_pref = self._ped_T_source_prefactor if self._ped_T_source_prefactor is not None else 2.0e10
        n_pref = self._ped_n_source_prefactor if self._ped_n_source_prefactor is not None else 2.0e8
        self._ped_ibc_snapshot = {'ped_rho': ped_rho, 'knee': knee, 'hwid': hwid,
                                   'blend': blend, 'transition': self._ped_transition_time,
                                   'windows': windows, 'T_pref': T_pref, 'n_pref': n_pref,
                                   'fields': {}}
        for field, height, right_bc, core_height in (
            ('T_e', self._ped_height_Te, self._Te_right_bc, self._core_height_Te),
            ('T_i', self._ped_height_Ti, self._Ti_right_bc, self._core_height_Ti),
            ('n_e', self._ped_height_ne, self._ne_right_bc, self._core_height_ne),
        ):
            if height is None:
                continue
            # Per-tm-time inner-edge (ped_top, slope) from the previous loop's evolved core
            # (_smooth_ibc_inner_edge); empty on the first IBC loop -> default_shape (user height
            # with a flat top, slope=0) is used at all times. Matching the evolved value+slope per
            # time removes the inner-edge kink (which would spike p' and fail TM).
            # Loop-1 default = flat core line at ped_top: the value drops to the foot through the
            # tanh edge band only (per the mtanh design above), a gentler seed than a full-width
            # linear ramp; loops 2+ replace it with the evolved per-time slope.
            shape_by_time = self._ped_smoothed.get(field, {})
            default_shape = (height, 0.0)
            # Foot (rho=1) tracks the edge BC; _pedestal_ibc evaluates it per ON time so the
            # IBC always agrees with right_bc at rho=1.
            spec = _pedestal_ibc(ped_rho, right_bc, shape_by_time, default_shape,
                                 self._ped_transition_time, windows, knee, hwid, blend)
            snap = {'shape_by_time': shape_by_time, 'default_shape': default_shape,
                    'right_bc': right_bc, 'smoothed': bool(shape_by_time), 'user_height': height,
                    'transition_active': False}

            # ── Core (0 -> ped_rho) L->H transition band (optional, per field via core_height) ──
            if core_height is not None:
                # Core target: user value on loop 1; the relaxed evolved on-axis value on loops 2+.
                core = self._ped_core_relaxed.get(field, core_height)
                # Smooth core shape over [0, ped_rho] that flattens INTO the pedestal top (no edge
                # roll-off -> no second pedestal). Ends at ped_height, so it lines up with the
                # pedestal band at ped_rho. exp_a/exp_b set the core curvature (>=2 = flat ends).
                hmode = _core_transition_profile(rhos_core, ped_rho, height, core,
                                                 self._core_exp_a, self._core_exp_b)
                hmode_target = {float(r): float(v) for r, v in zip(rhos_core, hmode)}
                # L-mode start of the ramp = the evolved core near lh (captured last loop); if not
                # available yet, fall back to a flat line at the pedestal top (degenerate but safe).
                lmode = self._ped_lmode_profile.get(field)
                if lmode is None:
                    lmode = {float(r): height for r in rhos_core}
                trans_loc = (0.0, float(rhos_core[-1]))   # strictly below ped_rho (no overlap)
                trans_spec = _transition_ibc(trans_loc, lmode, hmode_target, lh, transition)
                # Merge: the core band (0, core_hi) and the pedestal band (ped_rho, 1) are keyed by
                # different, disjoint (lo,hi) locations, so combine the per-time dicts (union).
                spec = _merge_ibc_specs(spec, trans_spec)
                snap.update({'transition_active': True, 'core': core,
                             'lmode_profile': lmode, 'hmode_target': hmode_target,
                             'trans_loc': trans_loc,
                             'core_source': 'relaxed' if field in self._ped_core_relaxed else 'user'})

            ibc[field] = spec
            self._ped_ibc_snapshot['fields'][field] = snap
        return ibc

    def _smooth_ibc_inner_edge(self):
        r'''! Match the imposed pedestal inner-edge value AND gradient to the TORAX-evolved core
                PER TM-TIME, so the next loop's IBC tracks the evolving profile time-point by
                time-point and joins the core smoothly (no kink/jog -> p' spike -> TM GS failure).
                Reads the evolved profile from self._data_tree (rho_norm grid): at each H-mode
                tm_time the inner-edge value is the core value just inside the edge and the gradient
                is a least-squares fit over a small core-side window (multi-cell fit is less noisy
                than a single-cell diff, esp. for density). Stored as {field: {tm_time: (value, slope)}}.

                NOT time-averaged: the few transition-region tm_times have anomalous steep slopes
                (density still filling in) that would poison a single average and flip its sign.
                L-mode tm_times are skipped (no IBC there). The fitted slope is clamped to <= 0
                (rho increases outward, so a physical pedestal falls from core to foot => slope<0):
                a positive slope would mean the imposed value rises above the pedestal top toward
                the edge, which is non-physical; flooring positives to 0 keeps the core->pedestal
                connection flat-or-falling.
        '''
        dt = getattr(self, '_data_tree', None)
        if dt is None:
            return
        # Ensure the H-mode mask is set (timing='input' has no detection state machine).
        self._set_confinement_from_window()
        ped_rho = 1.0 - self._ped_width
        # Fit window: core-side cells within PED_GRAD_WINDOW of the inner edge (min 3 cells).
        lo_rho = ped_rho - PED_GRAD_WINDOW
        heights = {'T_e': self._ped_height_Te, 'T_i': self._ped_height_Ti, 'n_e': self._ped_height_ne}
        for field, height in heights.items():
            if height is None:
                continue
            da = getattr(dt.profiles, field)
            rho = da.coords['rho_norm'].values
            hi = int(np.argmin(np.abs(rho - ped_rho)))           # IBC inner-edge cell
            lo = min(int(np.searchsorted(rho, lo_rho)), hi - 2)  # >= 3 cells in the fit
            if lo < 0:
                continue
            by_time = {}
            for i, t in enumerate(self._tm_times):
                if not self._state['confinement_mode'][i]:
                    continue   # skip pure-L-mode points (no IBC there)
                y = da.sel(time=t, method='nearest').values
                slope = float(np.polyfit(rho[lo:hi + 1], y[lo:hi + 1], 1)[0])

                # Clamp to <= 0: rho increases outward, so a physical pedestal falls from the core
                # to the foot (slope < 0). Floor any positive fitted slope to 0 so the imposed value
                # never rises above the pedestal top toward the edge.
                slope = min(slope, 0.0)
                by_time[float(t)] = (float(y[hi]), slope)
            if by_time:
                self._ped_smoothed[field] = by_time

            # Core transition: read the relaxed on-axis (rho~0) value after the transition band has
            # retreated and the core has relaxed, for the next loop's core target. dwell = one fit
            # window past lh+transition.
            core_height = {'T_e': self._core_height_Te, 'T_i': self._core_height_Ti,
                           'n_e': self._core_height_ne}[field]
            if core_height is not None and self._lh_time is not None:
                t_read = self._lh_time + self._ped_transition_time + PED_GRAD_WINDOW
                axis = int(np.argmin(rho))   # rho ~ 0 cell
                y_read = da.sel(time=t_read, method='nearest').values
                self._ped_core_relaxed[field] = float(y_read[axis])
                # Evolved core at ~lh = the start of next loop's transition ramp, on the same core
                # grid as _build_manual_ibc (ends one cell short of ped_rho, disjoint from pedestal).
                core_dx = ped_rho / PED_N_SAMPLE
                rhos_core = np.linspace(0.0, ped_rho - core_dx, PED_N_SAMPLE)
                y_lh = da.sel(time=self._lh_time, method='nearest').values
                prof = np.interp(rhos_core, rho, y_lh)
                self._ped_lmode_profile[field] = {float(r): float(v) for r, v in zip(rhos_core, prof)}

    def _run_torax(self, config):
        r'''! Run TORAX via prepare_simulation + run_loop so we keep the raw list[SimState]
                (torax.run_simulation() drops the pedestal confinement_mode). Returns the same
                (data_tree, history) that run_simulation() would, plus the SimState list.

                @param config TORAX config as a dict or a ToraxConfig.
                @return Tuple (data_tree, state_history, state_list).
        '''
        tx_config = config if isinstance(config, torax.ToraxConfig) else torax.ToraxConfig.from_dict(config)
        initial_state, post_processed, step_fn = run_sim_lib.prepare_simulation(tx_config)
        state_list, pp_history, sim_error = run_loop_lib.run_loop(
            initial_state=initial_state,
            initial_post_processed_outputs=post_processed,
            step_fn=step_fn,
            progress_bar=True,
        )
        history = output_lib.StateHistory(
            state_history=state_list,
            post_processed_outputs_history=pp_history,
            sim_error=sim_error,
            torax_config=tx_config,
        )
        return history.simulation_output_to_xr(), history, state_list

    def _capture_confinement_mode(self, state_list):
        r'''! Map per-step pedestal confinement_mode to an H-mode bool series, store it per
                tm_time in _state['confinement_mode'], and record the H-mode start/end times.
                H_MODE and TRANSITIONING_TO_H_MODE -> True; L_MODE and TRANSITIONING_TO_L_MODE
                -> False. No-op (leaves all False) if the L/H state machine is inactive.

                The confinement_mode can dither around a transition, so a transition is only
                accepted once the new state persists for at least the pedestal transition_time
                (dwell debounce): a brief L dip inside H-mode is not read as the back-transition,
                and a brief H blip in L-mode is not read as the L->H transition.

                @param state_list list[SimState] from _run_torax.
        '''
        h_codes = {int(ConfinementMode.H_MODE), int(ConfinementMode.TRANSITIONING_TO_H_MODE)}
        tx_times, is_h = [], []
        for s in state_list:
            pts = s.pedestal_transition_state
            if pts is None:
                continue
            tx_times.append(float(s.t))
            is_h.append(int(pts.confinement_mode) in h_codes)
        if not tx_times:
            return
        tx_times = np.array(tx_times)
        is_h = np.array(is_h, dtype=bool)
        dwell = self._ped_transition_time

        # First L->H that stays H for at least `dwell` (ignores brief H blips in L-mode).
        lh = self._first_sustained(tx_times, is_h, True, 0, dwell)
        if lh is None:
            return
        self._lh_time = float(tx_times[lh])
        # First H->L after LH that stays L for at least `dwell` (ignores brief L dips in H-mode);
        # if none, H-mode persists to the end and there is no back-transition.
        hl = self._first_sustained(tx_times, is_h, False, lh, dwell)
        self._hl_time = None if hl is None else float(tx_times[hl])

        self._set_confinement_from_window(end_fallback=float(tx_times[-1]))

        hl_str = f'{self._hl_time:.4g} s' if self._hl_time is not None else 'none (H-mode to t_final)'
        self._log(f'Loop {self._current_loop}: L->H at {self._lh_time:.4g} s, H->L at {hl_str}.')
        self._print(f'  TORAX: L->H transition at {self._lh_time:.4g} s, H->L back-transition at {hl_str}')

    def _set_confinement_from_window(self, end_fallback=None):
        r'''! Fill the per-tm-time H-mode flag from the known [lh, hl] window. Used by both the
                detection path and timing='input' (where there is no detection state machine, so
                the flag would otherwise stay all-False and the IBC smoothing would skip every
                point). No-op if no L->H time is known.'''
        if self._lh_time is None:
            return
        hl_t = self._hl_time
        if hl_t is None:
            hl_t = end_fallback if end_fallback is not None else self._t_final
        for i, t in enumerate(self._tm_times):
            self._state['confinement_mode'][i] = self._lh_time <= t <= hl_t

    @staticmethod
    def _first_sustained(times, flags, target, start, dwell):
        r'''! Index of the first step at/after `start` where `flags == target` and that value
                holds continuously for at least `dwell` seconds. Returns None if never sustained.'''
        n = len(flags)
        i = start
        while i < n:
            if flags[i] == target:
                j = i
                while j + 1 < n and flags[j + 1] == target:
                    j += 1
                if times[j] - times[i] >= dwell or j == n - 1:
                    return i
                i = j + 1   # too brief: skip this run and keep looking
            else:
                i += 1
        return None

    def _run_tx_relax(self, *, stage, eqdsk_path, prescribed_profiles):
        r'''! Short TORAX relax: initial run on the seed EQDSK, or inter-loop on TM i=0 EQDSK.
        
                Uses flattened user inputs (_apply_tx_set_overrides + _flatten_time_dependent).
                If prescribed_profiles is None, psi follows EQDSK initial_psi_mode='geometry'
                and n_e, T_e, T_i stay as already merged from base / loaded config and
                _apply_tx_set_overrides (user inputs). If a dict is passed, it must supply
                psi, n_e, T_e, T_i tuples (advanced; loop N relax uses None so kinetics are always user-specified).
        
                @param stage 'initial' or 'interloop' (logging / output names only).
                @param eqdsk_path Path to gEQDSK for geometry_configs at t_initial.
                @param prescribed_profiles None, or dict with keys psi, n_e, T_e, T_i (3-tuples).
                
        '''
        runtime = float(self._relax_duration)
        dt_relax = float(self._relax_dt)
        tag = 'Initial TORAX relax' if stage == 'initial' else f'Loop {self._current_loop} inter-loop relax'
        self._log(f'{tag}: building config ({runtime:.4g} s, dt={dt_relax})...')
        with self._loop0_coarse_tx_relax_scope(stage):

            self._n_e_init = None
            self._T_e_init = None
            self._T_i_init = None

            init_config = copy.deepcopy(BASE_CONFIG)
            if self._loaded_config is not None:
                self._tx_config_merge(init_config, self._loaded_config)

            self._apply_tx_set_overrides(init_config)

            use_path = eqdsk_path

            init_config['geometry'] = {
                'geometry_type': 'eqdsk',
                'geometry_directory': os.getcwd(),
                'last_surface_factor': self._last_surface_factor,
                'n_surfaces': 50,
                'Ip_from_parameters': True,
                'geometry_configs': {self._t_init: {'geometry_file': use_path, 'cocos': self._cocos}},
            }

            if self._tx_grid_type == 'n_rho':
                init_config['geometry']['n_rho'] = self._tx_grid
            elif self._tx_grid_type == 'face_centers':
                init_config['geometry']['face_centers'] = self._tx_grid

            init_config.setdefault('numerics', {})
            init_config['numerics']['t_initial'] = self._t_init
            init_config['numerics']['t_final'] = self._t_init + runtime
            init_config['numerics']['fixed_dt'] = dt_relax

            if not getattr(self, '_relax_kinetics', False):
                init_config['numerics']['evolve_current'] = True
                init_config['numerics']['evolve_density'] = False
                init_config['numerics']['evolve_ion_heat'] = False
                init_config['numerics']['evolve_electron_heat'] = False
            else:
                self._log(f'{tag}: relax_kinetics ON (evolve_* from set_evolve / config).')
                init_config['numerics'].setdefault('evolve_current', True)

            init_config.setdefault('profile_conditions', {})
            if prescribed_profiles is None:
                init_config['profile_conditions'].pop('psi', None)
                init_config['profile_conditions']['initial_psi_mode'] = 'geometry'
                init_config['profile_conditions']['initial_psi_from_j'] = False
            else:
                pc = init_config['profile_conditions']
                pc['psi'] = copy.deepcopy(prescribed_profiles['psi'])
                pc['initial_psi_mode'] = 'profile_conditions'
                pc['initial_psi_from_j'] = False
                pc['n_e'] = copy.deepcopy(prescribed_profiles['n_e'])
                pc['T_e'] = copy.deepcopy(prescribed_profiles['T_e'])
                pc['T_i'] = copy.deepcopy(prescribed_profiles['T_i'])

            self._flatten_time_dependent(init_config)

            pc_flat = init_config.get('profile_conditions', {})
            _user_ref_curves = {}
            for _k in ('n_e', 'T_e', 'T_i'):
                if _k in pc_flat:
                    _r, _y = TokaMaker_TORAX._relax_flat_profile_to_rho_y(pc_flat[_k])
                    if _r is not None:
                        _user_ref_curves[_k] = (_r, _y)
            if 'psi' in pc_flat:
                _r, _y = TokaMaker_TORAX._relax_flat_profile_to_rho_y(pc_flat['psi'])
                if _r is not None:
                    _user_ref_curves['psi'] = (_r, _y)

            if self._output_mode in ('normal', 'debug') and self._out_dir is not None:
                if stage == 'initial':
                    cfg_name = 'tx_config_relax000_initial.py'
                else:
                    cfg_name = f'tx_config_relax_inter_{self._current_loop:03d}.py'
                if self._output_file_tag is not None:
                    cfg_name = f'{self._output_file_tag}_{cfg_name}'
                config_filename = os.path.join(self._out_dir, cfg_name)
                with open(config_filename, 'w') as f:
                    f.write('# Torax configuration\n# TORAX relax ({stage})\n\n'.format(stage=stage))
                    f.write('tx_config = ')
                    f.write(pprint.pformat(self._numpy_to_plain_python(init_config), width=100, sort_dicts=False))

            self._log(f'{tag}: running TORAX...')
            data_tree, hist, _ = self._run_torax(init_config)

            if hist.sim_error != torax.SimError.NO_ERROR:
                raise ValueError(f'TORAX relax ({stage}) failed: {hist.sim_error}')

            t_final_relax = self._t_init + runtime

            psi_xr = data_tree.profiles.psi.sel(time=t_final_relax, method='nearest')
            rho_psi_arr = psi_xr.coords['rho_norm'].to_numpy()
            psi_arr = psi_xr.to_numpy()
            self._psi_init = ([self._t_init], rho_psi_arr.tolist(), [psi_arr.tolist()])

            if getattr(self, '_relax_kinetics', False):
                for _name, _attr in (('n_e', '_n_e_init'), ('T_e', '_T_e_init'), ('T_i', '_T_i_init')):
                    _xr = getattr(data_tree.profiles, _name).sel(time=t_final_relax, method='nearest')
                    _rho = _xr.coords['rho_norm'].to_numpy()
                    _arr = _xr.to_numpy()
                    setattr(self, _attr, ([self._t_init], _rho.tolist(), [_arr.tolist()]))
                self._log(f'{tag}: relaxed n_e, T_e, T_i will override coupling configs after set_*().')
            elif prescribed_profiles is not None:
                # Kinetics were held fixed during relax; carry those fixed inputs into the
                # main coupling config so steady-state inter-loop relax remains the final IC step.
                for _name, _attr in (('n_e', '_n_e_init'), ('T_e', '_T_e_init'), ('T_i', '_T_i_init')):
                    _v = prescribed_profiles.get(_name)
                    setattr(self, _attr, copy.deepcopy(_v) if _v is not None else None)

            if getattr(self, '_debug_mode', False) and self._out_dir is not None:
                if stage == 'initial':
                    _fig_name = 'tx_relax_profiles_initial.png'
                else:
                    _fig_name = f'tx_relax_profiles_inter_loop{self._current_loop:03d}.png'
                if self._output_file_tag is not None:
                    _fig_name = f'{self._output_file_tag}_{_fig_name}'
                _relax_plot_path = os.path.join(self._out_dir, _fig_name)
                try:
                    plot_tx_relax_profiles(
                        self,
                        data_tree,
                        stage=stage,
                        t_final_relax=t_final_relax,
                        user_ref_curves=_user_ref_curves,
                        save_path=_relax_plot_path,
                        display=False,
                    )
                    self._log(f'TORAX relax profile figure ({stage}): {_relax_plot_path}')
                except Exception as _e:
                    self._log(f'TORAX relax profile figure failed ({stage}): {_e}')


    def _run_tx(self):
        r'''! Run the TORAX transport simulation.
                @return Tuple (consumed_flux, consumed_flux_integral).
                
        '''

        if (self._current_loop >= 1 and self._relax
                and not self._coupling_iteration_is_first()):
            prev_lp = self._current_loop - 1
            i_eq = (len(self._tm_times) - 1) if self._steady_state_mode else 0
            tm_eq0 = os.path.join(self._eqdsk_dir, f'{prev_lp:03d}.{i_eq:03d}.eqdsk')
            relax_eq = tm_eq0
            if tm_eq0 in self._eqdsk_skip or not os.path.isfile(tm_eq0):
                seed_eqdsk = self._init_files[0]
                if self._test_eqdsk_tx_config(seed_eqdsk):
                    relax_eq = seed_eqdsk
                    self._log(
                        f'Loop {self._current_loop}: inter-loop relax: TM EQDSK (t_idx={i_eq}) from loop {prev_lp} '
                        f'unavailable ({os.path.basename(tm_eq0)}), using seed EQDSK.'
                    )
                else:
                    raise ValueError(
                        f'Loop {self._current_loop}: inter-loop relax needs {tm_eq0} but it is missing '
                        f'or skipped and seed EQDSK is not valid for TORAX: {seed_eqdsk}'
                    )
            self._print(
                f'  TORAX: Running relax ({self._relax_duration:g} s) simulation...'
            )
            prescribed_profiles = None
            if self._steady_state_mode and self._steady_state_tx_seed:
                prescribed_profiles = copy.deepcopy(self._steady_state_tx_seed)
                self._log(
                    f'Loop {self._current_loop}: steady_state_mode inter-loop relax seeded from '
                    f'previous loop t_final profiles on EQDSK {os.path.basename(relax_eq)}.'
                )
            # If prescribed_profiles is None: existing behavior (user kinetics, psi from geometry).
            # In steady_state_mode with a captured seed: relax starts from previous loop t_final
            # profiles on the final-previous-loop EQDSK, then relaxed outputs seed main TORAX.
            self._run_tx_relax(stage='interloop', eqdsk_path=relax_eq, prescribed_profiles=prescribed_profiles)

        with self._loop0_coarse_tx_main_scope():
            myconfig = self._get_tx_config()
            self._print('  TORAX: running simulation...')
            try:
                data_tree, hist, state_list = self._run_torax(myconfig)
            except Exception as e:
                self._print(f'  TORAX: config/init FAILED — {e}')
                raise

            self._data_tree = data_tree  # store for visualization at full TORAX resolution
            # (set even on sim_error so partial data up to the failure is plottable)

            if hist.sim_error != torax.SimError.NO_ERROR:
                self._print(f'  TORAX: sim FAILED ({hist.sim_error})')
                raise ValueError(f'TORAX failed to run the simulation: {hist.sim_error}')

            # Capture the L/H confinement state (active in ADAPTIVE_SOURCE + formation model).
            self._capture_confinement_mode(state_list)

            try:
                self._capture_relax_tx_profiles_from_datatree(data_tree, self._t_init)
                if self._relax_profiles_snapshot is not None:
                    self._relax_mainrun_profile_history.append(
                        {
                            'loop': int(self._current_loop),
                            'profiles': copy.deepcopy(self._relax_profiles_snapshot),
                        }
                    )
            except Exception as _e:
                self._log(f'Warning: could not snapshot TORAX profiles at t_init for next relax: {_e}')

            v_loops = np.zeros(len(self._tm_times))
            for i, t in enumerate(self._tm_times):
                self._tx_update(i, data_tree)
                v_loops[i] = data_tree.scalars.v_loop_lcfs.sel(time=t, method='nearest')

            # internal_manual pedestal: after an IBC-on loop, match the imposed inner-edge
            # value+gradient to the TORAX-evolved core so the next loop's IBC is smooth.
            if self._manual_ibc_active():
                self._smooth_ibc_inner_edge()

            if self._save_outputs:
                self._res_update(data_tree)

            # Flux consumption: positive Ip drives psi_lcfs down in TM-native.
            # Convention: consumed_flux > 0 means plasma has consumed flux from CS.
            consumed_flux = -2.0 * np.pi * (self._state['psi_lcfs_tx'][-1] - self._state['psi_lcfs_tx'][0])
            consumed_flux_integral = np.trapezoid(v_loops[0:], self._tm_times[0:])
            self._log(f"Loop {self._current_loop} TORAX: cflux={consumed_flux:.4f} Wb")
            self._print(f'  TORAX: done (cflux={consumed_flux:.4f} Wb)')
            if self._steady_state_mode:
                try:
                    self._capture_steady_state_tx_seed(data_tree)
                except Exception as _e:
                    self._log(f'steady_state_mode: could not capture t_final profiles for next loop: {_e}')
            return consumed_flux, consumed_flux_integral

    def _tx_update(self, i, data_tree, tx_time=None):
        r'''! Update the simulation state from TORAX results at timestep i.
        
                If sawtooth averaging is enabled, all profile and scalar extractions 
                use time-averaged methods to smooth sawtooth oscillations.
        
                @param i Timestep index (TokaMaker time index; TM comparison fields in plots use this key).
                @param data_tree Result object from Torax.
                @param tx_time If not None, use this time (s) for all TORAX extractions instead of
                       self._tm_times[i] (e.g. last TORAX output time after a failed run).
                
        '''
        t = float(tx_time) if tx_time is not None else self._tm_times[i]

        # ── Scalars ─────────────────────────────────────────────────────────
        self._state['Ip'][i]        = self._extract_tx_scalar(data_tree, 'Ip', t)
        self._state['Ip_tx'][i]     = self._state['Ip'][i]
        self._state['Ip_ni_tx'][i]  = self._extract_tx_scalar(data_tree, 'I_non_inductive', t)
        self._state['pax'][i]       = self._extract_tx_scalar_at_rho(data_tree, 'pressure_thermal_total', t, 0.0)
        self._state['beta_pol'][i]  = self._extract_tx_scalar(data_tree, 'beta_pol', t)
        self._state['beta_N_tx'][i] = self._extract_tx_scalar(data_tree, 'beta_N', t)
        self._state['vloop_tx'][i]  = self._extract_tx_scalar(data_tree, 'v_loop_lcfs', t)
        self._state['f_GW'][i]      = self._extract_tx_scalar(data_tree, 'fgw_n_e_line_avg', t)
        self._state['f_GW_vol'][i]  = self._extract_tx_scalar(data_tree, 'fgw_n_e_volume_avg', t)
        self._state['q95'][i]       = self._extract_tx_scalar(data_tree, 'q95', t)
        self._state['q0'][i]        = self._extract_tx_scalar_at_rho(data_tree, 'q', t, 0.0, rho_coord='rho_face_norm')

        # ── Source profiles for GS solve (FF', p', resistivity) ─────────────
        self._state['ffp_prof'][i]  = self._extract_tx_profile(data_tree, 'FFprime', t)
        self._state['pp_prof'][i]   = self._extract_tx_profile(data_tree, 'pprime',  t)
        self._state['p_prof_tx'][i] = self._extract_tx_profile(data_tree, 'pressure_thermal_total', t)

        self._state['ffp_prof_tx'][i] = self._extract_tx_profile(data_tree, 'FFprime', t)
        self._state['ffp_prof_tx'][i]['y'] *= -2.0 * np.pi  # TX → TM units

        self._state['pp_prof_tx'][i] = self._extract_tx_profile(data_tree, 'pprime', t)
        self._state['pp_prof_tx'][i]['y'] *= -2.0 * np.pi  # TX → TM units

        conductivity = self._extract_tx_profile(data_tree, 'sigma_parallel', t, load_into_state=None)
        self._state['eta_prof'][i] = {'x': self._psi_N.copy(), 'y': 1.0 / conductivity, 'type': 'linterp'}

        # ── Current density profiles ────────────────────────────────────────
        self._state['j_tot'][i]              = self._extract_tx_profile(data_tree, 'j_total',           t, profile_type='jphi-linterp')
        self._state['j_ohmic'][i]            = self._extract_tx_profile(data_tree, 'j_ohmic',           t, profile_type='jphi-linterp')
        self._state['j_ni'][i]               = self._extract_tx_profile(data_tree, 'j_non_inductive',   t, profile_type='jphi-linterp')
        self._state['j_bootstrap'][i]        = self._extract_tx_profile(data_tree, 'j_bootstrap',       t, profile_type='jphi-linterp')
        self._state['j_ecrh'][i]             = self._extract_tx_profile(data_tree, 'j_ecrh',            t, profile_type='jphi-linterp')
        self._state['j_external'][i]         = self._extract_tx_profile(data_tree, 'j_external',        t, profile_type='jphi-linterp')
        self._state['j_generic_current'][i]  = self._extract_tx_profile(data_tree, 'j_generic_current', t, profile_type='jphi-linterp')

        self._state['R_inv_avg_tx'][i] = self._extract_tx_profile(data_tree, 'gm9', t)

        ffp_ni = self._calc_tx_ffp_ni(i)
        self._state['ffp_ni_prof'][i] = {'x': self._psi_N.copy(), 'y': ffp_ni.copy(), 'type': 'linterp'}

        # ── Kinetic profiles ────────────────────────────────────────────────
        self._state['T_i'][i]  = self._extract_tx_profile(data_tree, 'T_i', t)
        self._state['T_e'][i]  = self._extract_tx_profile(data_tree, 'T_e', t)
        self._state['n_i'][i]  = self._extract_tx_profile(data_tree, 'n_i', t)
        self._state['n_e'][i]  = self._extract_tx_profile(data_tree, 'n_e', t)
        self._state['ptot'][i] = self._extract_tx_profile(data_tree, 'pressure_thermal_total', t)

        # ── q profile ──────────────────────────────────────────────────────
        self._state['q_prof_tx'][i] = self._extract_tx_profile(data_tree, 'q', t)

        # ── Psi profile (flux surfaces) ─────────────────────────────────────
        # TORAX stores psi in COCOS 11 (Wb, grows outward, positive-Ip sign).
        # TokaMaker native convention: Wb/rad, psi_axis > psi_lcfs (decreasing outward),
        # with overall sign determined by the equilibrium (can be positive or negative).
        # Transform: psi_TM = -psi_COCOS11 / (2pi)
        psi_tx = -self._extract_tx_profile(data_tree, 'psi', t, load_into_state=None) / (2.0 * np.pi)
        self._state['psi_tx'][i] = {'x': self._psi_N.copy(), 'y': psi_tx.copy(), 'type': 'linterp'}
        self._state['psi_lcfs_tx'][i] = self._state['psi_tx'][i]['y'][-1]
        self._state['psi_axis_tx'][i] = self._state['psi_tx'][i]['y'][0]

        # ── Transport coefficient profiles (chi, D) ─────────────────────────
        for chi_key in ['chi_neo_e', 'chi_neo_i', 'chi_etg_e', 'chi_itg_e', 'chi_itg_i',
                        'chi_tem_e', 'chi_tem_i', 'chi_turb_e', 'chi_turb_i']:
            try:
                self._state[chi_key][i] = self._extract_tx_profile(data_tree, chi_key, t)
            except (KeyError, AttributeError):
                pass

        for d_key in ['D_itg_e', 'D_neo_e', 'D_tem_e', 'D_turb_e']:
            try:
                self._state[d_key][i] = self._extract_tx_profile(data_tree, d_key, t)
            except (KeyError, AttributeError):
                pass

    def _calc_tx_ffp_ni(self, i):
        r'''! Calculate non-inductive FF' profile from TORAX current densities.
                
                The full GS relation is:
                    FF'_total = 2 * mu_0 * (j_tor + p' * <R>) / <1/R>
                
                To avoid double-counting p' when decomposing into inductive/non-inductive:
                    FF'_NI = 2 * mu_0 * j_NI / <1/R>
                    FF'_I  = 2 * mu_0 * (j_I + p' * <R>) / <1/R>
                
                @param i Time index
                @return FF'_NI profile array
                
        '''
        R_inv_avg = self._state['R_inv_avg_tx'][i]['y']

        j_ni = self._state['j_ni'][i]['y']
        ffp_ni = np.where(R_inv_avg != 0, mu_0 * j_ni / R_inv_avg, 0.0)

        return ffp_ni

    def _res_update(self, data_tree):
        r'''! Populate self._results from TORAX data with time-averaging.'''

        self._results['t_res'] = self._tm_times

        # ── Per-timepoint profiles (copy from _state to avoid redundant extraction) ──
        for i, t in enumerate(self._tm_times):
            self._results['T_e'][t] = self._state['T_e'][i]
            self._results['T_i'][t] = self._state['T_i'][i]
            self._results['n_e'][t] = self._state['n_e'][i]
            self._results['q'][t]   = self._state['q_prof_tx'][i]

        # ── Scalar time-series (with time-averaging) ────────────────────────
        self._results['E_fusion']      = self._extract_tx_scalar_timeseries(data_tree, 'E_fusion')
        self._results['Q']             = self._extract_tx_scalar_timeseries(data_tree, 'Q_fusion')
        self._results['Ip']            = self._extract_tx_scalar_timeseries(data_tree, 'Ip')
        self._results['B0']            = self._extract_tx_scalar_timeseries(data_tree, 'B_0')
        self._results['n_e_line_avg']  = self._extract_tx_scalar_timeseries(data_tree, 'n_e_line_avg')
        self._results['n_i_line_avg']  = self._extract_tx_scalar_timeseries(data_tree, 'n_i_line_avg')
        self._results['beta_N']        = self._extract_tx_scalar_timeseries(data_tree, 'beta_N')
        self._results['q95']           = self._extract_tx_scalar_timeseries(data_tree, 'q95')
        self._results['H98']           = self._extract_tx_scalar_timeseries(data_tree, 'H98')
        self._results['v_loop_lcfs']   = self._extract_tx_scalar_timeseries(data_tree, 'v_loop_lcfs')
        self._results['li3']           = self._extract_tx_scalar_timeseries(data_tree, 'li3')
        self._results['P_alpha_total'] = self._extract_tx_scalar_timeseries(data_tree, 'P_alpha_total')
        self._results['P_aux_total']   = self._extract_tx_scalar_timeseries(data_tree, 'P_aux_total')
        self._results['P_ohmic_e']     = self._extract_tx_scalar_timeseries(data_tree, 'P_ohmic_e')
        self._results['P_radiation_e'] = self._extract_tx_scalar_timeseries(data_tree, 'P_radiation_e', scale=-1.0)
        self._results['P_SOL_total']   = self._extract_tx_scalar_timeseries(data_tree, 'P_SOL_total')
        self._results['f_ni']          = self._extract_tx_scalar_timeseries(data_tree, 'f_non_inductive')
        self._results['I_ni']          = self._extract_tx_scalar_timeseries(data_tree, 'I_non_inductive')
        self._results['beta_pol']      = self._extract_tx_scalar_timeseries(data_tree, 'beta_pol')

        # ── Core values (profile at rho=0) ──────────────────────────────────
        self._results['n_e_core'] = self._extract_tx_scalar_at_rho_timeseries(data_tree, 'n_e', 0.0)
        self._results['n_i_core'] = self._extract_tx_scalar_at_rho_timeseries(data_tree, 'n_i', 0.0)
        self._results['T_e_core'] = self._extract_tx_scalar_at_rho_timeseries(data_tree, 'T_e', 0.0)
        self._results['T_i_core'] = self._extract_tx_scalar_at_rho_timeseries(data_tree, 'T_i', 0.0)

        # ── Optional scalars (may not exist in all TORAX versions) ──────────
        try:
            self._results['W_thermal'] = self._extract_tx_scalar_timeseries(data_tree, 'W_thermal')
        except AttributeError:
            pass
        try:
            self._results['tau_E'] = self._extract_tx_scalar_timeseries(data_tree, 'tau_E')
        except AttributeError:
            pass
        try:
            self._results['f_bootstrap'] = self._extract_tx_scalar_timeseries(data_tree, 'f_bootstrap')
        except AttributeError:
            pass

        # ── Greenwald fraction (already time-averaged in _tx_update) ────────
        self._results['f_GW'] = {
            'x': list(self._tm_times),
            'y': np.array(self._state['f_GW']),
        }

        # ── Peak values for quick access ────────────────────────────────────
        Q_arr = self._results.get('Q', {}).get('y', np.array([0]))
        self._results['Q_max'] = float(np.nanmax(Q_arr)) if len(Q_arr) > 0 else 0.0
        self._results['Q_avg_flattop'] = float(np.nanmean(Q_arr[self._flattop])) if np.any(self._flattop) and len(Q_arr) == len(self._flattop) else 0.0

        # ── TokaMaker state arrays (for visualization) ──────────────────────
        self._results['Ip_tm']       = {'x': list(self._tm_times), 'y': np.array(self._state['Ip_tm'])}
        self._results['Ip_tx']       = {'x': list(self._tm_times), 'y': np.array(self._state['Ip_tx'])}
        self._results['Ip_ni_tx']    = {'x': list(self._tm_times), 'y': np.array(self._state['Ip_ni_tx'])}
        self._results['psi_lcfs_tm'] = {'x': list(self._tm_times), 'y': np.array(self._state['psi_lcfs_tm'])}
        self._results['psi_axis_tm'] = {'x': list(self._tm_times), 'y': np.array(self._state['psi_axis_tm'])}
        self._results['psi_lcfs_tx'] = {'x': list(self._tm_times), 'y': np.array(self._state['psi_lcfs_tx'])}
        self._results['psi_axis_tx'] = {'x': list(self._tm_times), 'y': np.array(self._state['psi_axis_tx'])}
        self._results['vloop_tm']    = {'x': list(self._tm_times), 'y': np.array(self._state['vloop_tm'])}
        self._results['vloop_tx']    = {'x': list(self._tm_times), 'y': np.array(self._state['vloop_tx'])}
        self._results['beta_N_tm']   = {'x': list(self._tm_times), 'y': np.array(self._state['beta_N_tm'])}
        self._results['l_i_tm']      = {'x': list(self._tm_times), 'y': np.array(self._state['l_i_tm'])}
        self._results['q95_tm']      = {'x': list(self._tm_times), 'y': np.array(self._state['q95_tm'])}
        self._results['q0_tm']       = {'x': list(self._tm_times), 'y': np.array(self._state['q0_tm'])}
        self._results['pax']         = {'x': list(self._tm_times), 'y': np.array(self._state['pax'])}
        self._results['pax_tm']      = {'x': list(self._tm_times), 'y': np.array(self._state['pax_tm'])}


    # ─── TokaMaker (TM) Methods ─────────────────────────────────────────────────

    def _run_tm(self):
        r'''! Run the GS solve across n timesteps using TokaMaker.
                @return Tuple (consumed_flux, consumed_flux_integral).
                
        '''
        from tqdm import tqdm # creates progress bars
        self._log(f"Loop {self._current_loop} TokaMaker:")

        self._eqdsk_skip = []
        _loop_level_log = []

        n_tm = len(self._tm_times)
        _stride = int(getattr(self, '_loop0_tm_stride', 2))
        _eff_stride = _stride if self._current_loop == 0 else 1
        solve_idx_list = self._loop0_tm_solve_indices(n_tm, _eff_stride)
        solve_idx_set = set(solve_idx_list)
        prev_solve_for = {}
        _last_si = None
        for _si in solve_idx_list:
            prev_solve_for[_si] = _last_si
            _last_si = _si

        if self._current_loop == 0 and len(solve_idx_list) < n_tm:
            self._log(f'Loop 0 TokaMaker: subsampling {len(solve_idx_list)}/{n_tm} timesteps (stride={_stride}).')
            self._print(f'  TokaMaker: {len(solve_idx_list)}/{n_tm} GS solves (loop 0 subsampled, stride={_stride})...')
        else:
            self._print(f'  TokaMaker: solving {len(self._tm_times)} equilibria...')

        # ── Per-loop initialization (before timestep sweep) ──────────────────
        self._state['psi_grid_prev_tm'] = {}

        if (self._steady_state_mode
                and not self._coupling_iteration_is_first()
                and self._steady_state_tm_psi_seed is not None):
            self._psi_warm_start[0] = np.asarray(self._steady_state_tm_psi_seed, dtype=float).copy()

        # Seed coil regularization targets
        # First coupling pass: use zero targets (None) so the solver freely finds the correct
        # coil configuration without being biased toward the initial equilibrium.
        # Later coupling loops: seed from the last solve of the previous loop for warm-starting.
        if getattr(self, '_coil_reg_config', None):
            init_targets = None
            if self._current_loop > 0 and not self._coupling_iteration_is_first():
                try:
                    init_targets, _ = self._tm.get_coil_currents()
                except Exception as e:
                    self._log(f'TM: could not read initial equilibrium coil currents: {e}')
            self._apply_tm_coil_reg(targets=init_targets)

        # Debug: log coil bounds at the start of each loop
        if self._debug_mode and hasattr(self, '_coil_bounds') and self._coil_bounds:
            self._log('  TM coil bounds [A/turn] (turns * A/turn = total A-turns):')
            for cname, (lo, hi) in self._coil_bounds.items():
                n_turns = self._tm.coil_sets.get(cname, {}).get('net_turns', 1.0)
                self._log(f'    {cname}: [{lo:.3g}, {hi:.3g}] A/turn  x {n_turns:.0f} turns = [{lo*n_turns:.3g}, {hi*n_turns:.3g}] A-turns')


        if 0 in self._psi_warm_start and self._psi_warm_start[0] is not None:
            self._tm.set_psi_dt(psi0=self._psi_warm_start[0], dt=1.0e10)

        # Progress bar counts actual GS attempts (subsampled loop 0 → 13/13 not 13/25).
        _pbar_total = len(solve_idx_list)
        with tqdm(total=_pbar_total,
                  desc=f'  TM loop {self._current_loop}', unit='solve',
                  bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_inv_fmt}]{postfix}'
                  ) as _pbar:
            for i, t in enumerate(self._tm_times):
                if i not in solve_idx_set:
                    _loop_level_log.append({
                        'i': i, 't': t, 'skipped': True, 'succeeded': False,
                        'level': None, 'level_name': None, 'error': None,
                    })
                    continue

                prev_tm_idx = prev_solve_for[i]

                # Clear isoflux, flux, and saddle targets from previous timepoint
                self._tm.set_isoflux_constraints(None)
                self._tm.set_psi_constraints(None, None)
                self._tm.set_saddle_constraints(None)

                Ip_target = abs(self._state['Ip'][i])
                P0_target = abs(self._state['pax'][i])
        
                self._tm.set_targets(Ip=Ip_target, pax=P0_target) # using pax target with j_phi inputs 
                self._tm.set_resistivity(eta_prof=self._state['eta_prof'][i])
        
                ffp_prof = {'x': self._state['ffp_prof'][i]['x'].copy(),
                               'y': self._state['ffp_prof'][i]['y'].copy(),
                               'type': self._state['ffp_prof'][i]['type']}
                pp_prof = {'x': self._state['pp_prof'][i]['x'].copy(),
                              'y': self._state['pp_prof'][i]['y'].copy(),
                              'type': self._state['pp_prof'][i]['type']}
        
                lcfs = self._state['lcfs_geo'][i]

                # Set saddle-point (X-point) constraints during diverted phase
                use_x_points = (
                    self._x_point_targets is not None
                    and self._diverted_times is not None
                    and self._diverted_times[0] <= t <= self._diverted_times[1]
                )
                if use_x_points:
                    saddle_weights = self._x_point_weight * np.ones(self._x_point_targets.shape[0])
                    self._tm.set_saddle_constraints(self._x_point_targets, saddle_weights)

                    # trims lcfs targets near X-point(s)
                    if self._trim_lcfs:
                        perc_limit = self._trim_lcfs_perc_limit       # LCFS points above percentage limit* max(abs(Z)) are removed from isoflux targets
                        Z_max_abs = np.max(np.abs(lcfs[:, 1]))
                        Z_lim = perc_limit * Z_max_abs
                        if np.shape(self._x_point_targets)[0] == 1 and self._x_point_targets[0][1] > 0: # upper single null
                            lcfs = lcfs[lcfs[:, 1] <= Z_lim]
                        elif np.shape(self._x_point_targets)[0] == 1 and self._x_point_targets[0][1] < 0: # lower single null
                            lcfs = lcfs[lcfs[:, 1] >= -Z_lim]
                        elif np.shape(self._x_point_targets)[0] == 2: # double null
                            lcfs = lcfs[np.abs(lcfs[:, 1]) <= Z_lim]

                # When diverted, add manually specified strike points to isoflux targets
                if use_x_points and self._strike_point_targets is not None:
                    self._state['strike_pts'][i] = self._strike_point_targets
                    lcfs = np.vstack([lcfs, self._strike_point_targets])
                else:
                    self._state['strike_pts'][i] = np.empty((0, 2))

                isoflux_weights = LCFS_WEIGHT * np.ones(len(lcfs))
                lcfs_psi_target = self._state['psi_lcfs_tx'][i] # _state in Wb/rad, TM uses Wb/rad (AKA Wb-rad)

                # set_isoflux on all LCFS points for lcfs shape targets
                self._tm.set_isoflux_constraints(lcfs, isoflux_weights) # shape targets

                # Pick outboard midplane point (largest R at approx Z = Z_axis)
                z_axis = self._state['Z'][i]
                omp_idx = np.argmax(lcfs[:, 0] * np.exp(-0.5 * ((lcfs[:, 1] - z_axis) / (0.3 * self._state['a'][i]))**2))
                omp_point = lcfs[omp_idx:omp_idx+1, :]  # shape (1, 2)
                # Set lcfs psi value target (from TORAX) only at midplane outboard side of lcfs.
                self._tm.set_psi_constraints(omp_point, targets=np.array([lcfs_psi_target]),
                                             weights=np.array([LCFS_WEIGHT * 10.])) # psi value target
        
        
                self._tm.update_settings()

                if prev_tm_idx is not None:
                    _psi_prev = self._state['psi_grid_prev_tm'].get(prev_tm_idx)
                    if _psi_prev is not None:
                        self._tm.set_psi_dt(
                            psi0=_psi_prev,
                            dt=self._tm_times[i] - self._tm_times[prev_tm_idx],
                        )

                skip_coil_update = False
                eq_name = os.path.join(self._eqdsk_dir, f'{self._current_loop:03d}.{i:03d}.eqdsk')
                step_fail_msg = None  # set when TM GS succeeds but TORAX rejects all EQDSK resolutions

                solve_succeeded = False
                level_attempts = []

                ffp_prof_raw = copy.deepcopy(ffp_prof)
                pp_prof_raw  = copy.deepcopy(pp_prof)
        
                # Pre-calculate all level profiles
                level_profiles = []

                # Level 0: jphi
                # ffp_0 = self._state['j_tot'][i]
                # pp_0 = pp_prof
                # level_profiles.append({'ffp': ffp_0, 'pp': pp_0, 'name': 'lv0: jphi'})

                # Level 1: raw
                ffp_1, pp_1 = self._tm_prof_input_raw(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
                level_profiles.append({'ffp': ffp_1, 'pp': pp_1, 'name': 'raw tx profs'})
        
                # Level 2: sign flip
                ffp_2, pp_2 = self._tm_prof_input_sign_flip(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
                level_profiles.append({'ffp': ffp_2, 'pp': pp_2, 'name': 'sign_flip'})
        
                # Level 3: pedestal smoothing (takes p_profile as input) # TODO: read in actual n_rho_ped_top, have to add to state first
                ffp_3, pp_3 = self._tm_prof_input_ped_smooth(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw), copy.deepcopy(self._state['p_prof_tx'][i])) 
                level_profiles.append({'ffp': ffp_3, 'pp': pp_3, 'name': 'ped_smoothing'})
        
                # Level 4: power flux
                ffp_4, pp_4 = self._tm_prof_input_power_flux(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
                level_profiles.append({'ffp': ffp_4, 'pp': pp_4, 'name': 'analytic'})

                # Try each level
                for level_idx, level_prof in enumerate(level_profiles):
                    level_name = level_prof['name']
                    ffp_level = level_prof['ffp']
                    pp_level = level_prof['pp']

                    # Initialize psi from geometry parameters
                    self._tm.init_psi(self._state['R0_mag'][i], self._state['Z'][i], self._state['a'][i], self._state['kappa'][i], self._state['delta'][i])

                    if i in self._psi_warm_start and self._psi_warm_start[i] is not None:
                        self._tm.set_psi(self._psi_warm_start[i], update_bounds=True)
                    elif prev_tm_idx is not None and self._state['psi_grid_prev_tm'].get(prev_tm_idx) is not None:
                        self._tm.set_psi(self._state['psi_grid_prev_tm'][prev_tm_idx], update_bounds=True)

                    try:
                        self._tm.set_profiles(ffp_prof=ffp_level, pp_prof=pp_level, ffp_NI_prof=self._state['ffp_ni_prof'][i])

                        if self._output_mode == 'debug': # allows TM terminal outputs in debug mode
                            self._state['equil'][i] = self._tm.solve()
                        else:
                            with self._quiet_tm(): # silences TM terminal outputs
                                self._state['equil'][i] = self._tm.solve()

                        level_attempts.append({'level': level_idx, 'name': level_name,
                                              'ffp': ffp_level, 'pp': pp_level,
                                              'succeeded': True, 'error': None})
                        ffp_prof, pp_prof = ffp_level, pp_level
                        solve_succeeded = True
                        break
                    except Exception as e:
                        level_attempts.append({'level': level_idx, 'name': level_name,
                                              'ffp': ffp_level, 'pp': pp_level,
                                              'succeeded': False, 'error': str(e)})
                        if self._output_mode == 'debug' and self._out_dir is not None:
                            _ldiag_name = f'tm_diag_loop{self._current_loop:03d}_tidx{i:03d}_lv{level_idx:02d}.png'
                            if self._output_file_tag is not None:
                                _ldiag_name = f'{self._output_file_tag}_{_ldiag_name}'
                            try:
                                tm_diagnostic_plot(
                                    self, i, t, level_attempts, solve_succeeded=False,
                                    save_path=os.path.join(self._out_dir, _ldiag_name),
                                    display=False, tm_gs_ok=False, step_error=str(e),
                                )
                            except Exception as _le:
                                self._log(f'tm_diagnostic_plot (level {level_idx}) failed at i={i}: {_le}')

                if not solve_succeeded:
                    self._eqdsk_skip.append(eq_name)
                    skip_coil_update = True
                    self._log(f'\tTM: Solve failed at t={t} (all levels attempted).')
                    self._state['psi_grid_prev_tm'][i] = None  # if solve failed, set psi grid to None
                    # In debug mode per-level plots are already saved inside the except block above.
                    if self._out_dir is not None and self._output_mode != 'debug':
                        _last_err = level_attempts[-1].get('error') if level_attempts else None
                        diag_name = f'tm_diag_loop{self._current_loop:03d}_tidx{i:03d}.png'
                        if self._output_file_tag is not None:
                            diag_name = f'{self._output_file_tag}_{diag_name}'
                        _diag_path = os.path.join(self._out_dir, diag_name)
                        tm_diagnostic_plot(
                            self, i, t, level_attempts, solve_succeeded,
                            save_path=_diag_path, display=False,
                            tm_gs_ok=False, step_error=_last_err,
                        )

                tm_gs_ok = solve_succeeded
                if solve_succeeded:
                    torax_accepted = False
                    _n_attempts = len(EQDSK_SAVE_NR_NZ_SEQUENCE)
                    for _attempt_idx, nr_nz in enumerate(EQDSK_SAVE_NR_NZ_SEQUENCE):
                        quiet_test = (_attempt_idx < _n_attempts - 1)
                        with self._quiet_tm():
                            self._state['equil'][i].save_eqdsk(eq_name,
                                lcfs_pad=1-self._last_surface_factor, run_info='TokaMaker EQDSK',
                                cocos=self._cocos, nr=nr_nz, nz=nr_nz,
                                truncate_eq=self._truncate_eq)
                        if self._test_eqdsk_tx_config(eq_name, quiet=quiet_test):
                            torax_accepted = True
                            if _attempt_idx > 0 and self._output_mode == 'debug':
                                base_nz = EQDSK_SAVE_NR_NZ_SEQUENCE[0]
                                # self._print(
                                #     f'    EQDSK {os.path.basename(eq_name)}: base nr=nz={base_nz} rejected by TORAX; '
                                #     f'accepted at nr=nz={nr_nz}.'
                                # )
                            break
                    if not torax_accepted:
                        # Same coupling outcome as a failed TM solve: skip this time for TORAX geometry
                        # (TORAX interpolates neighbors in _get_tx_config) instead of aborting the whole fly.
                        self._eqdsk_skip.append(eq_name)
                        skip_coil_update = True
                        self._state['psi_grid_prev_tm'][i] = None
                        step_fail_msg = (
                            f'TORAX rejected EQDSK after save attempts nr=nz in '
                            f'{EQDSK_SAVE_NR_NZ_SEQUENCE}: {os.path.basename(eq_name)}'
                        )
                        self._log(f'\tTM: {step_fail_msg}')
                        solve_succeeded = False
                    else:
                        self._tm_update(i)

                        # Store diverted/limited flag for this timestep
                        if not hasattr(self, '_diverted_flags'):
                            self._diverted_flags = {}
                        self._diverted_flags[i] = self._state['equil'][i].diverted
                        # Log it: a timestep that should be diverted (inside the
                        # diverted window) but solves LIMITED means the X-point was
                        # not held -> shape targets missed for that frame.
                        _in_div_window = (
                            self._diverted_times is not None
                            and self._diverted_times[0] <= t <= self._diverted_times[1]
                        )
                        if _in_div_window and not self._state['equil'][i].diverted:
                            self._log(f'  [shape] t={t:.2f}s: in diverted window but '
                                      f'solved LIMITED (X-point not held).')

                        # Store psi on nodes for later movie generation
                        self._tm_psi_on_nodes.setdefault(self._current_loop, {})[i] = self._state['equil'][i].get_psi(normalized=False)

                _winning = next((a for a in level_attempts if a['succeeded']), None)
                _last_attempt = level_attempts[-1] if level_attempts else {}
                _log_err = step_fail_msg if step_fail_msg else (
                    _last_attempt.get('error') if not solve_succeeded else None
                )
                _loop_level_log.append({
                    'i': i, 't': t,
                    'succeeded': solve_succeeded,
                    'level': _winning['level'] if _winning else None,
                    'level_name': _winning['name'] if _winning else None,
                    'error': _log_err,
                })

                if self._output_mode in ('debug', 'normal') and self._out_dir is not None:
                    # In debug mode, per-level failure plots are saved inside the solve loop above;
                    # only emit the combined summary plot here on success.
                    if self._output_mode == 'debug' and solve_succeeded:
                        diag_name = f'tm_diag_loop{self._current_loop:03d}_tidx{i:03d}.png'
                        if self._output_file_tag is not None:
                            diag_name = f'{self._output_file_tag}_{diag_name}'
                        _diag_path = os.path.join(self._out_dir, diag_name)
                        try:
                            tm_diagnostic_plot(
                                self, i, t, level_attempts, solve_succeeded,
                                save_path=_diag_path, display=False,
                                tm_gs_ok=tm_gs_ok, step_error=_log_err,
                            )
                        except Exception as _e:
                            self._log(f'tm_diagnostic_plot failed at i={i}: {_e}')
                    if solve_succeeded:
                        prof_name = f'profile_loop{self._current_loop:03d}_tidx{i:03d}.png'
                        if self._output_file_tag is not None:
                            prof_name = f'{self._output_file_tag}_{prof_name}'
                        _prof_path = os.path.join(self._out_dir, prof_name)
                        try:
                            profile_plot(self, i, t, save_path=_prof_path, display=False)
                        except Exception as _e:
                            self._log(f'profile_plot failed at i={i}: {_e}')

                    # Manual-pedestal IBC diagnostic (debug mode, when the IBC is active).
                    if (self._output_mode == 'debug' and solve_succeeded
                            and self._ped_mode == 'internal_manual' and self._manual_ibc_active()):
                        ped_name = f'IBC_debug_loop{self._current_loop:03d}_tidx{i:03d}.png'
                        if self._output_file_tag is not None:
                            ped_name = f'{self._output_file_tag}_{ped_name}'
                        try:
                            IBC_debug(self, i, t,
                                      save_path=os.path.join(self._out_dir, ped_name))
                        except Exception as _e:
                            self._log(f'IBC_debug failed at i={i}: {_e}')

                # Update progress bar postfix; print FAIL messages above the bar
                if solve_succeeded:
                    lvl = _winning['level']
                    _pbar.set_postfix_str(f't={t:.2f}s OK(L{lvl})', refresh=False)
                else:
                    err_short = (step_fail_msg or _last_attempt.get('error') or 'unknown')[:60]
                    _fail_tag = 'TORAX EQDSK' if step_fail_msg else 'TM'
                    tqdm.write(f'    WARNING: {_fail_tag} FAIL at t={t:.2f}s — {err_short}')
                    self._log(f'    {_fail_tag} FAIL at t={t:.2f}s — {err_short}')
                    _pbar.set_postfix_str(f't={t:.2f}s FAIL', refresh=False)

                if not skip_coil_update and getattr(self, '_coil_reg_config', None):
                    prev_coil_targets, _ = self._state['equil'][i].get_coil_currents()
                    self._apply_tm_coil_reg(targets=prev_coil_targets)
                _pbar.update(1)

        # Flux consumption: positive Ip drives psi_lcfs down in TM convention
        # Convention: consumed_flux > 0
        consumed_flux = -(self._state['psi_lcfs_tm'][-1] - self._state['psi_lcfs_tm'][0]) * 2.0 * np.pi
        consumed_flux_integral = np.trapezoid(self._state['vloop_tm'][0:], self._tm_times[0:])

        n_ok = sum(1 for e in _loop_level_log if e.get('succeeded'))
        n_skip = sum(1 for e in _loop_level_log if e.get('skipped'))
        n_gs = len(solve_idx_list)
        self._print(f'  TokaMaker: {n_ok}/{n_gs} solved (cflux={consumed_flux:.4f} Wb)')

        # Compact level-usage summary for log
        from collections import Counter
        _lvl_counts = Counter(e['level_name'] for e in _loop_level_log if e.get('succeeded'))
        _lvl_summary = ', '.join(f'{name}: {cnt}' for name, cnt in sorted(_lvl_counts.items()))
        n_fail = sum(1 for e in _loop_level_log if not e.get('skipped') and not e.get('succeeded'))
        # Omit "Skipped: N" when timesteps were subsampled (already stated in the GS-solve preface).
        _skip_note = (
            f' Skipped: {n_skip}.'
            if n_skip and n_gs >= len(self._tm_times) else ''
        )
        self._print(f'\tTM summary: {n_ok}/{n_gs} solved. Levels: {_lvl_summary}.'
                  + _skip_note
                  + (f' Failures: {n_fail}.' if n_fail else ''))

        if self._debug_mode:
            summary_name = f'tm_summary_loop{self._current_loop:03d}.png'
            if self._output_file_tag is not None:
                summary_name = f'{self._output_file_tag}_{summary_name}'
            _summary_path = os.path.join(self._out_dir, summary_name)
            try:
                tm_loop_summary_plot(self, _loop_level_log, save_path=_summary_path, display=False)
            except Exception as _e:
                self._log(f'tm_loop_summary_plot failed: {_e}')

        if self._steady_state_mode:
            i_hi = len(self._tm_times) - 1
            pg = self._state['psi_grid_prev_tm'].get(i_hi)
            if pg is not None:
                self._steady_state_tm_psi_seed = np.asarray(pg, dtype=float).copy()

        return consumed_flux, consumed_flux_integral

    # ── Profile level functions ──────────────────────────────────────────
    # Each level takes (self, ffp_prof, pp_prof, i) and returns (ffp_prof, pp_prof).
    # All levels receive deep copies of the raw TORAX profiles (not cumulative).

    def _tm_prof_input_raw(self, ffp_prof, pp_prof):
        r'''! Raw TORAX profiles passed through unchanged.'''
        return ffp_prof, pp_prof

    def _tm_prof_input_sign_flip(self, ffp_prof, pp_prof):
        r'''! Sign-flip clipping: clip each profile to its dominant sign.'''
        def _clip(prof):
            y = prof['y']
            sign = 1 if np.sum(y > 0) >= np.sum(y < 0) else -1
            y_new = np.clip(y, 0, None) if sign > 0 else np.clip(y, None, 0)
            return {**prof, 'y': y_new}
        return _clip(ffp_prof), _clip(pp_prof)

    def _tm_prof_input_ped_smooth(self, ffp_prof, pp_prof, p_prof, transition_psi_N = 0.6, gauss_sigma=8, blend_width=0.02, sav_window=41, sav_order=3):
        r'''! Edge smoothing with Gaussian filter: smooth p profile and take derivative for pp_prof.'''
        
        # Extract pressure 'y' values and ensure they're 1D
        p = np.atleast_1d(p_prof['y'])
        
        # Handle case where input is empty or scalar
        if p.size == 0:
            return ffp_prof, pp_prof
        
        # First smooth entire profile
        p_smooth = gaussian_filter1d(p, gauss_sigma, mode='nearest')

        # Sigmoid blend weight: 0 = pure original, 1 = pure smoothed
        # Centered at edge_psi, width controlled by blend_width
        blend = 0.5 * (1 + np.tanh((self._psi_N - transition_psi_N) / blend_width))

        # blend original and smoothed profiles so the value and slope are continuous across transition
        p_new = (1 - blend) * p + blend * p_smooth

        pp_new = np.gradient(p_new, self._psi_N)
        pp_new_smooth = savgol_filter(pp_new, sav_window, sav_order)

        # Return modified pp_prof with smoothed values, ffp_prof unchanged
        return ffp_prof, {**pp_prof, 'y': pp_new_smooth}

    def _tm_prof_input_power_flux(self, ffp_prof, pp_prof):
        r'''! Generic power-flux shape.'''
        ffp_out = create_power_flux_fun(N_PSI, 1.5, 2.0)
        pp_out  = create_power_flux_fun(N_PSI, 4.0, 1.0)
        return ffp_out, pp_out

    def _tm_update(self, i):
        r'''! Update internal state and coil current results based on results of GS solver.
                @param i Timestep of the solve.
                
        '''
        eq_stats = self._state['equil'][i].get_stats(li_normalization='iter')
        self._state['Ip'][i] = eq_stats['Ip']
        self._state['Ip_tm'][i] = eq_stats['Ip']
        self._state['pax_tm'][i] = eq_stats['P_ax']
        self._state['beta_N_tm'][i] = eq_stats['beta_n']
        self._state['l_i_tm'][i] = eq_stats['l_i']

        self._state['psi_lcfs_tm'][i] = self._state['equil'][i].psi_bounds[0] # TM outputs in Wb/rad (AKA Wb-rad) which is how psi_lcfs is stored
        self._state['psi_axis_tm'][i] = self._state['equil'][i].psi_bounds[1]
        self._state['psi_tm'][i] = {'x': self._psi_N.copy(), 'y': self._state['psi_axis_tm'][i] + (self._state['psi_lcfs_tm'][i] - self._state['psi_axis_tm'][i]) * self._psi_N, 'type': 'linterp'}

        try:
            self._state['vloop_tm'][i] = self._state['equil'][i].calc_loopvoltage()
        except ValueError:
            self._log(f'WARNING: calc_loopvoltage failed at t-idx = {i} '
                      f'(likely Ip_ni > Ip); setting vloop_tm to 0.')
            self._state['vloop_tm'][i] = 0.0

        # store TokaMaker pressure profile from get_profiles()
        tm_psi, tm_f_prof, tm_fp_prof, tm_p_prof, tm_pp_prof = self._state['equil'][i].get_profiles(npsi=N_PSI)

        self._state['ffp_prof_tm'][i] = {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, tm_psi, tm_fp_prof*tm_f_prof), 'type': 'linterp'}
        self._state['pp_prof_tm'][i] =  {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, tm_psi, tm_pp_prof), 'type': 'linterp'}
        self._state['p_prof_tm'][i] =   {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, tm_psi, tm_p_prof), 'type': 'linterp'}
        self._state['f_prof_tm'][i] =   {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, tm_psi, tm_f_prof), 'type': 'linterp'}

        # pull geo profiles
        psi_geo, q_tm, geo, _, _, _ = self._state['equil'][i].get_q(npsi=N_PSI, psi_pad=1-self._last_surface_factor)

        self._state['q0_tm'][i] = q_tm[0] if len(q_tm) > 0 else np.nan
        self._state['q95_tm'][i] = np.interp(0.95, psi_geo, q_tm) if len(psi_geo) > 0 and len(q_tm) > 0 else np.nan
        self._state['q_prof_tm'][i] = {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, psi_geo, q_tm), 'type': 'linterp'}

        self._state['R_avg_tm'][i] =     {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, psi_geo, np.array(geo[0])), 'type': 'linterp'}
        self._state['R_inv_avg_tm'][i] = {'x': self._psi_N.copy(), 'y': np.interp(self._psi_N, psi_geo, np.array(geo[1])), 'type': 'linterp'}

        # Update Results
        coils, _ = self._state['equil'][i].get_coil_currents()
        if 'COIL' not in self._results:
            self._results['COIL'] = {coil: {} for coil in coils}
        for coil, current in coils.items():
            if coil not in self._results['COIL']:
                self._results['COIL'][coil] = {}
            self._results['COIL'][coil][self._tm_times[i]] = current

        # get psi to use in next timestep
        self._state['psi_grid_prev_tm'][i] = self._state['equil'][i].get_psi(normalized=False)
        self._psi_warm_start[i] = self._state['equil'][i].get_psi(normalized=False)  # persist across steps

        # Extract LCFS contour and X-points from the solved equilibrium
        try:
            lcfs_tm = self._state['equil'][i].trace_surf(1.0)    
        except Exception:
            try:
                lcfs_tm = self._state['equil'][i].trace_lcfs(0.99)  
            except Exception:
                self._state['lcfs_geo_tm'][i] = None

        self._state['lcfs_geo_tm'][i] = np.asarray(lcfs_tm) if lcfs_tm is not None else None    
        
        try:
            x_pts, _ = self._state['equil'][i].get_xpoints()
            self._state['x_pts_tm'][i] = np.asarray(x_pts) if x_pts is not None else None
        except Exception:
            self._state['x_pts_tm'][i] = None


    # ─── I/O & Logging ──────────────────────────────────────────────────────────

    def save_state(self, fname):
        r'''! Save intermediate simulation state to JSON.
                @param fname Filename to save to.
                
        '''
        with open(fname, 'w') as f:
            json.dump(self._state, f, cls=MyEncoder)

    def save_res(self):
        r'''! Save simulation results to JSON.'''
        if self._fname_out is not None:
            with open(self._fname_out, 'w') as f:
                json.dump(self._results, f, cls=MyEncoder)

    def _log(self, msg):
        r'''! Write message to log file only.'''
        if hasattr(self, '_log_file') and self._log_file is not None:
            with open(self._log_file, 'a') as f:
                print(msg, file=f)

    def _print(self, msg):
        r'''! Write message to both stdout and log file.'''
        print(msg)
        self._log(msg)

    def _quiet_tm(self):
        r'''! Context manager: redirect C/Fortran-level stdout+stderr to /dev/null.
                
        '''
        @contextmanager
        def _cm():
            target_fd = os.open(os.devnull, os.O_WRONLY)
            saved_out = os.dup(1)
            saved_err = os.dup(2)
            os.dup2(target_fd, 1)
            os.dup2(target_fd, 2)
            try:
                yield
            finally:
                sys.stdout.flush()
                sys.stderr.flush()
                os.dup2(saved_out, 1)
                os.dup2(saved_err, 2)
                os.close(saved_out)
                os.close(saved_err)
                os.close(target_fd)
        return _cm()

    def configure_redirect_to_log(self):
        r'''! Step 3/3 of setup to divert noisy outputs to log file.
                Captures INFO-level and above.
                
        '''
        if self._logging_configured or not self._log_file:
            return

        root_logger = logging.getLogger()
        file_handler = logging.FileHandler(self._log_file, mode='a')

        root_logger.setLevel(logging.INFO)
        file_handler.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s [%(name)-12s:%(levelname)-8s] %(message)s')
        file_handler.setFormatter(formatter)

        root_logger.addHandler(file_handler)
        self._logging_configured = True
        logging.info(f"File logging configured. All logs will be written to {self._log_file}")

    # =========================================================================
    #  fly — run simulation loop
    # =========================================================================

    def fly(self, run_name='tmp', convergence_threshold=-1.0, max_loop=3,
            output_mode=False, skip_bad_init_eqdsks=False,
            initial_relax=True, relax=False, relax_kinetics=False, relax_duration=1.0, relax_dt=0.1,
            t_ave_toggle='off', t_ave_window=0.5, t_ave_causal=True, t_ave_ignore_start=0.25,
            loop0=False, steady_state_mode=False): # TODO: separate steady_state_mode?
        r'''! Run TokaMaker_TORAX coupled pulse design loop.
        
                @param convergence_threshold Max fractional change in consumed flux between loops for convergence.
                @param max_loop Highest **counted** coupling index to run (inclusive): full-resolution passes
                       use indices 1 … max_loop. The optional cheap pass at index 0 (when loop0=True)
                       is always attempted first and does **not** count toward this limit. Example: max_loop=3
                       runs loop 0 (if enabled) then loops 1, 2, and 3 before stopping unless
                       convergence ends earlier.
                @param run_name Name for this run (used in output directory and log file).
                @param output_mode Output level selector: False (or None), True (alias for 'normal'),
                       'minimal', 'normal', or 'debug'. String values 'false'/'none'/'off'/'no' map to False.
                       When not False, artifacts go under ./TokaMaker_TORAX_outputs/ (tmp/ for
                       run_name='tmp', else {run_name}_{timestamp}/). Unless run_name='tmp', saved
                       filenames are prefixed {run_name}_{timestamp}_. A log file is always written
                       (into the output directory when one exists, otherwise the current working directory).
                       results.json path is set on the instance (self._fname_out) but is not
                       written automatically; call save_res() to persist it. In-memory self._results
                       is updated after each successful TORAX pass for any non-False mode.
                       End-of-run PNG/MP4 saves are skipped inside Jupyter notebooks.
                
                       False / None — No output directory. Log only (TokaMaker_TORAX_log_tmp.log or
                       TokaMaker_TORAX_log_{run_name}_{timestamp}.log in cwd). No plots or config files.
                       TokaMaker gEQDSK files use a temporary directory and are deleted at exit.
                
                       'minimal' — Per completed coupling loop: scalars_loop{N}.png,
                       PLH_components_loop{N}.png. At end of run (non-Jupyter): profile_evolution.png,
                       lcfs_evolution.png, movie_loop{N}.mp4 (N = last completed loop index).
                       On TokaMaker GS failure at a timestep: tm_diag_loop{N}_tidx{i}.png. On TORAX
                       failure: scalars_loop{N}_torax_failed.png and, if partial TORAX data exist,
                       profile_loop{N}_torax_failed_tfinal.png. No TORAX config .py files, no
                       per-timestep profile plots, no relax-profile figures, no persisted gEQDSK files
                       (temporary EQDSK dir removed at exit).
                
                       'normal' (also output_mode=True) — Per loop: scalars_loop{N}.png,
                       PLH_components_loop{N}.png; tx_config_loop{N}.py; initial / inter-loop relax
                       configs tx_config_relax000_initial.py and tx_config_relax_inter_{N}.py; each
                       successful TokaMaker timestep profile_loop{N}_tidx{i}.png. At end (non-Jupyter):
                       profile_evolution.png and movie_loop{N}.mp4 (no lcfs_evolution.png).
                       Same failure plots as minimal. No tm_diag on successful solves, no
                       relax-profile figures, no tm_summary_loop{N}.png, no persisted gEQDSK files.
                
                       'debug' — All normal artifacts plus lcfs_evolution.png at end of run. gEQDSK files
                       {loop:03d}.{i:03d}.eqdsk are kept in the output directory (not deleted).
                       Initial / inter-loop relax: tx_relax_profiles_initial.png,
                       tx_relax_profiles_inter_loop{N}.png. Every TokaMaker timestep (success or fail):
                       tm_diag_loop{N}_tidx{i}.png (successful solves also get profile_loop{N}_tidx{i}.png).
                       After each loop: tm_summary_loop{N}.png. Python logging (TORAX, JAX, etc.) is
                       redirected to the log file; per-loop wall time is printed.
                @param skip_bad_init_eqdsks If True, skip broken initial gEQDSK files instead of raising.
                @param initial_relax If True (default), run a short TORAX relax on the seed EQDSK before
                       the first coupling TM-TORAX pass (flattened user inputs; psi from geometry unless an
                       inter-loop relax already set psi). If False, skip and start from EQDSK psi / loaded config.
                       When relax is True, this is forced True so initial relax establishes psi before
                       later coupling iterations.
                @param relax If True, run an additional short TORAX relax before each coupling iteration
                       with index ≥1, on TM-solved EQDSK (previous_loop).000.eqdsk (fallback: seed), with **user**
                       n_e, T_e, T_i from loaded config and set_*() and psi from geometry on that EQDSK—
                       avoiding drift from previous TORAX outputs between loops.
                       Default False (backward compatible). Implies initial_relax=True.
                @param relax_kinetics If False (default), each relax only evolves current; density and
                       ion/electron heat stay fixed at the profiles present before that relax. If True,
                       relax uses the same evolve_* flags as set_evolve() / loaded config. When True,
                       relaxed n_e, T_e, T_i are injected into main TORAX after set_*() like psi.
                @param relax_duration Duration (s) of each relax simulation. Default 1.0 s.
                @param relax_dt Fixed timestep (s) for each initial / inter-loop TORAX relax run
                       (numerics.fixed_dt). Default 0.1 s.
                @param t_ave_toggle Time-averaging mode: 'off' (no averaging), 'flattop' (average only
                       during flat-top), or 'pulse' (average over the whole pulse).
                @param t_ave_window Averaging window size in seconds. Default 0.5 s.
                @param t_ave_causal If True, window is entirely behind the timepoint (backward-looking).
                       If False, window is centred on the timepoint.
                @param t_ave_ignore_start Ignore the first N seconds of the pulse when building the
                       averaging window (avoids numerical transients). Default 0.25 s.
                @param loop0 If True, run a first coupling pass at index 0 with reduced cost:
                       coarse TORAX radial grid (face_centers linspace, DEFAULT_LOOP0_TX_FACE_POINTS),
                       the same coarse grid on the optional initial relax, and subsampled TokaMaker times
                       (stride 2). If False (default), skip that pass and start coupling at index 1 at full resolution
                       (initial relax is unchanged and controlled only by initial_relax / relax).
                @param steady_state_mode If False (default), each coupling loop after the first reuses the
                       time-dependent TokaMaker EQDSK sequence from the previous loop (existing behavior).
                       If True, after each completed loop the next loop seeds TORAX with psi and kinetics from
                       the previous main run at t_final, uses the final TokaMaker EQDSK from the
                       previous loop for all TORAX geometry times (flat equilibrium shape in time), warm-starts
                       TokaMaker at i=0 from the previous loop's final psi grid, and runs inter-loop relax
                       on that final EQDSK when relax is True.
                
        '''
        import tempfile

        if relax:
            initial_relax = True

        self._steady_state_mode = bool(steady_state_mode)
        self._steady_state_tx_seed = None
        self._steady_state_tm_psi_seed = None

        self._fly_loop0 = bool(loop0)
        self._loop0_coarse_tx = self._fly_loop0
        self._loop0_relax_coarse_tx = self._fly_loop0
        self._loop0_tm_stride = 2 if self._fly_loop0 else 1
        self._loop0_tx_face_points = int(DEFAULT_LOOP0_TX_FACE_POINTS)

        # Disable JAX's persistent XLA compilation cache before any TORAX/JAX JIT
        # compilation occurs, was causing semaphore leaks.
        try:
            import jax
            jax.config.update('jax_enable_compilation_cache', False)
        except Exception:
            pass  # older JAX versions may not have this config key

        if output_mode is True:
            output_mode = 'normal'
        if isinstance(output_mode, str):
            output_mode = output_mode.strip().lower()
            if output_mode in ('false', 'none', 'off', 'no') or output_mode is None:
                output_mode = False
        if output_mode not in (None, False, 'minimal', 'normal', 'debug'):
            raise ValueError("Invalid output_mode. Use None, False, 'minimal', 'normal', or 'debug'.")

        self._output_mode = output_mode
        self._save_outputs = (self._output_mode is not False)
        self._debug_mode = (self._output_mode == 'debug')
        self._diagnostics = self._debug_mode
        self._skip_bad_init_eqdsks = skip_bad_init_eqdsks
        self._run_name = run_name
        self._relax_duration = float(relax_duration)
        _relax_dt = float(relax_dt)
        if _relax_dt <= 0:
            raise ValueError('relax_dt must be positive.')
        self._relax_dt = _relax_dt
        self._relax_kinetics = bool(relax_kinetics)
        self._relax = bool(relax)

        # Time-averaging settings for sawtooth smoothing
        self._t_ave_toggle = t_ave_toggle
        self._t_ave_window = t_ave_window
        self._t_ave_causal = t_ave_causal
        self._t_ave_ignore_start = t_ave_ignore_start

        dt_str = datetime.now().strftime('%Y-%m-%d_%H%M%S')
        _sim_start_time = time.time()

        self._run_timestamp = None if run_name == 'tmp' else dt_str
        self._output_file_tag = None if run_name == 'tmp' else f'{run_name}_{dt_str}'

        # ── Output directory ──
        if self._output_mode is not False:
            if run_name == 'tmp':
                self._out_dir = os.path.join('./TokaMaker_TORAX_outputs', 'tmp')
                if os.path.exists(self._out_dir):
                    shutil.rmtree(self._out_dir)
            else:
                self._out_dir = os.path.join('./TokaMaker_TORAX_outputs', self._output_file_tag)
            os.makedirs(self._out_dir, exist_ok=True)
            results_name = 'results.json'
            if self._output_file_tag is not None:
                results_name = f'{self._output_file_tag}_{results_name}'
            self._fname_out = os.path.join(self._out_dir, results_name)
        else:
            self._out_dir = None
            self._fname_out = None
            self._output_file_tag = None
            self._run_timestamp = None

        # ── Log file: saved in the output directory if there is one, otherwise
        #    in the script's run directory (cwd). ──
        log_dir = self._out_dir if self._out_dir is not None else '.'
        if run_name == 'tmp':
            self._log_file = os.path.abspath(os.path.join(log_dir, 'TokaMaker_TORAX_log_tmp.log'))
        else:
            self._log_file = os.path.abspath(os.path.join(log_dir, f'TokaMaker_TORAX_log_{run_name}_{dt_str}.log'))
        with open(self._log_file, 'w'):
            pass
        print(f'  Log file: {self._log_file}', flush=True)
        self._log(f'Log file: {self._log_file}')

        # In debug mode, attach file handler to Python logging so library
        # messages (TORAX, JAX, etc.) are captured in the log file.
        if self._debug_mode:
            self._logging_configured = False
            self.configure_redirect_to_log()

        # ── EQDSK directory: persisted only for debug mode ──
        if self._output_mode == 'debug':
            self._eqdsk_dir = self._out_dir
            self._eqdsk_dir_is_temp = False
        else:
            self._eqdsk_dir = tempfile.mkdtemp(prefix='TokaMaker_TORAX_equil_')
            self._eqdsk_dir_is_temp = True

        # ── Diverted / saddle-point configuration (set via set_x_points) ──
        if self._diverted_times is not None and self._x_point_targets is not None:
            t_div_start, t_div_end = self._diverted_times
            n_diverted = int(np.sum([(t_div_start <= t <= t_div_end) for t in self._tm_times]))
            self._log(f'Diverted window: t=[{t_div_start}, {t_div_end}] s '
                      f'({n_diverted}/{len(self._tm_times)} timesteps)')

        # ── Flattop detection ──
        Ip_arr = np.array(self._state['Ip'])
        Ip_max = np.max(Ip_arr)
        flattop_threshold = 0.95 * Ip_max
        above = Ip_arr >= flattop_threshold
        if np.any(above):
            ft_start = self._tm_times[np.argmax(above)]
            ft_end   = self._tm_times[len(above) - 1 - np.argmax(above[::-1])]
            self._flattop = np.array([(t >= ft_start and t <= ft_end) for t in self._tm_times])
        else:
            self._flattop = np.zeros(len(self._tm_times), dtype=bool)

        # ── Header ──
        self._print(f'\n{"="*60}\n TokaMaker_TORAX  \n run_name = {run_name} | t=[{self._t_init:.1f}, {self._t_final:.1f}] s '
                      f'| {len(self._tm_times)} timepoints | dt={self._tx_dt} s | max_loop={max_loop}')

        err = convergence_threshold + 1.0
        cflux_tx_prev = 0.0
        tm_cflux_psi = []
        tm_cflux_vloop = []
        tx_cflux_psi = []
        tx_cflux_vloop = []

        try:
            self._relax_mainrun_profile_history = []

            # ── Initial TORAX relax (optional) ──
            if initial_relax:
                self._print(f'\n{"="*60}\n  Initial TORAX relax\n{"="*60}')
                if self._relax_kinetics:
                    self._print('  Initial relax: relax_kinetics ON (evolve_* from set_evolve / config)')
                init_seed = self._init_files[0]
                if not self._test_eqdsk_tx_config(init_seed):
                    raise ValueError(f'Initial TORAX relax: first seed EQDSK not valid for TORAX: {init_seed}')
                self._run_tx_relax(stage='initial', eqdsk_path=init_seed, prescribed_profiles=None)
            else:
                self._psi_init = None
                self._n_e_init = None
                self._T_e_init = None
                self._T_i_init = None
                self._relax_profiles_snapshot = None

            self._current_loop = 0 if self._fly_loop0 else 1

            # ── Main loop ──
            # Counted loop indices: 1 ... max_loop (inclusive). Index 0 does not count toward max_loop.
            while err > convergence_threshold:
                if not (self._fly_loop0 and self._current_loop == 0):
                    if self._current_loop > max_loop:
                        break
                self._print(f'\n{"="*60}\n  Loop {self._current_loop}\n{"="*60}')

                _t_coupling_loop0 = time.perf_counter()
                try:
                    cflux_tx, cflux_tx_vloop = self._run_tx()
                except Exception as _tx_exc:
                    # On any TORAX failure (e.g. temperature collapse), save the scalar
                    # plot so the user has TM/TX time-series up to the failure point.
                    # TokaMaker entries plot as zeros / prior-loop values; that's expected.
                    # Save whenever fly() created an output directory (any non-False output_mode).
                    if self._out_dir is not None:
                        scalars_name = f'scalars_loop{self._current_loop:03d}_torax_failed.png'
                        if self._output_file_tag is not None:
                            scalars_name = f'{self._output_file_tag}_{scalars_name}'
                        _scalars_path = os.path.join(self._out_dir, scalars_name)
                        try:
                            plot_scalars(self, save_path=_scalars_path, display=False)
                            self._log(
                                f'TORAX failed at loop {self._current_loop} ({_tx_exc}); '
                                f'scalars plot saved to {_scalars_path}'
                            )
                            self._print(
                                '  TORAX failed; scalars plot saved to '
                                f'{_fmt_saved_artifact_link(_scalars_path)}'
                            )
                        except Exception as _e:
                            self._log(
                                f'plot_scalars failed after TORAX failure at loop '
                                f'{self._current_loop}: {_e}'
                            )
                            self._print(
                                f'  TORAX failed; scalars plot not saved: {_e}'
                            )
                        # Same layout as profile_plot at the last TORAX save (e.g. low-T stop).
                        try:
                            dt_fail = getattr(self, '_data_tree', None)
                            if (
                                dt_fail is not None
                                and hasattr(dt_fail, 'profiles')
                                and hasattr(dt_fail.profiles, 'psi')
                            ):
                                tx_times_fail = dt_fail.profiles.psi.coords['time'].values
                                if len(tx_times_fail) > 0:
                                    t_final_tx = float(tx_times_fail[-1])
                                    _tm_arr = np.asarray(self._tm_times, dtype=float)
                                    pp_tm = self._state.get('pp_prof_tm') or {}
                                    if pp_tm:
                                        i_prof = int(
                                            min(
                                                pp_tm.keys(),
                                                key=lambda idx: abs(float(_tm_arr[idx]) - t_final_tx),
                                            )
                                        )
                                    else:
                                        i_prof = int(np.argmin(np.abs(_tm_arr - t_final_tx)))
                                        if not _seed_tm_profiles_for_failure_profile_plot(self, i_prof):
                                            self._log(
                                                f'TORAX failed at loop {self._current_loop}: could not seed '
                                                f'TM profiles for failure profile plot (t_idx={i_prof}).'
                                            )
                                            self._print(
                                                '  TORAX failed; profile plot (last TORAX save) skipped '
                                                '(no TM profiles and no seed EQDSK profiles for this index).'
                                            )
                                            i_prof = None
                                    if i_prof is not None:
                                        self._tx_update(i_prof, dt_fail, tx_time=t_final_tx)
                                        prof_name = (
                                            f'profile_loop{self._current_loop:03d}_torax_failed_tfinal.png'
                                        )
                                        if self._output_file_tag is not None:
                                            prof_name = f'{self._output_file_tag}_{prof_name}'
                                        _prof_path = os.path.join(self._out_dir, prof_name)
                                        profile_plot(
                                            self, i_prof, t_final_tx,
                                            save_path=_prof_path, display=False,
                                        )
                                        self._log(
                                            f'TORAX failed at loop {self._current_loop} ({_tx_exc}); '
                                            f'profile plot (last TORAX time t={t_final_tx:.6g} s, TM idx {i_prof}) '
                                            f'saved to {_prof_path}'
                                        )
                                        self._print(
                                            '  TORAX failed; profile plot (last TORAX save) saved to '
                                            f'{_fmt_saved_artifact_link(_prof_path)}'
                                        )
                                else:
                                    self._log(
                                        f'TORAX failed at loop {self._current_loop}: no TORAX time coordinates '
                                        f'in data_tree; skipping failure profile plot.'
                                    )
                                    self._print(
                                        '  TORAX failed; profile plot (last TORAX save) skipped '
                                        '(no times in TORAX output).'
                                    )
                            else:
                                self._log(
                                    f'TORAX failed at loop {self._current_loop}: no partial data_tree; '
                                    f'skipping failure profile plot.'
                                )
                                self._print(
                                    '  TORAX failed; profile plot (last TORAX save) skipped '
                                    '(no partial TORAX data_tree).'
                                )
                        except Exception as _e:
                            self._log(
                                f'profile_plot after TORAX failure at loop {self._current_loop}: {_e}'
                            )
                            self._print(
                                f'  TORAX failed; profile plot (last TORAX save) not saved: {_e}'
                            )
                    raise

                cflux_tm, cflux_tm_vloop = self._run_tm()

                tm_cflux_psi.append(cflux_tm)
                tm_cflux_vloop.append(cflux_tm_vloop)
                tx_cflux_psi.append(cflux_tx)
                tx_cflux_vloop.append(cflux_tx_vloop)

                err = np.abs(cflux_tx - cflux_tx_prev) / cflux_tx_prev if cflux_tx_prev != 0 else convergence_threshold + 1.0
                cflux_diff = np.abs(cflux_tx - cflux_tm) / cflux_tm * 100.0 if cflux_tm != 0 else np.inf

                self._print(f'  Loop {self._current_loop} result: conv_err={err*100:.3f}% | '
                              f'TX-TM diff={cflux_diff:.4f}% | '
                              f'cflux_TX={cflux_tx:.4f} Wb | cflux_TM={cflux_tm:.4f} Wb')
                self._log(f'TX Convergence error = {err*100.0:.3f} %')
                self._log(f'Difference Convergence error = {cflux_diff:.4f} %')

                if self._output_mode in ('normal', 'minimal', 'debug') and self._out_dir is not None:
                    scalars_name = f'scalars_loop{self._current_loop:03d}.png'
                    if self._output_file_tag is not None:
                        scalars_name = f'{self._output_file_tag}_{scalars_name}'
                    _scalars_path = os.path.join(self._out_dir, scalars_name)
                    try:
                        plot_scalars(self, save_path=_scalars_path, display=False)
                    except Exception as _e:
                        self._log(f'plot_scalars failed at loop {self._current_loop}: {_e}')

                    plh_name = f'PLH_components_loop{self._current_loop:03d}.png'
                    if self._output_file_tag is not None:
                        plh_name = f'{self._output_file_tag}_{plh_name}'
                    try:
                        plot_PLH_components(self, save_path=os.path.join(self._out_dir, plh_name), display=False)
                    except Exception as _e:
                        self._log(f'plot_PLH_components failed at loop {self._current_loop}: {_e}')

                if self._output_mode == 'debug':
                    _loop_elapsed = time.perf_counter() - _t_coupling_loop0
                    _loop_mins, _loop_secs = divmod(_loop_elapsed, 60)
                    self._print(
                        f'  Loop {self._current_loop} wall time: '
                        f'{int(_loop_mins)}m {_loop_secs:.1f}s'
                    )

                cflux_tx_prev = cflux_tx
                self._current_loop += 1

        finally:
            # ── Cleanup temp EQDSK directory ──
            if getattr(self, '_eqdsk_dir_is_temp', False) and hasattr(self, '_eqdsk_dir') and os.path.exists(self._eqdsk_dir):
                try:
                    shutil.rmtree(self._eqdsk_dir)
                except OSError:
                    pass


        self._current_loop -= 1 # adjust back to last completed loop for reporting
        # ── End-of-run mode-specific outputs ──
        if self._output_mode is not False and not _in_jupyter() and self._out_dir is not None:
            profile_evo_name = 'profile_evolution.png'
            if self._output_file_tag is not None:
                profile_evo_name = f'{self._output_file_tag}_{profile_evo_name}'
            _profile_evo_path = os.path.join(self._out_dir, profile_evo_name)
            try:
                plot_profile_evolution(self, save_path=_profile_evo_path, display=False)
            except Exception as _e:
                self._log(f'plot_profile_evolution failed: {_e}')

            if self._output_mode in ('minimal', 'debug'):
                lcfs_evo_name = 'lcfs_evolution.png'
                if self._output_file_tag is not None:
                    lcfs_evo_name = f'{self._output_file_tag}_{lcfs_evo_name}'
                _lcfs_evo_path = os.path.join(self._out_dir, lcfs_evo_name)
                try:
                    plot_lcfs_evolution(self, save_path=_lcfs_evo_path, display=False)
                except Exception as _e:
                    self._log(f'plot_lcfs_evolution failed: {_e}')

            movie_name = f'movie_loop{self._current_loop:03d}.mp4'
            if self._output_file_tag is not None:
                movie_name = f'{self._output_file_tag}_{movie_name}'
            _movie_path = os.path.join(self._out_dir, movie_name)
            try:
                self.make_movie(save_path=_movie_path, display=False)
            except Exception as _e:
                self._log(f'make_movie failed: {_e}')

        # ── Summary table ──
        _sim_elapsed = time.time() - _sim_start_time
        n_completed = len(tx_cflux_psi)
        _coupling_loop0 = getattr(self, '_fly_loop0', True)
        _loop_label0 = 0 if _coupling_loop0 else 1
        converged = err <= convergence_threshold
        self._print(f'\n{"="*60}')
        if converged:
            self._print(f'  CONVERGED in {n_completed} loops (err={err*100:.3f}%)')
        else:
            self._print(f'  Max loop index ({max_loop}) reached (err={err*100:.3f}%)')

        # Print convergence history
        self._print(f'\n  {"Loop":<6} {"cflux TX [Wb]":<16} {"cflux TM [Wb]":<16} {"TX-TM diff %":<14}')
        self._print(f'  {"-"*52}')
        for s in range(len(tx_cflux_psi)):
            diff_pct = np.abs(tx_cflux_psi[s] - tm_cflux_psi[s]) / tm_cflux_psi[s] * 100 if tm_cflux_psi[s] != 0 else np.inf
            _idx = _loop_label0 + s
            self._print(f'  {_idx:<6} {tx_cflux_psi[s]:<16.4f} {tm_cflux_psi[s]:<16.4f} {diff_pct:<14.4f}')
        self._print(f'{"="*60}')

        # ── Elapsed time ──
        _mins, _secs = divmod(_sim_elapsed, 60)
        self._print(f'  Total sim time: {int(_mins)}m {_secs:.1f}s')

        if self._output_mode is not False:
            self._print(f'  Outputs saved to: {self._out_dir}')
        self._print(f'  Log file: {self._log_file}')

    # ─── Results & Visualization ────────────────────────────────────────────────

    @property
    def results(self):
        r'''! Access simulation results dict.'''
        return self._results
    
    @property
    def state(self):
        r'''! Access simulation state dict.'''
        return self._state

    def get_final_timepoint_results(self, eqdsk_save_dir=None):
        r'''! Profiles and scalars at the last TokaMaker timepoint of the last completed coupling loop.
        
                Call after fly() when TokaMaker has populated state['equil'] at every timestep index.
        
                Profiles use the usual TokaMaker_TORAX flux-surface dict form {'x', 'y', 'type'} (x =
                normalized poloidal flux when applicable). TokaMaker p_prime and FF_prime are
                P' and F F' from get_profiles(npsi=...) on the same grid as stored in state.
                Kinetic profiles and eta come from the TORAX-updated state at that timestep; current
                densities j_* are TORAX flux-surface profiles (sources for the GS solve).
        
                Scalars: fusion Q from TORAX (last save or live data tree); q95 and q0 from
                TokaMaker get_q (l_i likewise from get_stats); V_loop from TORAX
                v_loop_lcfs; Ip is the TokaMaker equilibrium value.
        
                @param eqdsk_save_dir If set (non-empty str or path-like), write the final TokaMaker
                       equilibrium gEQDSK with save_eqdsk into this directory (created if needed).
                       If None (default), no file is written.
                @return dict with time, simulation_loop, tm_time_index, profile keys, nested
                        j current profiles, scalar values, coil_currents, coil_currents_by_region,
                        tokamaker_equilibrium (TokaMaker equilibrium instance at this timestep — same
                        reference as state['equil'][tm_time_index]), and eqdsk_path (absolute path
                        string when saved, else None).
                
        '''
        def _snap(prof):
            if prof is None:
                return None
            out = {}
            for k, v in prof.items():
                if isinstance(v, np.ndarray):
                    out[k] = np.array(v, copy=True)
                else:
                    out[k] = copy.deepcopy(v)
            return out

        tm_times = getattr(self, '_tm_times', None)
        if tm_times is None or len(tm_times) == 0:
            raise RuntimeError('TokaMaker_TORAX has no tm_times; run fly() first.')
        i = len(tm_times) - 1
        t = float(tm_times[i])
        equil = self._state.get('equil', {}).get(i)
        if equil is None:
            raise RuntimeError(
                f'No TokaMaker equilibrium at final time index {i} (t={t} s); cannot build snapshot.'
            )

        coils, currents_reg = equil.get_coil_currents()

        eqdsk_path = None
        if eqdsk_save_dir:
            out_dir = os.path.abspath(os.path.expanduser(str(eqdsk_save_dir)))
            os.makedirs(out_dir, exist_ok=True)
            lp = int(getattr(self, '_current_loop', -1))
            fname = f'TokaMaker_TORAX_final_loop{lp:03d}_t{t:.6f}s.eqdsk'
            eqdsk_path = os.path.join(out_dir, fname)
            nr_nz = EQDSK_SAVE_NR_NZ_SEQUENCE[-1]
            with self._quiet_tm():
                equil.save_eqdsk(
                    eqdsk_path,
                    lcfs_pad=1.0 - self._last_surface_factor,
                    run_info='TokaMaker EQDSK (TokaMaker_TORAX final timepoint)',
                    cocos=self._cocos,
                    nr=nr_nz,
                    nz=nr_nz,
                    truncate_eq=self._truncate_eq,
                )

        Q_val = np.nan
        Q_ts = self._results.get('Q')
        if isinstance(Q_ts, dict) and Q_ts.get('x') is not None and len(Q_ts['x']) > 0:
            xq = np.asarray(Q_ts['x'], dtype=float)
            yq = np.asarray(Q_ts['y'], dtype=float)
            Q_val = float(yq[int(np.argmin(np.abs(xq - t)))])
        else:
            dt = getattr(self, '_data_tree', None)
            if dt is not None:
                try:
                    Q_val = float(self._extract_tx_scalar(dt, 'Q_fusion', t))
                except Exception:
                    Q_val = np.nan

        j_keys = [
            ('j_total', 'j_tot'),
            ('j_ohmic', 'j_ohmic'),
            ('j_non_inductive', 'j_ni'),
            ('j_bootstrap', 'j_bootstrap'),
            ('j_ecrh', 'j_ecrh'),
            ('j_external', 'j_external'),
            ('j_generic_current', 'j_generic_current'),
        ]
        j_out = {}
        for out_name, state_key in j_keys:
            j_out[out_name] = _snap(self._state.get(state_key, {}).get(i))

        return {
            'time': t,
            'simulation_loop': int(getattr(self, '_current_loop', -1)),
            'tm_time_index': i,
            'ne': _snap(self._state.get('n_e', {}).get(i)),
            'Te': _snap(self._state.get('T_e', {}).get(i)),
            'Ti': _snap(self._state.get('T_i', {}).get(i)),
            'q': _snap(self._state.get('q_prof_tm', {}).get(i)),
            'psi': _snap(self._state.get('psi_tm', {}).get(i)),
            'q_torax': _snap(self._state.get('q_prof_tx', {}).get(i)),
            'psi_torax': _snap(self._state.get('psi_tx', {}).get(i)),
            'eta': _snap(self._state.get('eta_prof', {}).get(i)),
            'p_prime': _snap(self._state.get('pp_prof_tm', {}).get(i)),
            'FF_prime': _snap(self._state.get('ffp_prof_tm', {}).get(i)),
            'j': j_out,
            'Q': Q_val,
            'q95': float(self._state['q95_tm'][i]),
            'q0': float(self._state['q0_tm'][i]),
            'l_i': float(self._state['l_i_tm'][i]),
            'V_loop': float(self._state['vloop_tx'][i]),
            'Ip': float(self._state['Ip_tm'][i]),
            'Ip_torax': float(self._state['Ip_tx'][i]),
            'coil_currents': copy.deepcopy(coils),
            'coil_currents_by_region': np.array(currents_reg, copy=True),
            'tm_equilibrium': equil,
            'eqdsk_path': eqdsk_path,
        }

    # ─── Visualization ──────────────────────────────────────────────────────────
    # Creates TokaMaker_TORAX methods that wrap standalone plotting scripts below.

    def make_movie(self, save_path=None, **kwargs):
        r'''! Generate pulse movie from stored psi snapshots.
                @param save_path Path to save MP4 file. If None, does not save.
                
        '''
        return make_movie(self, save_path=save_path, **kwargs)

    def plot_scalars(self, save_path=None, display=True, **kwargs):
        r'''! Plot scalar time traces (Ip, Q, Te, ne, power channels, etc.).
                @param save_path Path to save figure. If None, does not save.
                @param display Whether to show the plot.
                
        '''
        return plot_scalars(self, save_path=save_path, display=display, **kwargs)

    def plot_profiles(self, **kwargs):
        r'''! Interactive profile viewer (ipywidgets slider in Jupyter, static otherwise).'''
        return plot_profiles_interactive(self, **kwargs)

    def plot_profile_evolution(self, save_path=None, display=True, one_plot=False, **kwargs):
        r'''! Plot profile evolution over time.
                @param save_path Path to save figure. If None, does not save.
                @param display Whether to show the plot.
                @param one_plot If True, combine all pulse phases into one figure.
                
        '''
        return plot_profile_evolution(self, save_path=save_path, display=display, one_plot=one_plot, **kwargs)

    def plot_coils(self, save_path=None, display=True, **kwargs):
        r'''! Plot coil current traces over the pulse.
                @param save_path Path to save figure. If None, does not save.
                @param display Whether to show the plot.
                
        '''
        return plot_coils(self, save_path=save_path, display=display, **kwargs)

    def plot_lcfs_evolution(self, save_path=None, display=True, one_plot=False, **kwargs):
        r'''! Plot time evolution of the LCFS for each pulse phase (rampup, flattop, rampdown).
                @param save_path Path prefix to save figures. If None, does not save.
                @param display Whether to show the plots.
                @param one_plot If True, combine all pulse phases into one figure.
                
        '''
        return plot_lcfs_evolution(self, save_path=save_path, display=display, one_plot=one_plot, **kwargs)

    def summary(self, **kwargs):
        r'''! Print/display a physics summary of the simulation.'''
        return summary(self, **kwargs)


# =============================================================================
#  Visualization
# =============================================================================

# ── Style constants ──────────────────────────────────────────────────────────

# Formatting used throughout plots for consistency
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
MOVIE_EQUIL_PSI_PLASMA_NLEVELS = 16


# ── Helpers ──────────────────────────────────────────────────────────────────

# Number of rho points the imposed manual pedestal shape is sampled onto in
# [ped_rho, 1]. TORAX linearly interpolates between them, so the IBC follows the shape.
PED_N_SAMPLE = 24

# rho_norm span (core-side of the IBC inner edge) over which the evolved gradient is
# least-squares fit for inner-edge smoothing. Wider = smoother slope, esp. for density.
PED_GRAD_WINDOW = 0.03


def _mtanh_pedestal(rho, ped_top, foot, core_slope, knee, hwid, blend):
    r'''! Linear core gradient smoothly blended into a tanh pedestal.

        A line of slope core_slope through (ped_rho, ped_top) owns the core; a
        tanh dropping from ped_top to foot owns the edge; a sigmoid (centred at
        knee, width blend) hands off between them. Because the sigmoid weight is
        ~0 at the inner edge, the gradient there is exactly core_slope.

        @param rho Normalized rho array.
        @param ped_top Value at the inner (left) edge.
        @param foot Separatrix value at rho=1.
        @param core_slope d(value)/d(rho_norm) imposed at the inner edge.
        @param knee Rho where the core hands off to the pedestal.
        @param hwid Pedestal (tanh) half-width.
        @param blend Sharpness of the core->pedestal handoff.
    '''
    A = (ped_top - foot) / 2.0
    B = (ped_top + foot) / 2.0
    ped = A * np.tanh((knee - rho) / hwid) + B
    core_line = ped_top + core_slope * (rho - rho[0])
    w = 1.0 / (1.0 + np.exp(-(rho - knee) / blend))
    return (1.0 - w) * core_line + w * ped


def _bc_at(right_bc, t):
    r'''! Evaluate an edge (rho=1) BC at time t: a scalar as-is, or a {time: value} map
        linearly interpolated in time (constant-extrapolated past the ends).'''
    if not isinstance(right_bc, dict):
        return float(right_bc or 0.0)
    ts = sorted(right_bc)
    return float(np.interp(t, ts, [right_bc[k] for k in ts]))


def _ped_top_slope_at(shape_by_time, t, default):
    r'''! Interpolate (ped_top, slope) in time from a {tm_time: (ped_top, slope)} dict
        (constant-extrapolated past the ends). Returns `default` if the dict is empty.'''
    if not shape_by_time:
        return default
    ts = sorted(shape_by_time)
    tops = [shape_by_time[k][0] for k in ts]
    slopes = [shape_by_time[k][1] for k in ts]
    return float(np.interp(t, ts, tops)), float(np.interp(t, ts, slopes))


def _pedestal_ibc(ped_rho, right_bc, shape_by_time, default_shape, transition, windows,
                  knee, hwid, blend):
    r'''! Build a TORAX internal-boundary-condition SparseTimeVaryingArray for one
        field: {time: {(ped_rho, 1.0): {rho: value}}}.

        Three band states:
          - OFF: all-zero -> ignored by TORAX (L-mode, TORAX evolves freely to the edge).
          - L-MODE EDGE: flat at the edge BC across the whole band -> agrees with
            n_e_right_bc at rho=1 and imposes no interior gradient.
          - ON: the mtanh pedestal, foot = edge BC at that time.
        The same (ped_rho, 1.0) location is present at every listed time (the "zero-off
        trick") so TORAX does not constant-extrapolate the band on backwards.

        PER-TM-TIME shape: the mtanh ped_top and inner-edge slope are taken from
        `shape_by_time` ({tm_time: (ped_top, slope)}, from the previous loop's evolved
        profile at each tm_time) and the band is listed at EVERY tm_time in the H-mode
        window, so the imposed pedestal tracks the evolving profile time-point by
        time-point. (A single time-averaged shape is wrong: the few transition-region
        tm_times have anomalous steep slopes that poison the average.) Before any
        smoothing exists, `default_shape` (ped_top, slope) is used at all times.

        The transition ramps L-MODE EDGE -> ON (NOT OFF -> ON), so the foot stays at the
        edge BC throughout the ramp (only the pedestal height above the foot grows in);
        ramping from zeros instead would spike rho=1 and notch just inside.

        @param ped_rho Inner edge rho of the imposed pedestal.
        @param right_bc Edge (rho=1) BC: scalar or {time: value} map; sets the foot per time.
        @param shape_by_time {tm_time: (ped_top, slope)} per-time mtanh shape (may be empty).
        @param default_shape (ped_top, slope) used when shape_by_time is empty / outside it.
        @param transition Ramp duration [s] at each L<->H edge.
        @param windows List of (lh_time, hl_time) H-mode windows; hl_time may be None (stays on).
        @param knee/hwid/blend mtanh shape parameters.
    '''
    rhos = np.linspace(ped_rho, 1.0, PED_N_SAMPLE)
    loc = (ped_rho, 1.0)
    off = {float(r): 0.0 for r in rhos}
    # Times at which to (re)list the ON shape inside the window: every smoothed tm_time,
    # plus any edge-BC breakpoints (so the foot tracks a ramping right_bc).
    shape_times = sorted(shape_by_time)
    bc_times = sorted(right_bc) if isinstance(right_bc, dict) else []

    def lmode_edge(t):
        # Flat at the edge BC across the band: agrees with n_e_right_bc at rho=1, no gradient.
        return {float(r): _bc_at(right_bc, t) for r in rhos}

    def on(t):
        ped_top, core_slope = _ped_top_slope_at(shape_by_time, t, default_shape)
        foot = _bc_at(right_bc, t)
        return {float(r): float(v) for r, v in
                zip(rhos, _mtanh_pedestal(rhos, ped_top, foot, core_slope, knee, hwid, blend))}

    eps = 1.0e-6   # tiny step so OFF->flat-edge is abrupt (band genuinely OFF through L-mode)
    spec = {0.0: {loc: off}}
    for lh, hl in windows:
        on_start = lh + transition
        on_end = hl if hl is not None else float('inf')
        spec[lh - eps] = {loc: off}                  # OFF right up to lh: TORAX free in L-mode
        spec[lh] = {loc: lmode_edge(lh)}             # flat at edge BC, then ramp UP to the pedestal
        spec[on_start] = {loc: on(on_start)}
        # List the ON shape at every smoothed tm_time and edge-BC breakpoint inside the window,
        # so the imposed pedestal tracks the per-time evolved shape and the ramping edge BC.
        for tt in sorted(set(shape_times) | set(bc_times)):
            if on_start < tt < on_end:
                spec[tt] = {loc: on(tt)}
        if hl is not None:
            spec[hl] = {loc: on(hl)}
            spec[hl + transition] = {loc: lmode_edge(hl + transition)}   # ramp DOWN to flat edge BC
            spec[hl + transition + eps] = {loc: off}                     # then OFF: TORAX free again
    return spec


def _core_transition_profile(rhos, ped_rho, ped, core, exp_a=2.0, exp_b=2.0):
    r'''! Smooth CORE-only H-mode-like profile over [0, ped_rho] that joins the pedestal top
        WITHOUT its own edge roll-off (so it doesn't create a second pedestal):

            f(rho) = ped + (core - ped) * (1 - (rho/ped_rho)^a)^b

        With a>=2, b>=2 this has ~zero gradient at rho=0 (flat core), monotonically falls through
        mid-radius, and FLATTENS into ped at rho=ped_rho (zero gradient there too), so it hands off
        smoothly to the pedestal band (whose inner-edge slope is ~0). Unlike Hmode_profiles it has
        NO tanh edge inside the core grid -- the drop ends exactly at ped_rho = ped_height, so the
        only pedestal is the [ped_rho,1] band. Looks like the H-mode core (steepest mid-radius,
        gentle at both ends).

        @param rhos Normalized rho array over [0, ped_rho] to sample onto.
        @param ped_rho Inner pedestal edge (where the core meets the pedestal top).
        @param ped Pedestal-top value (= ped_height); the core endpoint at ped_rho.
        @param core On-axis (rho=0) value.
        @param exp_a/exp_b Core-shape exponents (a>=2, b>=2 for flat ends). Default 2/2.
        @return core profile evaluated on `rhos`.
    '''
    x = np.clip(np.asarray(rhos, dtype=float) / ped_rho, 0.0, 1.0)
    return ped + (core - ped) * (1.0 - x ** exp_a) ** exp_b


def _merge_ibc_specs(spec_a, spec_b):
    r'''! Merge two single-field IBC SparseTimeVaryingArrays that use DIFFERENT band locations
        (e.g. the pedestal band (ped_rho,1) and the core band (0,core_hi)) into one
        {time: {(lo,hi): {rho: value}}}.

        TORAX groups IBC entries BY LOCATION (interpolated_param_2d): each band (lo,hi) becomes its
        own TimeVaryingArray from only the times IT was listed at, interpolated independently of the
        other band. So we just take the per-time union of whatever each spec actually lists -- a band
        absent at some time is NOT carried forward here (it doesn't need to be; TORAX interpolates it
        from its own listed times). This keeps each band sparse: e.g. the core band stays at its ~5
        natural times even when the pedestal band is listed at every H-mode tm_time. Listing the
        absent band at every union time would bloat the spec with redundant all-zero entries and
        slow TORAX, with no effect on the result.

        @param spec_a/spec_b {time: {(lo,hi): {rho: value}}} specs with distinct band keys.
        @return merged {time: {(lo,hi): {rho: value}}}.
    '''
    merged = {}
    for t in sorted(set(spec_a) | set(spec_b)):
        entry = {}
        if t in spec_a:
            entry.update(spec_a[t])
        if t in spec_b:
            entry.update(spec_b[t])
        merged[t] = entry
    return merged


def _transition_ibc(loc, lmode_profile, hmode_target, lh, transition):
    r'''! Build a CORE (rho 0 -> ped_rho) internal-boundary-condition SparseTimeVaryingArray for
        one field, active ONLY across the L->H transition window [lh, lh+transition]:
        {time: {loc: {rho: value}}} where loc = (0.0, ped_rho).

        The band ramps (TORAX linearly interpolates in time) from the TORAX-evolved L-mode core
        at t=lh to the H-mode core at t=lh+transition, then retreats (OFF). The pedestal band owns
        [ped_rho,1] at all times, so the two bands meet at ped_rho with no overlap. NO H->L entry.

        @param loc Band location (0.0, ped_rho).
        @param lmode_profile {rho: value} evolved L-mode core shape at ~lh (start of the ramp).
        @param hmode_target {rho: value} H-mode core profile (end of the ramp).
        @param lh L->H transition time [s].
        @param transition Ramp duration [s].
        @return {time: {loc: {rho: value}}} spec.
    '''
    off = {r: 0.0 for r in lmode_profile}
    eps = 1.0e-6
    return {
        0.0:                    {loc: off},          # OFF from t=0 (zero-off trick anchor)
        lh - eps:               {loc: off},          # OFF right up to lh: TORAX free in L-mode
        lh:                     {loc: lmode_profile}, # ramp starts at the evolved L-mode core
        lh + transition:        {loc: hmode_target},  # ramp ends at the H-mode core
        lh + transition + eps:  {loc: off},           # band retreats; pedestal band stays on
    }


def _in_jupyter():
    r'''! Return True if running inside a Jupyter notebook.'''
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        return shell in ('ZMQInteractiveShell',)
    except Exception:
        return False


def _save_or_display(fig, save_path=None, display=True):
    r'''! Save figure and/or display it, then close.'''
    if save_path is not None:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    if display:
        plt.show()
    else:
        plt.close(fig)


def _style(ax):
    r'''! Apply consistent styling to a matplotlib Axes.'''
    ax.grid(True, alpha=GRID_ALPHA)
    ax.tick_params(labelsize=TICK_FS)


def _tx_scalar(tt, var_name, scale=1.0):
    r'''! Return (times, values) arrays from data_tree.scalars at full TORAX resolution.'''
    dt = getattr(tt, '_data_tree', None)
    if dt is None:
        return None, None
    if not hasattr(dt.scalars, var_name):
        return None, None
    var = getattr(dt.scalars, var_name)
    return var.coords['time'].values, var.to_numpy() * scale


def _tx_profile_at_rho(tt, var_name, rho_val, rho_coord='rho_norm', scale=1.0):
    r'''! Return (times, values) for a profile variable at a fixed rho, full resolution.'''
    dt = getattr(tt, '_data_tree', None)
    if dt is None:
        return None, None
    var = getattr(dt.profiles, var_name)
    sliced = var.sel(**{rho_coord: rho_val}, method='nearest')
    return sliced.coords['time'].values, sliced.to_numpy() * scale


def _prof(state_dict, idx):
    r'''! Return (x, y) arrays from a profile state dict, or (None, None).'''
    d = state_dict.get(idx)
    if d is None:
        return None, None
    return d['x'], d['y']


def _make_temp_dir_viz():
    r'''! Create temp directory for visualization temporary files: RAM-backed on Linux, OS default elsewhere.'''
    if platform.system() == 'Linux' and os.path.isdir('/dev/shm'):
        return tempfile.mkdtemp(prefix='TokaMaker_TORAX_viz_', dir='/dev/shm')
    return tempfile.mkdtemp(prefix='TokaMaker_TORAX_viz_')


def _fmt_saved_artifact_link(path):
    r'''! Absolute path as a file:// URI when possible (terminal / IDE hyperlink), else abs path.'''
    abs_path = os.path.abspath(path)
    try:
        return Path(abs_path).as_uri()
    except ValueError:
        return abs_path


def _seed_tm_profiles_for_failure_profile_plot(tt, i):
    r'''! Fill TM comparison fields for profile_plot from seed EQDSK/GS inputs when TM has not run at index i.'''
    s = tt._state
    if i in s.get('pp_prof_tm', {}):
        return True
    if i not in s.get('pp_prof', {}):
        return False
    s.setdefault('pp_prof_tm', {})
    s.setdefault('ffp_prof_tm', {})
    s.setdefault('p_prof_tm', {})
    s.setdefault('psi_tm', {})
    s.setdefault('q_prof_tm', {})
    s['pp_prof_tm'][i] = copy.deepcopy(s['pp_prof'][i])
    s['ffp_prof_tm'][i] = copy.deepcopy(s['ffp_prof'][i])
    p_axis = float(s['pax'][i])
    p_tm = copy.deepcopy(s['p_prof_eqdsk'][i])
    p_tm['y'] = np.asarray(p_tm['y'], dtype=float) * max(p_axis, 1e-300)
    s['p_prof_tm'][i] = p_tm
    psi_a = float(s['psi_axis_tm'][i])
    psi_l = float(s['psi_lcfs_tm'][i])
    s['psi_tm'][i] = {
        'x': tt._psi_N.copy(),
        'y': psi_a + (psi_l - psi_a) * tt._psi_N,
        'type': 'linterp',
    }
    s['q_prof_tm'][i] = copy.deepcopy(s['q_prof_eqdsk'][i])
    return True


def _x_points_active(tt, i, t=None):
    r'''! Return True when X-point targets should be applied at timestep index i.'''
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


# ── Profile plot (per-timestep diagnostic) ───────────────────────────────────

def profile_plot(tt, i, t, save_path=None, display=True):
    r'''! Detailed profile comparison at a single timestep.'''
    s = tt._state
    psi_N = tt._psi_N

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
    ax.axhline(1.0, color='k', ls='-', lw=1, label='q=1')
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
    handles, _labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=9)
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
    handles, _labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, ncol=2)
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
    handles, _labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=9)
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
    handles, _labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=9)
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
    handles, _labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=9)
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
    handles, _labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, hspace=0.3, wspace=0.35)
    _save_or_display(fig, save_path, display)


def IBC_debug(tt, i, t, save_path=None):
    r'''! Debug figure for the internal_manual IBC technique (debug mode, one per TM solve).
            Always saved to save_path, never displayed. Three rows (T_e, T_i, n_e), all in rho_norm,
            zoomed on the pedestal region (or the whole domain during the L->H transition):
              - the TORAX-evolved profile (solid),
              - the IBC target (soft) imposed at this time (x markers),
              - the edge BC point at rho=1 (square),
              - the inner-edge slope fit from THIS loop's evolved profile (dashed line over the
                fit window),
              - the slope used to BUILD this loop's IBC, i.e. last loop's smoothed slope (dotted
                line extending beyond the fit window),
              - during the L->H transition: the core IBC target (rho 0->ped_rho) and core target point,
              - text of both slopes, the ped-top value, and the IBC stiffness prefactor.
            NOTE: TORAX applies IBCs as a SOFT adaptive source/sink with stiffness
            adaptive_{T,n}_source_prefactor, not a hard Dirichlet BC. The evolved profile relaxes
            toward the target with that stiffness, so a target<->evolved gap (esp. for T, whose
            default stiffness is weaker relative to heat transport) is expected behaviour, not a
            bug; raise T_source_prefactor in set_pedestal to track the target more exactly.
            Shows nothing useful unless the IBC is active (timing known); harmless otherwise.
    '''
    snap = getattr(tt, '_ped_ibc_snapshot', None)
    dt = getattr(tt, '_data_tree', None)
    if snap is None or dt is None:
        return
    ped_rho = snap['ped_rho']
    ped_x0 = max(0.0, ped_rho - 0.10)   # zoom: a bit inside the inner edge out to rho=1
    fields = [('T_e', 'T_e [keV]', 1.0), ('T_i', 'T_i [keV]', 1.0),
              ('n_e', r'n_e [m$^{-3}$]', 1.0)]
    # During the L->H transition the full-domain IBC is active; zoom out to the whole domain so the
    # imposed core profile + constraints are visible. Detect it from any field's transition band.
    lh = (snap['windows'][0][0] if snap.get('windows') else None)
    in_transition = (lh is not None and lh <= t <= lh + snap['transition']
                     and any(f.get('transition_active') for f in snap['fields'].values()))
    x0 = 0.0 if in_transition else ped_x0

    # IBCs are SOFT penalties (TORAX adaptive_{T,n}_source_prefactor): the evolved profile relaxes
    # toward the imposed target with this stiffness, so a gap target<->evolved is expected, not a bug.
    fig, axes = plt.subplots(3, 1, figsize=(10, 13), sharex=True)
    plt.suptitle(f'Manual pedestal IBC debug (soft target) — loop {tt._current_loop} - '
                 f't-idx {i}/{len(tt._tm_times)-1} - t = {t:.2f} s\n'
                 f'IBC stiffness: T_source_prefactor={snap["T_pref"]:.2g}, '
                 f'n_source_prefactor={snap["n_pref"]:.2g}', fontsize=12)

    for ax, (field, ylabel, scale) in zip(axes, fields):
        fcfg = snap['fields'].get(field)
        pref = snap['n_pref'] if field == 'n_e' else snap['T_pref']
        da = getattr(dt.profiles, field, None)
        if da is not None:
            rho = da.coords['rho_norm'].values
            y = da.sel(time=t, method='nearest').values * scale
            m = rho >= x0
            ax.plot(rho[m], y[m], 'k-', lw=1.8, label='TORAX evolved (relaxes to target)')
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.axvline(ped_rho, color='grey', ls=':', alpha=0.7, label=f'inner edge (rho={ped_rho:.3f})')

        if fcfg is None:
            ax.set_title(f'{field}: no pedestal height set')
            ax.legend(fontsize=8)
            continue

        # Per-time IBC shape (ped_top, slope) interpolated at THIS time t (same as the IBC build).
        used_top, used_slope = _ped_top_slope_at(fcfg['shape_by_time'], t, fcfg['default_shape'])
        # IBC sample points imposed at this time (rebuild the same mtanh samples).
        spec = _pedestal_ibc(ped_rho, fcfg['right_bc'], fcfg['shape_by_time'], fcfg['default_shape'],
                             snap['transition'], snap['windows'],
                             snap['knee'], snap['hwid'], snap['blend'])
        # Evaluate the spec at time t (linear-in-time per location, same as TORAX).
        ibc_pts = _eval_ibc_at_time(spec, t)
        if ibc_pts is not None:
            rr, vv = ibc_pts
            ax.plot(rr, np.array(vv) * scale, 'x', color='tab:orange', ms=7, mew=1.6,
                    label='IBC target (soft)')
        # Edge BC point at rho=1.
        bc = _bc_at(fcfg['right_bc'], t) * scale
        ax.plot(1.0, bc, 's', color='tab:red', ms=8, label='edge BC (rho=1)')

        # Core L->H transition: overlay the imposed core IBC target (rho 0 -> ped_rho) + constraints.
        if fcfg.get('transition_active'):
            trans_spec = _transition_ibc(fcfg['trans_loc'], fcfg['lmode_profile'],
                                         fcfg['hmode_target'], lh, snap['transition'])
            tpts = _eval_ibc_at_time(trans_spec, t)
            if tpts is not None:   # band active (ramping) at this time
                rr, vv = tpts
                ax.plot(rr, np.array(vv) * scale, '.', color='tab:purple', ms=6,
                        label='transition IBC target (core, soft)')
                core_src = fcfg.get('core_source', 'user')
                ax.plot(0.0, fcfg['core'] * scale, 'D', color='tab:purple', ms=8,
                        label=f'core target = {fcfg["core"]:.3g} ({core_src})')

        # Slope used to BUILD this IBC at THIS time (prev-loop smoothed slope, interpolated to t),
        # as a line through (ped_rho, ped_top) extending beyond the fit window so it's visible.
        xs = np.array([ped_rho - 0.06, ped_rho + 0.04])
        ax.plot(xs, (used_top + used_slope * (xs - ped_rho)) * scale,
                ls=':', color='tab:blue', lw=1.6,
                label=f'slope used (prev loop) = {used_slope:.3g}')

        # Slope freshly fit from THIS loop's evolved profile near the inner edge (what next loop
        # will use at this time), over the fit window only. Should match 'slope used' once converged.
        if da is not None:
            hi = int(np.argmin(np.abs(rho - ped_rho)))
            lo = min(int(np.searchsorted(rho, ped_rho - PED_GRAD_WINDOW)), hi - 2)
            if lo >= 0 and hi > lo:
                fit_slope = float(np.polyfit(rho[lo:hi + 1], y[lo:hi + 1], 1)[0])
                xf = rho[lo:hi + 1]
                yf = y[hi] + fit_slope * (xf - rho[hi])
                ax.plot(xf, yf, '--', color='tab:green', lw=2.0,
                        label=f'slope fit (this loop) = {fit_slope:.3g}')

        tag = 'smoothed' if fcfg['smoothed'] else 'user-height fallback'
        ax.set_title(f'{field}:  ped_top={used_top:.3g} ({tag}),  '
                     f'edge BC={bc:.3g},  stiffness={pref:.2g}', fontsize=10)
        ax.legend(fontsize=8, loc='best')

    axes[-1].set_xlabel(r'$\rho_\mathrm{norm}$')
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.18)
    # Always save, never display (this figure is only written during the sim).
    _save_or_display(fig, save_path, display=False)


def _eval_ibc_at_time(spec, t):
    r'''! Evaluate a single-field IBC SparseTimeVaryingArray spec at time t -> (rhos, values),
            linearly interpolating each rho's value in time between listed times (matching TORAX).
            Returns None if the band is OFF (all zero) at t.'''
    loc = next(iter(next(iter(spec.values()))))   # the (lo, hi) location key
    times = sorted(spec)
    lo_t = max([x for x in times if x <= t], default=times[0])
    hi_t = min([x for x in times if x >= t], default=times[-1])
    dlo = spec[lo_t][loc]
    dhi = spec[hi_t][loc]
    w = 0.0 if hi_t == lo_t else (t - lo_t) / (hi_t - lo_t)
    rhos = sorted(dlo)
    vals = [(1.0 - w) * dlo[r] + w * dhi[r] for r in rhos]
    if all(v == 0.0 for v in vals):
        return None
    return np.array(rhos), vals


# ── TokaMaker diagnostic plot ─────────────────────────────────────────────────

def tm_diagnostic_plot(tt, i, t, level_attempts, solve_succeeded, save_path=None, display=True,
                       tm_gs_ok=None, step_error=None):
    r'''! TokaMaker input/output diagnostic plot for a single timestep.

        @param solve_succeeded Whether the full TM timestep succeeded including EQDSK validation
               for TORAX (same flag as after `_run_tm` updates).
        @param tm_gs_ok Whether the Grad-Shafranov solve succeeded (independent of EQDSK).
               If None, inferred from level_attempts.
        @param step_error Optional failure string (e.g. TORAX EQDSK rejection); overrides parsing
               the last level attempt when provided.
    '''
    s = tt._state

    _winning = next((a for a in level_attempts if a['succeeded']), None)
    _last = level_attempts[-1] if level_attempts else {}
    _gs_ok = tm_gs_ok if tm_gs_ok is not None else (_winning is not None)
    fail_msg = (
        step_error if step_error is not None
        else (_last.get('error') if not solve_succeeded else None)
    )

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

    fig = plt.figure(figsize=(22, 12))
    gs_layout = fig.add_gridspec(3, 6, hspace=0.70, wspace=0.55)

    ax_ffp_tx = fig.add_subplot(gs_layout[0, 0:2])
    ax_pp_tx = fig.add_subplot(gs_layout[1, 0:2])
    ax_eta = fig.add_subplot(gs_layout[2, 0:2])
    ax_ffp_tm = fig.add_subplot(gs_layout[0, 2:4])
    ax_pp_tm = fig.add_subplot(gs_layout[1, 2:4])
    ax_q_tm = fig.add_subplot(gs_layout[2, 2:4])
    ax_tbl1 = fig.add_subplot(gs_layout[0, 4:6])
    gs_tbl_eq = gs_layout[1:3, 4:6].subgridspec(1, 2, width_ratios=[1.15, 0.85], wspace=0.12)
    ax_tbl2 = fig.add_subplot(gs_tbl_eq[0, 0])
    ax_mini_eq = fig.add_subplot(gs_tbl_eq[0, 1])

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
        _beta_pol_tm_ok_str = '\u2014'
        try:
            _equil_snap_ok = s.get('equil', {}).get(i)
            if _equil_snap_ok is not None:
                _ok_stats = _equil_snap_ok.get_stats()
                if 'beta_pol' in _ok_stats:
                    _beta_pol_tm_ok_str = f"{float(_ok_stats['beta_pol']) / 100.0:.4f}"
        except Exception:
            pass

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
            ['beta_pol', '\u2014', f'{beta_pol_tx:.4f}', _beta_pol_tm_ok_str],
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

        ax_q_tm.axhline(1.0, color='k', ls='-', lw=1, label='q=1')
        try:
            psi_geo, q_tm_vals, _, _, _, _ = tt._tm.get_q(npsi=len(tt._psi_N), psi_pad=1-tt._last_surface_factor)
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
        _tm_nl_tol = getattr(getattr(tt._tm, 'settings', None), 'nl_tol', None)
        _tol_str = f'{_tm_nl_tol:.2e}' if _tm_nl_tol is not None else '\u2014'
        _raw_err = _last.get('error') or (fail_msg if fail_msg else 'Unknown')
        # Strip the generic "Error in solve: " prefix so only the specific reason shows.
        _fail_reason_str = _raw_err.removeprefix('Error in solve: ') if _raw_err else 'Unknown'

        _beta_pol_tm_str = '\u2014'
        try:
            _tm_stats = tt._tm.get_stats()
            if 'beta_pol' in _tm_stats:
                _beta_pol_tm_str = f"{float(_tm_stats['beta_pol']) / 100.0:.4f}"
        except Exception:
            pass

        _psi_lcfs_tm_achieved_str = '\u2014'
        try:
            _pb = tt._tm.psi_bounds  # [psi_axis, psi_lcfs] in TM-native sign
            _psi_lcfs_tm_achieved_str = f'{float(_pb[1]):.4f} Wb/rad'
        except Exception:
            pass

        input_rows = [
            ['Parameter', 'Init EQDSK', 'TORAX'],
            ['Ip', f'{Ip_seed / 1e6:.3f} MA', f'{Ip_tx / 1e6:.3f} MA'],
            ['pax', f'{pax_seed / 1e3:.2f} kPa', f'{pax_tx / 1e3:.2f} kPa'],
            ['psi_lcfs tgt', f'{psi_lcfs_seed:.4f} Wb/rad', f'{psi_lcfs_tx:.4f} Wb/rad'],
            ['psi_lcfs TM', _psi_lcfs_tm_achieved_str, ''],
        ]
        diag_rows = [
            ['Parameter', 'Init EQDSK', 'TORAX'],
            ['q95', f'{q95_seed_eq:.3f}', f'{q95_tx:.3f}'],
            ['q0', f'{q0_seed_eq:.3f}', f'{q0_tx:.3f}'],
            ['v_loop', '\u2014', f'{vloop_tx:.3f} V'],
            ['beta_pol (TX)', '\u2014', f'{beta_pol_tx:.4f}'],
            ['beta_pol (TM)', '\u2014', _beta_pol_tm_str],
            ['beta_N', '\u2014', f'{beta_n_tx:.4f}'],
            ['TM nl_tol', _tol_str, ''],
        ]

        for blank_ax in [ax_ffp_tm, ax_pp_tm, ax_q_tm]:
            blank_ax.axis('off')
            blank_ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=11, color='gray',
                          transform=blank_ax.transAxes, fontweight='bold')

        render_table(ax_tbl1, input_rows, 'Scalar Inputs (Init EQDSK vs TORAX)')
        render_table(ax_tbl2, diag_rows, 'TORAX Diagnostics & Failure Info')
        ax_tbl2.text(
            0.5, -0.04,
            f'Fail: {_fail_reason_str}',
            transform=ax_tbl2.transAxes, fontsize=8, ha='center', va='top',
            color='darkred', fontweight='bold', wrap=True,
        )

        if _gs_ok:
            plt.suptitle(
                f'TM Diagnostic \u2014 loop {tt._current_loop}, t-idx {i}/{len(tt._tm_times) - 1}, t = {t:.2f} s'
                f'  |  TM GS OK \u2014 coupling failed (e.g. TORAX rejected EQDSK)',
                fontsize=13, color='darkorange', fontweight='bold',
            )
        else:
            plt.suptitle(
                f'TM Diagnostic \u2014 loop {tt._current_loop}, t-idx {i}/{len(tt._tm_times) - 1}, t = {t:.2f} s'
                f'  |  TokaMaker: FAILED',
                fontsize=13, color='darkred', fontweight='bold',
            )

    equil_snap = s.get('equil', {}).get(i)
    if _gs_ok and equil_snap is not None:
        tt._tm.plot_machine(
            fig, ax_mini_eq, equilibrium=equil_snap, coil_colormap='seismic', coil_symmap=False,
            coil_scale=1.E-6, coil_clabel=None,
        )
        tt._tm.plot_constraints(fig, ax_mini_eq, equilibrium=equil_snap)
        tt._tm.plot_psi(fig, ax_mini_eq, equilibrium=equil_snap, xpoint_color='r', vacuum_nlevels=3)
        x_pt = getattr(tt, '_x_point_targets', None)
        if x_pt is not None and _x_points_active(tt, i, t=t):
            ax_mini_eq.plot(x_pt[:, 0], x_pt[:, 1], 'x', color='purple', markersize=4, markeredgewidth=1)
        sp = s.get('strike_pts', {}).get(i)
        if sp is not None and len(sp) > 0:
            ax_mini_eq.plot(sp[:, 0], sp[:, 1], 'g^', markersize=4, markeredgewidth=1)
        ax_mini_eq.set_aspect('equal')
        ax_mini_eq.set_title('TM equilibrium', fontsize=9)
        ax_mini_eq.tick_params(labelsize=6)
    elif not _gs_ok:
        # TM failed to converge — current psi state may still be instructive
        _psi_plotted = False
        try:
            tt._tm.plot_machine(
                fig, ax_mini_eq, coil_colormap='seismic', coil_symmap=False,
                coil_scale=1.E-6, coil_clabel=None,
            )
            tt._tm.plot_psi(fig, ax_mini_eq, xpoint_color='r', vacuum_nlevels=3)
            ax_mini_eq.set_aspect('equal')
            ax_mini_eq.tick_params(labelsize=6)
            ax_mini_eq.text(
                0.02, 0.98, 'unconverged psi', transform=ax_mini_eq.transAxes, fontsize=7,
                ha='left', va='top', color='darkred', fontweight='bold',
            )
            _psi_plotted = True
        except Exception:
            pass
        if not _psi_plotted:
            ax_mini_eq.axis('off')
            ax_mini_eq.text(
                0.5, 0.5, 'TokaMaker did not converge\n(no psi available)',
                transform=ax_mini_eq.transAxes, fontsize=8,
                ha='center', va='center', color='darkred', fontweight='bold',
            )
        ax_mini_eq.set_title('TM psi (unconverged)', fontsize=9)
    else:
        ax_mini_eq.axis('off')
        if not solve_succeeded:
            _mini_msg = 'TM converged; EQDSK rejected by TORAX'
            _mini_color = 'darkorange'
        else:
            _mini_msg = 'No equilibrium snapshot'
            _mini_color = 'gray'
        ax_mini_eq.text(
            0.5, 0.5, _mini_msg, transform=ax_mini_eq.transAxes, fontsize=8,
            ha='center', va='center', color=_mini_color, fontweight='bold',
        )
        ax_mini_eq.set_title('TM equilibrium', fontsize=9)

    _save_or_display(fig, save_path, display)


# ── TM loop summary plot ──────────────────────────────────────────────────────

def tm_loop_summary_plot(tt, loop_level_log, save_path=None, display=True):
    r'''! Summary figure showing per-timestep GS solve outcomes.'''
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


# ── Profile evolution plot ────────────────────────────────────────────────────

def plot_profile_evolution(tt, save_path=None, display=True, one_plot=False):
    r'''! Plot profile evolution by pulse phase by default, or as one figure when one_plot=True.'''
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

        fig, axes = plt.subplots(2, 5, figsize=(25, 10))
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

        ax = axes[0, 3]
        ax.set_title(r'$\eta$')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$\eta$ [Ohm m]')
        plotted_eta = False
        for i_t in indices:
            if i_t in s.get('eta_prof', {}):
                color = cmap(norm(times[i_t]))
                ax.plot(s['eta_prof'][i_t]['x'], s['eta_prof'][i_t]['y'], color=color, linewidth=1.5, alpha=0.8)
                plotted_eta = True
        ax.set_xlim([0, 1])
        if plotted_eta:
            ax.set_yscale('log')

        ax = axes[1, 3]
        ax.set_title(r'$\psi$')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$\psi$ [Wb/rad]')
        for i_t in indices:
            if i_t in s.get('psi_tx', {}):
                color = cmap(norm(times[i_t]))
                ax.plot(s['psi_tx'][i_t]['x'], s['psi_tx'][i_t]['y'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_xlim([0, 1])

        ax = axes[0, 4]
        ax.set_title(r'$j_{\mathrm{tot}}$')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$j_{\mathrm{tot}}$ [MA/m$^2$]')
        for i_t in indices:
            if i_t in s.get('j_tot', {}):
                color = cmap(norm(times[i_t]))
                ax.plot(
                    s['j_tot'][i_t]['x'],
                    s['j_tot'][i_t]['y'] / 1e6,
                    color=color,
                    linewidth=1.5,
                    alpha=0.8,
                )
        ax.set_xlim([0, 1])

        ax = axes[1, 4]
        ax.set_title(r'$j_{\mathrm{NI}}$')
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$j_{\mathrm{NI}}$ [MA/m$^2$]')
        for i_t in indices:
            if i_t in s.get('j_ni', {}):
                color = cmap(norm(times[i_t]))
                ax.plot(
                    s['j_ni'][i_t]['x'],
                    s['j_ni'][i_t]['y'] / 1e6,
                    color=color,
                    linewidth=1.5,
                    alpha=0.8,
                )
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


# ── TORAX relax radial profiles (output_mode debug) ─────────────────────────

def plot_tx_relax_profiles(
    tt,
    data_tree,
    *,
    stage,
    t_final_relax,
    user_ref_curves=None,
    save_path=None,
    display=False,
):
    r'''! Plot psi, n_e, T_e, T_i vs rho_norm after a short TORAX relax.
    
        Shows main-run history (tt._relax_mainrun_profile_history), **user** kinetic
        profiles from flattened profile_conditions, **init EQDSK** psi at t_initial
        when psi is not prescribed (initial_psi_mode='geometry'), and the profile at
        t_final_relax. Used when TokaMaker_TORAX.fly(..., output_mode='debug').
        
    '''
    t_init = tt._t_init
    _sel_kw = dict(time=t_final_relax, method='nearest')
    _sel_init = dict(time=t_init, method='nearest')
    fig, axs = plt.subplots(2, 2, figsize=(8.5, 6.5))
    _hist = getattr(tt, '_relax_mainrun_profile_history', None) or []
    _uref = user_ref_curves or {}
    if stage == 'initial':
        _relax_loop = 1
    else:
        _relax_loop = tt._current_loop
    _stage_lbl = f'Loop {_relax_loop} relax'
    fig.suptitle(
        f'TORAX {_stage_lbl}: history, user inputs, init EQDSK psi, after relax',
        fontsize=11,
    )
    _cmap_hist = plt.get_cmap('tab10')
    _user_color = 'tab:green'
    _eqdsk_psi_color = 'tab:purple'

    def _one_panel(ax, var_name, ylabel):
        for _entry in _hist:
            lp = _entry['loop']
            profs = _entry['profiles']
            if var_name not in profs:
                continue
            tpr = profs[var_name]
            r_h = np.asarray(tpr[1], dtype=float)
            y_h = np.asarray(tpr[2][0], dtype=float)
            _c = _cmap_hist(((lp - 1) % 10) / 9.0)
            ax.plot(
                r_h,
                y_h,
                lw=1.1,
                ls=':',
                color=_c,
                alpha=0.95,
                label=f'loop {lp} main @ t_init',
            )

        if var_name in _uref:
            r_u, y_u = _uref[var_name]
            ax.plot(
                r_u,
                y_u,
                lw=1.25,
                ls='-.',
                color=_user_color,
                alpha=0.95,
                label=f'user {var_name} (config)',
            )
        elif var_name in ('n_e', 'T_e', 'T_i'):
            # If profile_conditions could not be parsed, show TORAX’s actual initial slice.
            try:
                xr_ic = getattr(data_tree.profiles, var_name).sel(**_sel_init)
                r_ic = xr_ic.coords['rho_norm'].to_numpy()
                y_ic = np.asarray(xr_ic.to_numpy())
                ax.plot(
                    r_ic,
                    y_ic,
                    lw=1.25,
                    ls='-.',
                    color='tab:cyan',
                    alpha=0.95,
                    label=f'{var_name} at t_init (TORAX IC)',
                )
            except Exception:
                pass
        elif var_name == 'psi':
            try:
                xr_geo = getattr(data_tree.profiles, 'psi').sel(**_sel_init)
                r_g = xr_geo.coords['rho_norm'].to_numpy()
                y_g = np.asarray(xr_geo.to_numpy())
                ax.plot(
                    r_g,
                    y_g,
                    lw=1.25,
                    ls='-.',
                    color=_eqdsk_psi_color,
                    alpha=0.95,
                    label=r'init EQDSK $\psi$ (t_init)',
                )
            except Exception:
                pass

        xr_post = getattr(data_tree.profiles, var_name).sel(**_sel_kw)
        r_post = xr_post.coords['rho_norm'].to_numpy()
        y_post = np.asarray(xr_post.to_numpy())
        ax.plot(
            r_post,
            y_post,
            lw=2.0,
            color='tab:blue',
            zorder=5,
            label='after this relax',
        )

        ax.set_xlabel(r'$\rho_{\mathrm{norm}}$')
        ax.set_ylabel(ylabel)
        ax.set_title(var_name)
        ax.set_xlim(0.0, 1.0)
        ax.legend(loc='best', fontsize=6, framealpha=0.92)

    _one_panel(axs[0, 0], 'psi', r'$\psi$ [Wb]')
    _one_panel(axs[0, 1], 'n_e', r'$n_e$ [m$^{-3}$]')
    _one_panel(axs[1, 0], 'T_e', r'$T_e$ [keV]')
    _one_panel(axs[1, 1], 'T_i', r'$T_i$ [keV]')
    fig.tight_layout()
    _save_or_display(fig, save_path, display)


# ── Scalar time-series plot ───────────────────────────────────────────────────

def plot_scalars(tt, save_path=None, display=True):
    r'''! Plot 4x3 grid of time-series scalars plus a bottom row: power channels, sources, P_LH.'''
    s = tt._state
    times = tt._tm_times

    _row_h = 3.0
    fig = plt.figure(figsize=(16, 5 * _row_h + 0.5 * _row_h))
    gs = fig.add_gridspec(5, 3, height_ratios=[1, 1, 1, 1, 1], hspace=0.35, wspace=0.3)
    axes = np.empty((4, 3), dtype=object)
    for i in range(4):
        for j in range(3):
            axes[i, j] = fig.add_subplot(gs[i, j])
    ax_power = fig.add_subplot(gs[4, 0])
    ax_sources = fig.add_subplot(gs[4, 1])
    ax_plh = fig.add_subplot(gs[4, 2])

    # (0,0): Ip
    ax = axes[0, 0]
    ax.set_title('Ip [A]')
    ax.plot(times, s['Ip_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='Ip TM')
    t_Ip, y_Ip = _tx_scalar(tt, 'Ip')
    if t_Ip is not None:
        ax.plot(t_Ip, y_Ip, color=COLOR_TX, ls='-', lw=1, label='Ip TX')
    t_Ini, y_Ini = _tx_scalar(tt, 'I_non_inductive')
    if t_Ini is not None:
        ax.plot(t_Ini, y_Ini, color=COLOR_TX, ls='--', lw=1, label='Ip NI TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Ip [A]')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # (0,1): psi_lcfs and psi_axis (left) + flux consumption (right, shares lcfs reference)
    ax = axes[0, 1]
    ax.set_title(r'$\psi_{lcfs}$ & $\psi_{axis}$; flux consumption')
    t_psi_lcfs, y_psi_lcfs = _tx_profile_at_rho(tt, 'psi', 1.0, scale=-1.0/(2.0*np.pi)) # TODO: just pull from state
    t_psi_axis, y_psi_axis = _tx_profile_at_rho(tt, 'psi', 0.0, scale=-1.0/(2.0*np.pi))
    ax.plot(times, s['psi_lcfs_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label=r'$\psi_{lcfs}$ TM')
    ax.plot(times, s['psi_axis_tm'], color=COLOR_TM, ls='--', marker='o', ms=MK_SZ, lw=1, label=r'$\psi_{axis}$ TM')
    if t_psi_lcfs is not None:
        ax.plot(t_psi_lcfs, y_psi_lcfs, color=COLOR_TX, ls='-', lw=1, label=r'$\psi_{lcfs}$ TX')
    if t_psi_axis is not None:
        ax.plot(t_psi_axis, y_psi_axis, color=COLOR_TX, ls='--', lw=1, label=r'$\psi_{axis}$ TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'$\psi$ [Wb/rad]')
    ax.grid(True, alpha=0.3)
    ax2_psi = ax.twinx()
    psi_lcfs_tm_arr = np.array(s['psi_lcfs_tm'])
    ax2_psi.plot(
        times,
        -(psi_lcfs_tm_arr - psi_lcfs_tm_arr[0]) * 2 * np.pi,
        color='darkorange',
        ls='-.',
        marker='o',
        ms=MK_SZ,
        lw=1,
        label='Flux TM',
    )
    if t_psi_lcfs is not None:
        flux_tx = -(y_psi_lcfs - y_psi_lcfs[0]) * 2 * np.pi
        ax2_psi.plot(t_psi_lcfs, flux_tx, color='seagreen', ls='-.', lw=1, label='Flux TX')
    ax2_psi.set_ylabel('Flux consumption [Wb]')
    ax2_psi.tick_params(axis='y')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2_psi.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc='upper left', ncol=2)

    # (0,2): V_loop
    ax = axes[0, 2]
    ax.set_title('V_loop (TM vs TX) [V]')
    ax.plot(times, s['vloop_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='TokaMaker')
    t_vl, y_vl = _tx_scalar(tt, 'v_loop_lcfs')
    if t_vl is not None:
        ax.plot(t_vl, y_vl, color=COLOR_TX, ls='-', lw=1, label='TORAX')
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
    ne_edge_y = [s['n_e'][ii]['y'][-1] if ii in s.get('n_e', {}) else np.nan for ii in range(len(times))]
    ax.plot(times, ne_edge_y, color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='n_e edge')
    t_ne_avg, y_ne_avg = _tx_scalar(tt, 'n_e_line_avg')
    t_ne_core, y_ne_core = _tx_profile_at_rho(tt, 'n_e', 0.0)
    if t_ne_avg is not None:
        ax.plot(t_ne_avg, y_ne_avg, color=COLOR_TX, ls='-', lw=1, label='n_e line avg')
    if t_ne_core is not None:
        ax.plot(t_ne_core, y_ne_core, color=COLOR_TX, ls='--', lw=1, label='n_e core')
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
    te_edge_y = [s['T_e'][ii]['y'][-1] if ii in s.get('T_e', {}) else np.nan for ii in range(len(times))]
    ax.plot(times, te_edge_y, color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='T_e edge')
    t_Te_vol_avg, y_Te_vol_avg = _tx_scalar(tt, 'T_e_volume_avg')
    t_Ti_vol_avg, y_Ti_vol_avg = _tx_scalar(tt, 'T_i_volume_avg')
    t_Te, y_Te = _tx_profile_at_rho(tt, 'T_e', 0.0)
    t_Ti, y_Ti = _tx_profile_at_rho(tt, 'T_i', 0.0)
    if t_Te_vol_avg is not None:
        ax.plot(t_Te_vol_avg, y_Te_vol_avg, color='darkorange', ls='-.', lw=1, label=r'$\langle T_e \rangle$')
    if t_Ti_vol_avg is not None:
        ax.plot(t_Ti_vol_avg, y_Ti_vol_avg, color='purple', ls='-.', lw=1, label=r'$\langle T_i \rangle$')
    if t_Te is not None:
        ax.plot(t_Te, y_Te, color=COLOR_TX, ls='-', lw=1, label='T_e core')
    if t_Ti is not None:
        ax.plot(t_Ti, y_Ti, color=COLOR_TX, ls='--', lw=1, label='T_i core')
    ax.set_xlabel('Time [s]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (2,0): scalar resistance R = eta * plasma cross-sectional area + TORAX v_loop
    ax = axes[2, 0]
    ax.set_title(r'Scalar resistance $R$ & $V_{loop}$ (TX)')
    area_tm = np.pi * np.asarray(s['a'], dtype=float) ** 2 * np.asarray(s['kappa'], dtype=float)
    eta_tm_scalar = np.full(len(times), np.nan, dtype=float)
    for ii in range(len(times)):
        eta_entry = s.get('eta_prof', {}).get(ii)
        if eta_entry is None:
            continue
        x_eta = np.asarray(eta_entry.get('x', []), dtype=float)
        y_eta = np.asarray(eta_entry.get('y', []), dtype=float)
        if y_eta.size == 0:
            continue
        if x_eta.size == y_eta.size and x_eta.size > 0:
            eta_tm_scalar[ii] = y_eta[np.argmin(np.abs(x_eta))]
        else:
            eta_tm_scalar[ii] = y_eta[0]
    r_tm_scalar = eta_tm_scalar * area_tm
    if np.any(np.isfinite(r_tm_scalar)):
        ax.plot(
            times, r_tm_scalar, color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1,
            label=r'$R$ TM ($\psi_N \approx 0$)',
        )
    tx_eta_names = ('eta', 'eta_parallel')
    for tx_eta_name in tx_eta_names:
        t_eta_tx, y_eta_tx = _tx_scalar(tt, tx_eta_name)
        if t_eta_tx is None:
            continue
        area_tx = np.interp(t_eta_tx, times, area_tm)
        r_tx = np.asarray(y_eta_tx, dtype=float) * area_tx
        ax.plot(
            t_eta_tx, r_tx, color=COLOR_TX, ls='-', lw=1.1, label=r'$R$ TX',
        )
        break
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'$R$ [$\Omega$]')
    r_positive = np.isfinite(r_tm_scalar) & (r_tm_scalar > 0)
    if np.any(r_positive):
        ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax2_eta = ax.twinx()
    t_vl_tx, y_vl_tx = _tx_scalar(tt, 'v_loop_lcfs')
    if t_vl_tx is not None:
        ax2_eta.plot(t_vl_tx, y_vl_tx, color='tab:green', ls='--', lw=1, label=r'$V_{loop}$ TX')
    ax2_eta.set_ylabel(r'$V_{loop}$ [V]')
    ax2_eta.tick_params(axis='y')
    h1_eta, l1_eta = ax.get_legend_handles_labels()
    h2_eta, l2_eta = ax2_eta.get_legend_handles_labels()
    if h1_eta or h2_eta:
        ax.legend(h1_eta + h2_eta, l1_eta + l2_eta, fontsize=8, loc='upper left')

    # (2,1): beta_N
    ax = axes[2, 1]
    ax.set_title('beta_N')
    ax.plot(times, s['beta_N_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='beta_N TM')
    t_bN, y_bN = _tx_scalar(tt, 'beta_N')
    if t_bN is not None:
        ax.plot(t_bN, y_bN, color=COLOR_TX, ls='-', lw=1, label='beta_N TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('beta_N')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (2,2): l_i
    ax = axes[2, 2]
    ax.set_title('l_i (li3)')
    ax.plot(times, s['l_i_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='l_i TM')
    t_li, y_li = _tx_scalar(tt, 'li3')
    if t_li is not None:
        ax.plot(t_li, y_li, color=COLOR_TX, ls='-', lw=1, label='l_i TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('l_i')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (3,0): q95 and q0
    ax = axes[3, 0]
    ax.set_title('Safety Factor q')
    t_q95, y_q95 = _tx_scalar(tt, 'q95')
    t_q0, y_q0 = _tx_profile_at_rho(tt, 'q', 0.0, rho_coord='rho_face_norm')
    ax.axhline(1.0, color='k', ls='-', lw=1) # q=1 horizontal line
    ax.plot(times, s['q95_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='q95 TM')
    ax.plot(times, s['q0_tm'], color=COLOR_TM, ls='--', marker='o', ms=MK_SZ, lw=1, label='q0 TM')
    if t_q95 is not None:
        ax.plot(t_q95, y_q95, color=COLOR_TX, ls='-', lw=1, label='q95 TX')
    if t_q0 is not None:
        ax.plot(t_q0, y_q0, color=COLOR_TX, ls='--', lw=1, label='q0 TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Safety Factor')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (3,1): pax
    ax = axes[3, 1]
    ax.set_title('pax [Pa]')
    ax.plot(times, s['pax_tm'], color=COLOR_TM, ls='-', marker='o', ms=MK_SZ, lw=1, label='pax TM')
    t_pax, y_pax = _tx_profile_at_rho(tt, 'pressure_thermal_total', 0.0)
    if t_pax is not None:
        ax.plot(t_pax, y_pax, color=COLOR_TX, ls='-', lw=1, label='pax TX')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('pax [Pa]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (3,2): tau_E (TORAX scalar) and H98 = H98(y,2) vs ITER98y2 (TORAX post_processing)
    ax = axes[3, 2]
    ax.set_title(r'$\tau_E$ & $H_{98(y,2)}$ (TX)')
    t_tau, y_tau = _tx_scalar(tt, 'tau_E')
    t_h98, y_h98 = _tx_scalar(tt, 'H98')
    if t_tau is not None:
        ax.plot(t_tau, y_tau, color=COLOR_TX, ls='-', lw=1, label=r'$\tau_E$')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'$\tau_E$ [s]')
    ax.grid(True, alpha=0.3)
    ax2_te = ax.twinx()
    if t_h98 is not None:
        ax2_te.plot(t_h98, y_h98, color='crimson', ls='--', lw=1, label=r'$H_{98(y,2)}$')
    ax2_te.set_ylabel(r'$H_{98(y,2)}$')
    ax2_te.tick_params(axis='y')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2_te.get_legend_handles_labels()
    if h1 or h2:
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left')

    # Bottom row (left): power channels; P_LH overlaid here and shown on log axis at right.
    ax = ax_power
    ax.set_title('Power channels [W]')
    power_rows = [('P_ohmic_e', 'P_ohmic_e', 'r-'), ('P_radiation_e', 'P_radiation_e', 'm--'),
                  ('P_SOL_total', 'P_SOL_total', 'c--'), ('P_alpha_total', 'P_alpha_total', 'g-.'),
                  ('P_aux_total', 'P_aux_total', 'y-.')]
    y_non_lh = []
    for tx_name, label, fmt in power_rows:
        scale = -1.0 if tx_name == 'P_radiation_e' else 1.0
        tx_t, tx_y = _tx_scalar(tt, tx_name, scale=scale)
        if tx_t is not None:
            ax.plot(tx_t, tx_y, fmt, linewidth=1, label=label)
            y_non_lh.append(tx_y)
    if y_non_lh:
        y_all = np.concatenate([np.asarray(a, dtype=float).ravel() for a in y_non_lh])
        y_all = y_all[np.isfinite(y_all)]
        if y_all.size:
            lo, hi = float(np.min(y_all)), float(np.max(y_all))
            span = hi - lo
            pad = 0.05 * span if span > 0 else 0.05 * max(abs(hi), abs(lo), 1.0)
            ax.set_ylim(lo - pad, hi + pad)
    t_plh_pc, y_plh_pc = _tx_scalar(tt, 'P_LH')
    t_plh_pc_d, y_plh_pc_d = _tx_scalar(tt, 'P_LH_delabie')
    if t_plh_pc is not None:
        ax.plot(t_plh_pc, y_plh_pc, 'k-', linewidth=1, label='P_LH_martin', zorder=5)
    if t_plh_pc_d is not None:
        ax.plot(t_plh_pc_d, y_plh_pc_d, 'k--', linewidth=1, label='P_LH_delabie', zorder=4)
    ax.set_xlabel('Time [s]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Bottom row (right): L–H transition power threshold, showing the high- and
    # low-density branches of the Martin and Delabie scalings. The selected P_LH
    # rides whichever branch is active; the crossover is at n_e_min_P_LH.
    ax = ax_plh
    ax.set_title(r'$P_{\mathrm{LH}}$ [W]')
    t_plh, y_plh = _tx_scalar(tt, 'P_LH')
    if t_plh is not None:
        ax.plot(t_plh, y_plh, 'k-', linewidth=1.4, label=r'$P_{\mathrm{LH}}$ Martin', zorder=6)
        t_m_hd, y_m_hd = _tx_scalar(tt, 'P_LH_high_density')
        t_m_ld, y_m_ld = _tx_scalar(tt, 'P_LH_low_density')
        if t_m_hd is not None:
            ax.plot(t_m_hd, y_m_hd, color='tab:blue', ls='--', lw=1, label='Martin high-n')
        if t_m_ld is not None:
            ax.plot(t_m_ld, y_m_ld, color='tab:cyan', ls=':', lw=1, label='Martin low-n')
        t_d, y_d = _tx_scalar(tt, 'P_LH_delabie')
        t_d_hd, y_d_hd = _tx_scalar(tt, 'P_LH_delabie_high_density')
        t_d_ld, y_d_ld = _tx_scalar(tt, 'P_LH_delabie_low_density')
        if t_d is not None:
            ax.plot(t_d, y_d, color='tab:red', ls='-', lw=1.4, label='Delabie', zorder=5)
        if t_d_hd is not None:
            ax.plot(t_d_hd, y_d_hd, color='tab:orange', ls='--', lw=1, label='Delabie high-n')
        if t_d_ld is not None:
            ax.plot(t_d_ld, y_d_ld, color='tab:pink', ls=':', lw=1, label='Delabie low-n')
        ax.legend(fontsize=7, loc='upper left', ncol=2)
    else:
        ax.text(
            0.5, 0.5, r'$P_{\mathrm{LH}}$ not in TORAX output',
            transform=ax.transAxes, ha='center', va='center', fontsize=9, color='gray',
        )
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Power [W]')
    ax.grid(True, alpha=0.3)

    def _tx_scalar_if_nonzero(name, scale=1.0):
        t, y = _tx_scalar(tt, name, scale=scale)
        if t is None or y is None:
            return None, None
        ya = np.asarray(y, dtype=float)
        if not np.any(np.isfinite(ya) & (np.abs(ya) > 0)):
            return None, None
        return t, y

    ax_src = ax_sources
    ax_src.set_title('Sources')
    heat_cfg = [
        ('P_aux_total', 'Aux', 'tab:blue', '-'),
        ('P_ohmic_e', 'Ohmic', 'tab:red', '-'),
        ('P_alpha_total', 'Alpha', 'tab:green', '-'),
    ]
    plotted_heat = False
    for tx_name, label, clr, ls in heat_cfg:
        tx_t, tx_y = _tx_scalar_if_nonzero(tx_name)
        if tx_t is None:
            continue
        ax_src.plot(tx_t, tx_y, color=clr, ls=ls, lw=1.2, label=label)
        plotted_heat = True
    ax_src.set_ylabel('Heating power [W]')
    ax_src.set_xlabel('Time [s]')
    ax_src.grid(True, alpha=0.3)

    fuel_cfg = [
        ('S_gas_puff', 'Gas puff $S$', 'tab:brown', '-'),
        ('S_pellet', 'Pellet $S$', 'tab:purple', '--'),
    ]
    ax2_src = None
    for tx_name, label, clr, ls in fuel_cfg:
        tx_t, tx_y = _tx_scalar_if_nonzero(tx_name)
        if tx_t is None:
            continue
        if ax2_src is None:
            ax2_src = ax_src.twinx()
        ax2_src.plot(tx_t, tx_y, color=clr, ls=ls, lw=1.2, label=label)
    if ax2_src is not None:
        ax2_src.set_ylabel(r'Particle source $S$ [s$^{-1}$]')
        ax2_src.tick_params(axis='y')

    h1s, l1s = ax_src.get_legend_handles_labels()
    h2s, l2s = ax2_src.get_legend_handles_labels() if ax2_src is not None else ([], [])
    if h1s or h2s:
        ax_src.legend(h1s + h2s, l1s + l2s, fontsize=8, loc='upper left', ncol=2)
    elif not plotted_heat:
        ax_src.text(
            0.5, 0.5, 'No non-zero source scalars in TORAX output',
            transform=ax_src.transAxes, ha='center', va='center', fontsize=9, color='gray',
        )

    plt.suptitle('Scalars', fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    _save_or_display(fig, save_path, display)


# ── P_LH component plot ───────────────────────────────────────────────────────
def plot_PLH_components(tt, save_path=None, display=True):
    r'''! Decompose the L->H transition power into every ingredient, for the
    current loop's TORAX run. 3x3 grid (all from tt._data_tree at full TORAX
    resolution; B_T is ~constant and printed in the title):

      (0,0) line-avg n_e vs n_e_min (Ryter) — which density BRANCH is active.
      (0,1) Ip [MA]            (drives n_min ∝ Ip^0.34).
      (0,2) S = g0_face[-1] [m^2] (LCFS surface area; P_LH ∝ S^0.941).
      (1,0) a_minor, R_major [m] (set n_min; a also via R/a).
      (1,1) the multiplicative P_LH factors: B^0.803, n̄^0.717, S^0.941, (2/Meff).
      (1,2) P_LH branch decomposition (selected / high-n / low-n n^-2 / floor).
      (2,0) f_GW (line-avg)    — headroom to the density limit.
      (2,1) P_SOL vs P_LH      — the transition gate (H-mode where P_SOL>P_LH).
      (2,2) the governing equations (LaTeX).

    Isolates why the L->H transition can fire in one loop but not another even
    when geometry/Ip are nearly identical: the branch crossover is set by
    line-avg n_e vs n_min, and the low-density n^-2 branch is hypersensitive there.
    '''
    if getattr(tt, '_data_tree', None) is None:
        return

    # Scaling exponents (torax scaling_laws._P_LH_SCALING_PARAMS). Prefactors
    # (Martin 0.0488 / Delabie 0.0441) only appear in the LaTeX panel.
    # Martin:  B^0.803 n^0.717 (2/Meff)^1.0   S^0.941
    # Delabie: B^0.580 n^1.08  (2/Meff)^0.975 S^1.0   (x divertor factor D)
    aB, an, aM, aS = 0.803, 0.717, 1.0, 0.941          # Martin
    aB_d, an_d, aM_d, aS_d = 0.580, 1.08, 0.975, 1.0   # Delabie
    Meff = 2.0  # (2/Meff)^aM factor; main-ion mass [amu], ~2-2.5 for D-T.
    try:
        Meff = float(np.asarray(getattr(tt._data_tree.scalars, 'A_i').to_numpy()).mean())
    except Exception:
        pass

    def geom_face(name, scale=1.0):
        # LCFS (rho=1) value of a profile geometry field vs time.
        try:
            return _tx_profile_at_rho(tt, name, 1.0, scale=scale)
        except Exception:
            return None, None

    t_n,  y_n  = _tx_scalar(tt, 'n_e_line_avg', scale=1e-20)   # 1e20 m^-3
    t_nm, y_nm = _tx_scalar(tt, 'n_e_min_P_LH', scale=1e-20)
    t_Ip, y_Ip = _tx_scalar(tt, 'Ip', scale=1e-6)             # MA
    t_a,  y_a  = _tx_scalar(tt, 'a_minor')
    t_R,  y_R  = _tx_scalar(tt, 'R_major')
    t_B,  y_B  = _tx_scalar(tt, 'B_0')
    if t_B is None:
        t_B, y_B = geom_face('R_major_profile')  # fallback; B handled below
    t_S,  y_S  = geom_face('g0')                              # surface area [m^2]
    t_fg, y_fg = _tx_scalar(tt, 'fgw_n_e_line_avg')
    B0 = float(np.nanmean(y_B)) if (y_B is not None) else float(getattr(tt, '_B0', np.nan))

    fig, axes = plt.subplots(3, 3, figsize=(17, 13))

    # (0,0) density vs branch threshold
    ax = axes[0, 0]
    if t_n is not None:
        ax.plot(t_n, y_n, color=COLOR_TX, ls='-', marker='o', ms=MK_SZ, lw=1, label=r'$\bar n_e$ (line-avg)')
    if t_nm is not None:
        ax.plot(t_nm, y_nm, color='tab:red', ls='--', lw=1.2, label=r'$n_{e,\min}$ (Ryter)')
    if t_n is not None and t_nm is not None:
        ax.fill_between(t_n, y_n, y_nm, where=(np.asarray(y_n) < np.asarray(y_nm)),
                        color='tab:red', alpha=0.15, label=r'$\bar n_e<n_{e,\min}$ ($n^{-2}$ branch)')
    ax.set_title(r'density vs branch threshold'); ax.set_ylabel(r'$n_e$ [$10^{20}$ m$^{-3}$]'); ax.legend(fontsize=8)

    # (0,1) Ip
    ax = axes[0, 1]
    if t_Ip is not None:
        ax.plot(t_Ip, y_Ip, color=COLOR_TX, ls='-', marker='o', ms=MK_SZ, lw=1)
    ax.set_title(r'$I_p$  (sets $n_{e,\min}\propto I_p^{0.34}$)'); ax.set_ylabel('[MA]')

    # (0,2) S = surface area
    ax = axes[0, 2]
    if t_S is not None:
        ax.plot(t_S, y_S, color=COLOR_TX, ls='-', marker='o', ms=MK_SZ, lw=1)
    ax.set_title(r'$S$ = g0_face[-1]  ($P_{LH}\propto S^{0.941}$)'); ax.set_ylabel(r'[m$^2$]')

    # (1,0) a_minor, R_major
    ax = axes[1, 0]
    if t_a is not None:
        ax.plot(t_a, y_a, color='tab:blue', ls='-', marker='o', ms=MK_SZ, lw=1, label=r'$a$ (minor)')
    if t_R is not None:
        ax.plot(t_R, y_R, color='tab:purple', ls='-', marker='s', ms=MK_SZ, lw=1, label=r'$R$ (major)')
    ax.set_title(r'$a$, $R$  (enter $n_{e,\min}$)'); ax.set_ylabel('[m]'); ax.legend(fontsize=8)

    # (1,1) the multiplicative P_LH factors (Martin solid, Delabie dashed)
    ax = axes[1, 1]
    if t_n is not None:
        ax.plot(t_n, (np.asarray(y_n)) ** an,   'C0-',  lw=1.2, label=r'M $\bar n_{e,20}^{0.717}$')
        ax.plot(t_n, (np.asarray(y_n)) ** an_d, 'C0--', lw=1.0, label=r'D $\bar n_{e,20}^{1.08}$')
    if t_S is not None:
        ax.plot(t_S, (np.asarray(y_S)) ** aS,   'C1-',  lw=1.2, label=r'M $S^{0.941}$')
        ax.plot(t_S, (np.asarray(y_S)) ** aS_d, 'C1--', lw=1.0, label=r'D $S^{1.0}$')
    ax.axhline(B0 ** aB,   color='C2', ls='-',  lw=1, label=fr'M $B_T^{{0.803}}={B0**aB:.1f}$')
    ax.axhline(B0 ** aB_d, color='C2', ls='--', lw=1, label=fr'D $B_T^{{0.580}}={B0**aB_d:.1f}$')
    ax.axhline((2.0 / Meff) ** aM, color='C3', ls=':', lw=1, label=fr'$(2/M_{{eff}})\approx{(2.0/Meff)**aM:.2f}$')
    ax.set_title('mult. $P_{LH}$ factors (M solid / D dashed)'); ax.set_ylabel('factor'); ax.set_yscale('log'); ax.legend(fontsize=7, ncol=2)

    # (1,2) P_LH branch decomposition (log) — Martin AND Delabie
    ax = axes[1, 2]
    for name, color, ls, lw, lab in [
        ('P_LH',                       'k',          '-',  1.6, r'$P_{LH}$ Martin (sel.)'),
        ('P_LH_high_density',          'C0',         '--', 1.0, 'M high-n'),
        ('P_LH_low_density',           'C1',         ':',  1.0, r'M low-n ($n^{-2}$)'),
        ('P_LH_min',                   'C2',         '-.', 1.0, 'M floor'),
        ('P_LH_delabie',               'tab:red',    '-',  1.6, r'$P_{LH}$ Delabie (sel.)'),
        ('P_LH_delabie_high_density',  'tab:orange', '--', 1.0, 'D high-n'),
        ('P_LH_delabie_low_density',   'tab:pink',   ':',  1.0, r'D low-n ($n^{-2}$)'),
        ('P_LH_delabie_min',           'tab:brown',  '-.', 1.0, 'D floor'),
    ]:
        t, y = _tx_scalar(tt, name, scale=1e-6)
        if t is not None:
            ax.plot(t, y, color=color, ls=ls, lw=lw, label=lab)
    ax.set_title('P_LH branch decomposition (Martin + Delabie)'); ax.set_ylabel('[MW]'); ax.set_yscale('log'); ax.legend(fontsize=7, ncol=2)

    # (2,0) Greenwald fraction
    ax = axes[2, 0]
    if t_fg is not None:
        ax.plot(t_fg, y_fg, color=COLOR_TX, ls='-', marker='o', ms=MK_SZ, lw=1)
    ax.axhline(1.0, color='tab:red', ls='--', lw=1, label='Greenwald limit')
    ax.set_title(r'$f_{GW}$ (line-avg)'); ax.set_ylabel(r'$f_{GW}$'); ax.legend(fontsize=8)

    # (2,1) P_SOL vs P_LH gate (Delabie is what the formation model uses here)
    ax = axes[2, 1]
    t_lh, y_lh = _tx_scalar(tt, 'P_LH', scale=1e-6)
    t_lhd, y_lhd = _tx_scalar(tt, 'P_LH_delabie', scale=1e-6)
    t_sol, y_sol = _tx_scalar(tt, 'P_SOL_total', scale=1e-6)
    if t_lh is not None:
        ax.plot(t_lh, y_lh, 'k--', lw=1, label=r'$P_{LH}$ Martin')
    if t_lhd is not None:
        ax.plot(t_lhd, y_lhd, color='tab:red', ls='-', marker='o', ms=MK_SZ, lw=1.4, label=r'$P_{LH}$ Delabie (used)')
    if t_sol is not None:
        ax.plot(t_sol, y_sol, color='tab:green', ls='-', marker='o', ms=MK_SZ, lw=1, label=r'$P_{SOL}$')
    # H-mode shading vs the Delabie threshold (the one the formation model gates on)
    t_gate, y_gate = (t_lhd, y_lhd) if t_lhd is not None else (t_lh, y_lh)
    if t_gate is not None and t_sol is not None:
        ax.fill_between(t_sol, y_sol, y_gate, where=(np.asarray(y_sol) > np.asarray(y_gate)),
                        color='tab:green', alpha=0.2, label=r'$P_{SOL}>P_{LH}$ (H-mode)')
    ax.set_title(r'transition gate: $P_{SOL}$ vs $P_{LH}$'); ax.set_ylabel('[MW]'); ax.legend(fontsize=8)

    # (2,2) the governing equations (LaTeX)
    ax = axes[2, 2]; ax.axis('off')
    # NOTE: matplotlib mathtext (not a full LaTeX engine) — no \mathbf, \begin{cases},
    # \! or \qquad. Use plain-text headers and split the piecewise into two lines.
    eqn_lines = [
        ('High-density branch:', False),
        (r'$P_{LH}^{hi}[MW] = C\, B_T^{a_B}\, \bar{n}_{e,20}^{a_n}\, (2/M_{eff})^{a_M}\, D\, S^{a_S}$', True),
        ('  Martin:  C=0.0488, a_B=0.803, a_n=0.717, a_M=1.0, a_S=0.941, D=1', False),
        ('  Delabie: C=0.0441, a_B=0.580, a_n=1.08, a_M=0.975, a_S=1.0', False),
        ('           D=1 (HT) / 1.93 (VT)', False),
        ('', False),
        ('Ryter density minimum:', False),
        (r'$n_{e,min} = 0.7\,(I_p/MA)^{0.34}\,a^{-0.95}\,B_0^{0.62}\,(R/a)^{0.4}\times 10^{19}$', True),
        ('', False),
        ('Branch selection:', False),
        (r'$\bar{n}_e > n_{e,min}:\ \ P_{LH} = P_{LH}^{hi}(\bar{n}_e)$', True),
        (r'$\bar{n}_e < n_{e,min}:\ \ P_{LH} = P_{LH}^{hi}(n_{e,min})\,(n_{e,min}/\bar{n}_e)^{2}$  ($n^{-2}$)', True),
        ('', False),
        (r'Transition: H-mode when  $P_{SOL} > P_{LH}\cdot$prefactor', True),
    ]
    y = 0.99
    for txt, is_math in eqn_lines:
        ax.text(0.0, y, txt, transform=ax.transAxes, va='top', ha='left',
                fontsize=10 if is_math else 9.5,
                fontweight=('normal' if is_math else 'bold') if txt.endswith(':') else 'normal')
        y -= 0.072 if txt else 0.04

    for ax in axes.flat:
        if ax is not axes[2, 2]:
            ax.grid(True, alpha=0.3); ax.set_xlabel('Time [s]')

    _loop = getattr(tt, '_current_loop', '?')
    plt.suptitle(fr'P_LH components — loop {_loop}   ($B_T \approx {B0:.2f}$ T,  $M_{{eff}}\approx{Meff:.2f}$)', fontsize=14)
    plt.tight_layout()
    _save_or_display(fig, save_path, display)


# ── Coil current plot ─────────────────────────────────────────────────────────

def _coil_net_turns(tt, cname):
    r'''! Return net_turns for a coil set name, defaulting to 1.'''
    return tt._tm.coil_sets.get(cname, {}).get('net_turns', 1.0)


def _coil_aturn_clim(tt):
    r'''! Return (min, max) colormap limits in A-turns from _coil_bounds (A/turn * turns).'''
    coil_bounds = getattr(tt, '_coil_bounds', {})
    if not coil_bounds:
        return -1.0, 1.0
    vals = []
    for cname, (lo, hi) in coil_bounds.items():
        n = _coil_net_turns(tt, cname)
        vals.extend([lo * n, hi * n])
    return min(vals), max(vals)


def plot_coils(tt, save_path=None, display=True):
    r'''! Plot coil current traces in MA-turns with limit bands.'''
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


# ── LCFS evolution plot ───────────────────────────────────────────────────────

_LCFS_EVO_PSI_TRY = (1.0, 0.99, 0.98, 0.97)


def _trace_lcfs_for_evolution_plot(equil):
    r'''! Try successive normalized psi surfaces until tracing returns a non-empty contour.'''
    for psi in _LCFS_EVO_PSI_TRY:
        try:
            lcfs = equil.trace_surf(psi)
        except Exception:
            lcfs = None
        if lcfs is not None and len(lcfs) > 0:
            return lcfs
    return None


def plot_lcfs_evolution(tt, save_path=None, display=True, one_plot=False):
    r'''! Plot time evolution of the last closed flux surface for each phase.
    
        Produces phase-split figures by default (rampup, flattop, rampdown),
        or one combined figure when one_plot=True.
        
    '''
    s = tt._state
    times = np.array(tt._tm_times)
    lcfs_tm_data = s.get('lcfs_geo_tm', {})
    equil_data = s.get('equil', {})
    if not lcfs_tm_data and not equil_data:
        return

    ft = getattr(tt, '_flattop', np.zeros(len(times), dtype=bool)).astype(bool)

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
            fig.suptitle(f'LCFS Evolution \u2014 {phase_name} (loop {tt._current_loop})', fontsize=TITLE_FS)

        for lc in lim_contours:
            lc = np.asarray(lc)
            if mirror_mode:
                ax.plot(lc[:, 1], lc[:, 0], color='k', linewidth=1.5, zorder=5)
            else:
                ax.plot(lc[:, 0], lc[:, 1], color='k', linewidth=1.5, zorder=5)

        for i in indices:
            lcfs = lcfs_tm_data.get(i)
            if lcfs is None:
                equil = equil_data.get(i)
                if equil is None:
                    continue
                lcfs = _trace_lcfs_for_evolution_plot(equil)
                if lcfs is None or len(lcfs) == 0:
                    if getattr(tt, '_output_mode', None) == 'debug' and hasattr(tt, '_print'):
                        psi_str = ', '.join(str(p) for p in _LCFS_EVO_PSI_TRY)
                        tt._print(
                            f'  LCFS evolution: no contour traced at time index {i} '
                            f'(t={times[i]:.4g} s, tried psi_N={psi_str})'
                        )
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


# ── Movie generation ──────────────────────────────────────────────────────────

def make_movie(tt, save_path=None, display=True, speed_factor=1.0, loop=None, notebook_mode=None):
    r'''! Create pulse movie from simulation data.
    
        Renders equilibrium plots from stored equilibrium snapshots, generates
        composite frames in a temp directory, encodes to MP4, then cleans up.
    
        Parameters
        ----------
        tt : TokaMaker_TORAX
            Fully-populated TokaMaker_TORAX object (after fly()).
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
        
    '''
    if loop is None:
        loop = tt._current_loop

    tmp_dir = _make_temp_dir_viz()
    try:
        times = tt._tm_times
        n = len(times)

        psi_lcfs_tm = np.array(tt._state['psi_lcfs_tm'])
        psi_lcfs_tx = np.array(tt._state['psi_lcfs_tx'])
        flux_con_tm = -(psi_lcfs_tm - psi_lcfs_tm[0]) * 2.0 * np.pi
        flux_con_tx = -(psi_lcfs_tx - psi_lcfs_tx[0]) * 2.0 * np.pi

        equil_dir = os.path.join(tmp_dir, 'equil')
        os.makedirs(equil_dir, exist_ok=True)
        _render_equil_frames(tt, loop, equil_dir)

        for idx in range(n):
            fpath = os.path.join(tmp_dir, f'frame_{idx:04d}.png')
            _render_frame(tt, loop, idx, times[idx], times,
                          flux_con_tm, flux_con_tx, fpath, tt._run_name, equil_dir)
            plt.close('all')

        total_time = times[-1] - times[0] if len(times) > 1 else 1.0
        fps = max(1, speed_factor * n / total_time)

        if save_path is None:
            save_path = os.path.join(os.getcwd(), f'TokaMaker_TORAX_pulse_loop{loop:03d}.mp4')

        _encode_video_from_dir(tmp_dir, save_path, fps=fps)

        show_in_nb = notebook_mode if notebook_mode is not None else (display and _in_jupyter())
        if show_in_nb:
            _embed_video_jupyter(save_path)

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _render_equil_frames(tt, loop, equil_dir):
    r'''! Render equilibrium plots from stored equilibrium snapshots to PNG files.'''
    equil_data = tt._state.get('equil', {})

    for i in range(len(tt._tm_times)):
        out_path = os.path.join(equil_dir, f'equil_{loop:03d}.{i:03d}.png')
        equil = equil_data.get(i)
        if equil is None:
            fig, ax = plt.subplots(1, 1, figsize=(11, 12))
            ax.text(0.5, 0.5, 'TokaMaker failed\nto converge',
                    transform=ax.transAxes, fontsize=16,
                    ha='center', va='center', color='darkred', fontweight='bold')
            ax.axis('off')
            fig.savefig(out_path, dpi=MOVIE_DPI, bbox_inches='tight', pad_inches=0.0)
            plt.close(fig)
            continue

        fig, ax = plt.subplots(1, 1, figsize=(11, 12))
        min_bound, max_bound = _coil_aturn_clim(tt)
        cb = tt._tm.plot_machine(fig, ax, equilibrium=equil, coil_colormap='seismic', coil_symmap=False,
                                  coil_scale=1.E-6, coil_clabel=r'$I_C$ [MA-turns]')
        tt._tm.plot_constraints(fig, ax, equilibrium=equil)
        if cb is not None:
            cb.mappable.set_clim(min_bound * 1e-6, max_bound * 1e-6)
        tt._tm.plot_psi(
            fig, ax, equilibrium=equil, xpoint_color='r', vacuum_nlevels=4,
            plasma_nlevels=MOVIE_EQUIL_PSI_PLASMA_NLEVELS,
        )
        x_pt = getattr(tt, '_x_point_targets', None)
        if x_pt is not None and _x_points_active(tt, i, t=tt._tm_times[i]):
            ax.plot(x_pt[:, 0], x_pt[:, 1], 'x', color='purple', markersize=10, markeredgewidth=2, label='Saddle point targets')
        sp = tt._state.get('strike_pts', {}).get(i)
        if sp is not None and len(sp) > 0:
            ax.plot(sp[:, 0], sp[:, 1], 'g^', markersize=10, markeredgewidth=2, label='Strike point targets')
        ax.set_aspect('equal')
        handles, _labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc='upper right', fontsize=12)
        ax.tick_params(labelsize=11)
        if cb is not None:
            cb.ax.tick_params(labelsize=11)
        fig.savefig(out_path, dpi=MOVIE_DPI, bbox_inches='tight', pad_inches=0.0)
        plt.close(fig)


def _render_frame(tt, loop, idx, t_now, times, flux_con_tm, flux_con_tx, out_path, run_name, equil_dir):
    r'''! Create a single movie frame.'''
    fig = plt.figure(figsize=(MOVIE_FIG_W, MOVIE_FIG_H), dpi=MOVIE_DPI)
    gs = GridSpec(6, 3, figure=fig, width_ratios=[1.2, 1.0, 1.0],
                  wspace=0.28, hspace=0.18,
                  left=0.045, right=0.94, top=0.96, bottom=0.05)

    ax_text = fig.add_subplot(gs[0:1, 0])
    ax_equil = fig.add_subplot(gs[1:5, 0])
    _draw_info(ax_text, tt, loop, run_name)
    _draw_equil_from_dir(ax_equil, tt, loop, idx, t_now, equil_dir)

    sax = [fig.add_subplot(gs[j, 1]) for j in range(6)]
    _draw_scalars_movie(sax, tt, times, t_now, flux_con_tm, flux_con_tx)
    for ax in sax[:-1]:
        ax.tick_params(labelbottom=False)
    sax[-1].set_xlabel('Time [s]', fontsize=LABEL_FS)

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
    grid_type = getattr(tt, '_tx_grid_type', None)
    gr = getattr(tt, '_tx_grid', None)
    if grid_type == 'n_rho':
        n_rho = int(gr) if np.isscalar(gr) else len(gr)
        grid_label = f'n_rho = {n_rho}'
    elif grid_type == 'face_centers':
        grid_label = f'face_centers = {len(np.atleast_1d(gr))}'
    elif np.isscalar(gr):
        grid_label = f'n_rho = {int(gr)}'
    else:
        grid_label = f'face_centers = {len(np.atleast_1d(gr))}'
    lines = [
        f'Run:   {run_name}            loop:  {loop}',
        f'{grid_label}          dt:    {dt_val} s',
        f'time range:     [{tt._t_init}, {tt._t_final}] s          times: {len(tt._tm_times)}',
        f'LSF:   {tt._last_surface_factor}',
        f't_ave: {t_ave}    window: {getattr(tt, "_t_ave_window", 0):.2f} s',
        f'causal: {getattr(tt, "_t_ave_causal", True)}    '
        f'ignore_start: {getattr(tt, "_t_ave_ignore_start", 0):.2f} s',
    ]
    ax.text(0.05, 0.95, '\n'.join(lines), transform=ax.transAxes, fontsize=INFO_FS,
            va='top', fontfamily='monospace')


def _draw_equil_from_dir(ax, tt, loop, idx, t_now, equil_dir):
    r'''! Draw equilibrium from pre-rendered PNG in equil_dir.'''
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
    r'''! Draw the 6 scalar time-series panels for a movie frame.'''
    s = tt._state

    ax = axes[0]
    ax.plot(times, np.array(s['Ip_tm']) / 1e6, color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label='Ip TM')
    ax.plot(times, np.array(s['Ip_tx']) / 1e6, color=COLOR_TX, ls=LS_PRI, lw=LW, label='Ip TX')
    ax.plot(times, np.array(s['Ip_ni_tx']) / 1e6, color=COLOR_TX, ls=LS_SEC, lw=LW, label='Ip_ni TX')
    ax.set_ylabel('Ip [MA]', fontsize=LABEL_FS)
    ax.set_title('Scalars', fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_FS, loc='upper left')
    _style(ax)

    ax = axes[1]
    ax.plot(times, s['vloop_tx'], color=COLOR_TX, ls=LS_PRI, lw=LW, label='Vloop TX')
    ax.set_ylabel('V_loop [V]', fontsize=LABEL_FS)
    ax.legend(fontsize=LEGEND_FS, loc='upper left')
    _style(ax)
    ax2 = ax.twinx()
    if s.get('l_i_tm', None) is not None:
        ax2.plot(times, s['l_i_tm'], color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label='l_i TM')
    if s.get('l_i_tx', None) is not None:
        ax2.plot(times, s['l_i_tx'], color=COLOR_TX, ls=LS_SEC, lw=LW, label='l_i TX')
    ax2.set_ylabel('l_i', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    ax2.legend(fontsize=LEGEND_FS, loc='upper right')

    ax = axes[2]
    pkeys = [
        ('P_ohmic_e', 'Ohmic', COLORS_MULTI[0]),
        ('P_aux_total', 'Aux', COLORS_MULTI[1]),
        ('P_alpha_total', 'Fusion', COLORS_MULTI[2]),
        ('P_radiation_e', 'Radiation', COLORS_MULTI[3]),
        ('P_SOL_total', 'SOL', COLORS_MULTI[4]),
        ('P_LH', 'P_LH', COLORS_MULTI[5]),
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
    ax2.legend(fontsize=LEGEND_FS, loc='upper right')

    ax = axes[3]
    ax.plot(times, s['psi_axis_tm'], color=COLOR_TM, ls=LS_PRI, lw=LW, marker=MK_TM, ms=MK_SZ, label='\u03c8_axis TM')
    ax.plot(times, s['psi_axis_tx'], color=COLOR_TX, ls=LS_PRI, lw=LW, label='\u03c8_axis TX')
    ax.set_ylabel('\u03c8 [Wb/rad]', fontsize=LABEL_FS)
    ax.plot(times, s['psi_lcfs_tm'], color=COLOR_TM, ls=LS_SEC, lw=LW, marker=MK_TM, ms=MK_SZ, label='\u03c8_lcfs TM')
    ax.plot(times, s['psi_lcfs_tx'], color=COLOR_TX, ls=LS_SEC, lw=LW, label='\u03c8_lcfs TX')
    ax.tick_params(labelsize=TICK_FS)
    _style(ax)
    ax.legend(fontsize=LEGEND_FS, loc='upper left')
    ax2 = ax.twinx()
    ax2.plot(times, flux_con_tm, color='darkorange', ls='-.', lw=LW, marker=MK_TM, ms=MK_SZ, label='Flux TM')
    ax2.plot(times, flux_con_tx, color='seagreen', ls='-.', lw=LW, label='Flux TX')
    ax2.set_ylabel('Flux Consumed [Wb]', fontsize=LABEL_FS)
    ax2.tick_params(labelsize=TICK_FS)
    ax2.legend(fontsize=LEGEND_FS, loc='upper right')

    coil_data = tt._results.get('COIL', {})
    cs_coils = []
    other_coils = []
    for cname, cvals in sorted(coil_data.items()):
        cname_tokens = [tok for tok in cname.upper().replace('-', '_').split('_') if tok]
        if any(tok.startswith('CS') for tok in cname_tokens):
            cs_coils.append((cname, cvals))
        else:
            other_coils.append((cname, cvals))

    ax = axes[4]
    if cs_coils:
        cs_colors = plt.cm.tab10(np.linspace(0, 1, max(len(cs_coils), 1)))
        for ci, (cname, cvals) in enumerate(cs_coils):
            ct = sorted(cvals.keys())
            n_turns = _coil_net_turns(tt, cname)
            ci_vals = [cvals[t_v] * n_turns * 1e-6 for t_v in ct]
            ax.plot(ct, ci_vals, ls=LS_PRI, lw=LW * 0.8, color=cs_colors[ci], label=cname)
        ax.legend(fontsize=LEGEND_FS - 2, loc='upper left', ncol=2)
    else:
        ax.text(0.5, 0.5, 'No CS coils', transform=ax.transAxes,
                ha='center', va='center', fontsize=LABEL_FS)
    ax.set_ylabel('I_coil [MA-turns]', fontsize=LABEL_FS)
    _style(ax)

    ax = axes[5]
    if other_coils:
        oth_colors = plt.cm.tab20(np.linspace(0, 1, max(len(other_coils), 1)))
        for ci, (cname, cvals) in enumerate(other_coils):
            ct = sorted(cvals.keys())
            n_turns = _coil_net_turns(tt, cname)
            ci_vals = [cvals[t_v] * n_turns * 1e-6 for t_v in ct]
            ax.plot(ct, ci_vals, ls=LS_PRI, lw=LW * 0.75, color=oth_colors[ci], label=cname)
        ax.legend(fontsize=LEGEND_FS - 3, loc='upper left', ncol=2)
    else:
        ax.text(0.5, 0.5, 'No PF/other coils', transform=ax.transAxes,
                ha='center', va='center', fontsize=LABEL_FS)
    ax.set_ylabel('I_coil [MA-turns]', fontsize=LABEL_FS)
    _style(ax)

    for ax in np.ravel(axes):
        ax.axvline(t_now, color=VLINE_COLOR, ls=VLINE_LS, lw=VLINE_LW, zorder=10)


def _draw_profiles_movie(axes, tt, idx):
    r'''! Draw the 6 radial profile panels for a movie frame.'''
    s = tt._state

    ax = axes[0]
    ax.axhline(1.0, color='k', ls='-', lw=1, label='q=1')
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

    ax = axes[5]
    x, y = _prof(s.get('eta_prof', {}), idx)
    if x is not None:
        ax.plot(x, y, color=COLOR_TX, ls=LS_PRI, lw=LW, label='\u03b7 TX')
    ax.set_ylabel('\u03b7 [\u03a9\u00b7m]', fontsize=LABEL_FS)
    ax.set_yscale('log')
    ax.legend(fontsize=LEGEND_FS)
    _style(ax)


def _encode_video_from_dir(frame_dir, out_path, fps=2):
    r'''! Encode frames from a directory into an MP4 file.'''
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
    r'''! Embed an MP4 video in a Jupyter notebook.'''
    try:
        from IPython.display import Video, display
        display(Video(video_path, embed=True, mimetype='video/mp4'))
    except ImportError:
        print(f'Video saved to: {video_path}')


# ── Interactive widgets (Jupyter) ─────────────────────────────────────────────

def plot_profiles_interactive(tt):
    r'''! Interactive profile viewer with ipywidgets slider to scrub through timesteps.'''
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

            ax = axes[0, 0]
            if i in s.get('n_e', {}):
                ax.plot(s['n_e'][i]['x'], s['n_e'][i]['y'], 'b-')
            ax.set_title('n_e')
            ax.set_ylabel(r'$n_e$ [m$^{-3}$]')

            ax = axes[0, 1]
            if i in s.get('T_e', {}):
                ax.plot(s['T_e'][i]['x'], s['T_e'][i]['y'], 'r-')
            ax.set_title('T_e')
            ax.set_ylabel('T_e [keV]')

            ax = axes[0, 2]
            if i in s.get('q_prof_tm', {}):
                ax.plot(s['q_prof_tm'][i]['x'], s['q_prof_tm'][i]['y'], color=COLOR_TM, label='TM')
            if i in s.get('q_prof_tx', {}):
                ax.plot(s['q_prof_tx'][i]['x'], s['q_prof_tx'][i]['y'], color=COLOR_TX, label='TX')
            ax.set_title('q')
            ax.set_ylabel('q')
            ax.legend(fontsize=8)

            ax = axes[1, 0]
            for key, label, clr in [('j_tot', 'j_tot', 'k'), ('j_ohmic', 'j_ohm', 'r'),
                                     ('j_ni', 'j_ni', 'b'), ('j_bootstrap', 'j_BS', 'g')]:
                if i in s.get(key, {}):
                    ax.plot(s[key][i]['x'], s[key][i]['y'] / 1e6, color=clr, label=label)
            ax.set_title('Current densities')
            ax.set_ylabel('j [MA/m²]')
            ax.legend(fontsize=7)

            ax = axes[1, 1]
            if i in s.get('ptot', {}):
                ax.plot(s['ptot'][i]['x'], s['ptot'][i]['y'], color=COLOR_TX, label='p TX')
            if i in s.get('p_prof_tm', {}):
                ax.plot(s['p_prof_tm'][i]['x'], s['p_prof_tm'][i]['y'], color=COLOR_TM, ls='--', label='p TM')
            ax.set_title('Pressure')
            ax.set_ylabel('p [Pa]')
            ax.legend(fontsize=8)

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
    r'''! Equilibrium viewer — widget slider in notebook, saved PNGs otherwise.
    
        Parameters
        ----------
        notebook_mode : bool or None
            True  → ipywidgets slider (Jupyter).
            False → save PNG files to save_path.
            None  → auto-detect (widget if in Jupyter, save files otherwise).
        save_path : str, optional
            Directory for PNG files when notebook_mode=False.
            Defaults to ./equil_loop<N>/ in the current directory.
        
    '''
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


# ── Physics summary ───────────────────────────────────────────────────────────

def summary(tt):
    r'''! Print and return key physics figures of merit from the simulation.'''
    s = tt._state
    times = np.array(tt._tm_times)

    ft = getattr(tt, '_flattop', np.zeros(len(times), dtype=bool))
    ft_mask = ft.astype(bool)

    out = {}

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

    # out['Ip_max_MA'] = float(np.max(np.abs(s['Ip_tm'])) / 1e6)

    _, y_bN = _tx_scalar(tt, 'beta_N')
    if y_bN is not None:
        out['beta_N_max'] = float(np.nanmax(y_bN))
    # out['beta_N_tm_max'] = float(np.nanmax(s['beta_N_tm']))

    t_H98, y_H98 = _tx_scalar(tt, 'H98')
    if y_H98 is not None:
        out['H98_max'] = float(np.nanmax(y_H98))
        if np.any(ft_mask):
            ft_start, ft_end = times[ft_mask][0], times[ft_mask][-1]
            ft_h_mask = (t_H98 >= ft_start) & (t_H98 <= ft_end)
            out['H98_flattop_avg'] = float(np.nanmean(y_H98[ft_h_mask])) if np.any(ft_h_mask) else None

    _, y_Te = _tx_profile_at_rho(tt, 'T_e', 0.0)
    _, y_Ti = _tx_profile_at_rho(tt, 'T_i', 0.0)
    if y_Te is not None:
        out['T_e_core_max_keV'] = float(np.nanmax(y_Te))
    if y_Ti is not None:
        out['T_i_core_max_keV'] = float(np.nanmax(y_Ti))

    _, y_ne = _tx_scalar(tt, 'n_e_line_avg')
    if y_ne is not None:
        out['n_e_line_avg_max'] = float(np.nanmax(y_ne))

    out['f_GW_max'] = float(np.nanmax(s['f_GW']))
    out['q95_min'] = float(np.nanmin(s['q95_tm'][s['q95_tm'] > 0])) if np.any(s['q95_tm'] > 0) else None
    out['q0_min'] = float(np.nanmin(s['q0_tm'][s['q0_tm'] > 0])) if np.any(s['q0_tm'] > 0) else None

    psi_lcfs_tm = np.array(s['psi_lcfs_tm'])
    out['flux_consumed_Wb'] = float(-(psi_lcfs_tm[-1] - psi_lcfs_tm[0]) * 2 * np.pi)

    _, y_Pa = _tx_scalar(tt, 'P_alpha_total')
    _, y_Po = _tx_scalar(tt, 'P_ohmic_e')
    if y_Pa is not None:
        out['P_fusion_max_MW'] = float(np.nanmax(y_Pa) / 1e6 * 5)
    if y_Po is not None:
        out['P_ohmic_max_MW'] = float(np.nanmax(y_Po) / 1e6)

    out['l_i_flattop_avg'] = float(np.nanmean(s['l_i_tm'][ft_mask])) if np.any(ft_mask) else None

    if np.any(ft_mask):
        # out['vloop_tm_flattop_avg_V'] = float(np.nanmean(np.array(s['vloop_tm'])[ft_mask]))
        out['vloop_tx_flattop_avg_V'] = float(np.nanmean(np.array(s['vloop_tx'])[ft_mask]))

    print(f"\n{'=' * 55}")
    print("  TokaMaker_TORAX Physics Summary")
    print(f"{'=' * 55}")
    for key, val in out.items():
        if val is not None:
            if isinstance(val, float):
                print(f"  {key:30s}  {val:.4g}")
            else:
                print(f"  {key:30s}  {val}")
    print(f"{'=' * 55}\n")
    return out

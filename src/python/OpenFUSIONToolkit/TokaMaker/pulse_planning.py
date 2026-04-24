import copy
import json
import logging
import os
import platform
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

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, create_power_flux_fun

LCFS_WEIGHT = 100.0
N_PSI = 1000
_NBI_W_TO_MA = 1/16e6
mu_0 = 4.0 * np.pi * 1e-7

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
        'bootstrap_current': {},
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
        'smoothing_width': 0.3,
        'DV_effective': True,
        'include_ITG': True,
        'include_TEM': True,
        'include_ETG': True,
        'avoid_big_negative_s': False,
    },
    'solver': {
        'solver_type': 'newton_raphson',
        'use_predictor_corrector': True,
        'n_corrector_steps': 10,
        'chi_pereverzev': 30,
        'D_pereverzev': 15,
        'use_pereverzev': True,
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

class TokTox:
    '''! TokaMaker + TORAX Coupled Pulse Simulation Code'''


    # ─── Initialization ─────────────────────────────────────────────────────────

    def __init__(self, t_init, t_final, eqtimes, g_eqdsk_arr, tx_dt=0.1, tm_times=None, last_surface_factor=0.99, cocos=2, oft_env=None, oft_threads=4, truncate_eq=False):
        r'''! Initialize the Coupled TokaMaker + TORAX object.
        @param t_init Start time (s).
        @param t_final End time (s).
        @param eqtimes Time points of each gEQDSK file.
        @param g_eqdsk_arr Filenames of each gEQDSK file.
        @param tx_dt Time step (s) of TORAX simulation.
        @param tm_times Time points where TokaMaker solves equilibrium.
        @param last_surface_factor Last surface factor for Torax.
        @param n_rho Number of grid cells for torax.
        @param cocos COCOS version of input EQDSK.
        @param oft_env OFT environment.
        @param truncate_eq Whether to truncate equilibrium when saving TokaMaker output to EQDSK.
        @param oft_threads Number of threads for OFT.
        '''
        if oft_env is not None:
            self._oftenv = oft_env
        else:
            self._oftenv = OFT_env(nthreads=oft_threads)
        self._tm = TokaMaker(self._oftenv)
        self._cocos = cocos

        self._state = {}
        self._eqtimes = eqtimes
        self._results = {}
        self._init_files = g_eqdsk_arr
        self._t_init = t_init
        self._t_final = t_final
        self._tx_dt = tx_dt # TORAX timestep
        self._last_surface_factor = last_surface_factor
        self._n_rho = 50 # resolution of TORAX grid, default 50, changed with set_tx_grid()
        self._psi_N = np.linspace(0.0, 1.0, N_PSI) # standardized psi_N grid all values should be mapped onto
        self._truncate_eq = truncate_eq

        self._current_loop = 0

        if tm_times is None:
            self._tm_times = sorted(eqtimes)
        else:
            self._tm_times = sorted(tm_times)

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

        def interp_lcfs(lcfs, n=1000):
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

            lcfs.append(interp_lcfs(g['rzout']))

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
        
        self._Ip = None
        self._Zeff = None

        self._nbi_heating = None
        self._ecrh_heating = None
        self._ecrh_loc = None
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
        # load_config().
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

        self._loaded_config = None   # set by load_config()
        self._tx_grid_type = None
        self._tx_grid = None

        self._diverted_times = None
        self._x_point_targets = None
        self._x_point_weight = 100.0
        self._strike_point_targets = None

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

    def output_dir(self):
        r'''! Return output directory for the latest run, or None when output_mode=False.'''
        return self._out_dir

    # ─── Static Utilities ───────────────────────────────────────────────────────

    @staticmethod
    def _to_plain_python(obj):
        r'''! Recursively convert numpy scalars/arrays to plain Python types.

        Used before pformat-saving config dicts so the saved .py files are
        loadable without numpy (no ``array([...])`` references).
        '''
        if isinstance(obj, dict):
            return {TokTox._to_plain_python(k): TokTox._to_plain_python(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            converted = [TokTox._to_plain_python(v) for v in obj]
            return type(obj)(converted)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    @staticmethod
    def _config_merge(base, override):
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
                TokTox._config_merge(base[key], val)
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
                    TokTox._flatten_time_dependent(val)
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


    # ─── Setup & Configuration ──────────────────────────────────────────────────

    def load_config(self, config):
        r'''! Load a TORAX config dict.

        The loaded config is deep-merged on top of BASE_CONFIG when the
        simulation config is built.  Any key present in the loaded config
        will override the corresponding BASE_CONFIG key; keys only in
        BASE_CONFIG are kept as-is.  Geometry is always overwritten by
        TokTox (eqdsk-based).

        Explicit ``set_*()`` calls made AFTER ``load_config()`` will
        override both the base and the loaded config.

        @param config Dictionary (TORAX config format).
        '''
        self._loaded_config = copy.deepcopy(config)

    def set_tx_grid(self, grid_type, grid):
        r'''! Set TORAX grid type and grid points.
        @param grid_type Grid type ('n_rho' or 'face_centers').
        @param grid Grid points (integer or np.array).
        '''
        self._tx_grid_type = grid_type
        self._tx_grid = grid
        if grid_type not in ['n_rho', 'face_centers']:
            raise ValueError(f'Invalid grid type: {type}. Must be "n_rho" or "face_centers".')

    def initialize_tm(self, mesh, R0_geo, weights=None, vsc=None):
        r'''! Initialize GS Solver Object.
        @param mesh Filename of reactor mesh.
        @param R0_geo Major radius of machine geometric center.
        @param vsc Vertical Stability Coil.
        '''
        mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh(mesh)
        self._tm.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        self._tm.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
        self._tm.setup(order = 2, F0 = R0_geo*self._state['B0'][0])

        self._tm.settings.maxits = 100

        if vsc is not None:
            self._tm.set_coil_vsc({vsc: 1.0})

    def set_coil_reg(self, coil_bounds=None, updownsym=False,
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
        self._apply_coil_reg(targets=None)

    def _apply_coil_reg(self, targets=None):
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
        r'''! Set plasma current (Amps).
        @param ip Plasma current.
        '''
        self._Ip = Ip

    def set_ne(self, n_e):
        r'''! Set density profiles.
        @param n_e Electron density (m^-3).
        '''
        self._n_e = n_e

    def set_Te(self, T_e):
        r'''! Set electron temperature profiles (keV).
        @param T_e Electron temperature.
        '''
        self._T_e = T_e

    def set_Ti(self, T_i):
        r'''! Set ion temperature profiles (keV).
        @param T_i ion temperature.
        '''
        self._T_i = T_i

    def set_Zeff(self, Zeff):
        r'''! Set plasma effective charge.
        @param z_eff Effective charge.
        '''
        self._Zeff = Zeff

    def set_plasma_composition(self, main_ion=None, impurity=None):
        r'''! Set plasma composition (fuel and impurity species).

        Must be called together with set_Zeff() — set_Zeff() provides the
        Z_eff target which TORAX uses to determine the impurity density for the
        species specified here.

        @param main_ion Main ion species dict, e.g. {'D': 0.5, 'T': 0.5} for DT.
        @param impurity Impurity species string, e.g. 'Ne', 'Ar', 'W'.
        '''
        if main_ion is not None:
            self._main_ion = main_ion
        if impurity is not None:
            self._impurity = impurity

    def set_sources(self, fusion=True, ei_exchange=True):
        r'''! Enable standard TORAX physics sources using TORAX defaults.

        Each flag adds the corresponding source with an empty config dict so
        TORAX uses its built-in defaults.  Only call the ones you need — empty
        dicts pull in TORAX's machine-specific defaults (e.g. ITER geometry for
        fusion), which may not be appropriate for every simulation.
        
        Ohmic heating is enabled by default in the base_config. 
        Currently there is no way to disable ohmic heating, other than changing base_config.

        @param fusion    Enable fusion alpha heating.
        @param ei_exchange Enable electron-ion energy exchange.
        '''
        self._enable_fusion = fusion
        self._enable_ei_exchange = ei_exchange

    def set_nbar(self, nbar, normalize_to_nbar=True):
        r'''! Set line averaged density over time.
        @param nbar Density (m^-3).
        @param normalize_to_nbar Whether to normalize initial n_e profile to match nbar.
        '''
        self._nbar = nbar
        self._normalize_to_nbar = normalize_to_nbar # when True, initial n_e profile will be normalized to match nbar, but time evolution will not be constrained to match nbar

    def set_right_bc(self, ne_right_bc=None, Te_right_bc=None, Ti_right_bc=None):
        if ne_right_bc:
            self._ne_right_bc = ne_right_bc
        if Te_right_bc:
            self._Te_right_bc = Te_right_bc
        if Ti_right_bc:
            self._Ti_right_bc = Ti_right_bc

    def set_heating(self, generic_heat=None, generic_heat_loc=None, generic_heat_width=0.25, nbi_current=False, ecrh=None, ecrh_loc=None, ecrh_width=0.1, ohmic=None):
        r'''! Set heating sources for Torax.

        Ohmic heating is always enabled (it is on by default in BASE_CONFIG).
        @param generic_heat Generic heating (dictionary of {time: power_in_watts}).
        @param generic_heat_loc Generic heating deposition location (normalized rho).
        @param generic_heat_width Generic heating deposition width (normalized rho).
        @param nbi_current Use NBI current.
        @param ecrh ECRH heating (dictionary of {time: power_in_watts}).
        @param ecrh_loc ECRH deposition location (normalized rho).
        @param ecrh_width ECRH deposition width (normalized rho).
        @param nbi_current Whether to include NBI current drive, uses _NBI_W_TO_MA = 1/16e6 to convert heating to current driven.
        '''
        if generic_heat is not None and generic_heat_loc is not None:
            self._generic_heat = generic_heat
            self._generic_heat_loc = generic_heat_loc
            self._generic_heat_width = generic_heat_width
        if ecrh is not None and ecrh_loc is not None:
            self._ecrh_heating = ecrh
            self._ecrh_loc = ecrh_loc
            self._ecrh_width = ecrh_width
        if ohmic is not None:
            self._ohmic_power = ohmic
        
        self._use_nbi_current = nbi_current

    def set_pedestal(self, set_pedestal=False, T_i_ped=None, T_e_ped=None, n_e_ped=None, ped_top=0.95):
        r'''! Set pedestals for ion/electron temperatures and density (legacy API).

        This preserves the original behavior:
        - set_pedestal=True uses model_name='set_T_ped_n_ped'
        - set_pedestal=False uses model_name='no_pedestal'

        @param set_pedestal Toggle pedestal model on/off.
        @param T_i_ped Ion temperature pedestal (time-varying scalar allowed).
        @param T_e_ped Electron temperature pedestal (time-varying scalar allowed).
        @param n_e_ped Electron density pedestal (time-varying scalar allowed).
        @param ped_top Pedestal-top location rho_norm_ped_top.
        '''
        # Using legacy setter disables full pedestal dict replacement.
        self._pedestal_config = None

        self._set_pedestal = set_pedestal
        if T_i_ped is not None:
            self._T_i_ped = T_i_ped
        if T_e_ped is not None:
            self._T_e_ped = T_e_ped
        if n_e_ped is not None:
            self._n_e_ped = n_e_ped
        self._ped_top = ped_top

    def load_pedestal_config(self, pedestal_config):
        r'''! Load pedestal config dict and fully replace 'pedestal' in TORAX config.

        If provided, this dict is copied to 'myconfig['pedestal']' directly,
        replacing the full pedestal section from BASE_CONFIG and load_config().

        @param pedestal_config Dictionary in TORAX pedestal config format, or
                               None to clear and fall back to set_pedestal().
        '''
        if pedestal_config is None:
            self._pedestal_config = None
            return
        if not isinstance(pedestal_config, dict):
            raise TypeError('pedestal_config must be a dictionary or None.')
        self._pedestal_config = copy.deepcopy(pedestal_config)

    def set_evolve(self, density=True, Ti=True, Te=True, current=True):
        r'''! Set variables as either prescribed (False) or evolved (True).
        @param density Evolve density.
        @param Ti Evolve ion temperature.
        @param Te Evolve electron temperature.
        @param current Evolve current.
        '''
        self._evolve_density = density
        self._evolve_current = current
        self._evolve_Ti = Ti
        self._evolve_Te = Te

    def set_gas_puff(self, S_total=None, decay_length=None):
        r'''! Set gas puff particle source.
        @param S_total Particle source (particles/s).
        @param decay_length Decay length from edge (normalized rho coordinates).
        '''
        self._gp_s = S_total
        self._gp_dl = decay_length

    def set_pellet(self, pellet_deposition_location=None, pellet_width=None, S_total=None):
        r'''! Set pellet fueling particle source.
        @param pellet_deposition_location Pellet deposition location (normalized rho coordinates).
        @param pellet_width Pellet deposition width (normalized rho coordinates).
        @param S_total Particle source (particles/s).
        '''
        self._pellet_deposition_location = pellet_deposition_location
        self._pellet_width = pellet_width
        self._pellet_s_total = S_total

    def set_chi(self, chi_min=None, chi_max=None):
        if chi_min is not None:
            self._chi_min = chi_min
        if chi_max is not None:
            self._chi_max = chi_max

    def set_De(self, De_min=None, De_max=None):
        if De_min is not None:
            self._De_min = De_min
        if De_max is not None:
            self._De_max = De_max

    def set_Ve(self, Ve_min=None, Ve_max=None):
        if Ve_min is not None:
            self._Ve_min = Ve_min
        if Ve_max is not None:
            self._Ve_max = Ve_max

    def set_x_points(self, diverted_times=None, x_point_targets=None, x_point_weight=100.0):
        r'''! Configure diverted window and X-point targets for TM saddle constraints.

        @param diverted_times Tuple (t_start, t_end) defining the diverted plasma window.
        @param x_point_targets X-point target locations, shape (n_xpoints, 2) with [R, Z] pairs.
        @param x_point_weight Weight for saddle-point constraints.
        '''
        if diverted_times is not None and len(diverted_times) != 2:
            raise ValueError('diverted_times must be a (t_start, t_end) tuple.')

        self._diverted_times = diverted_times
        self._x_point_targets = None if x_point_targets is None else np.atleast_2d(x_point_targets)
        self._x_point_weight = x_point_weight

    def set_strike_points(self, strike_point_targets):
        r'''! Manually specify strike point locations to add as isoflux targets during the diverted phase.

        Strike points are where the separatrix legs intersect the divertor/limiter surface.
        They are added to the isoflux constraint array during the diverted window alongside
        the LCFS shape targets. Automatic detection from the EQDSK boundary is not possible
        because rzout traces only the closed plasma boundary (inside the vessel).

        @param strike_point_targets Strike point locations, shape (n_points, 2) with [R, Z] pairs,
                                    or None to disable.
        '''
        self._strike_point_targets = None if strike_point_targets is None else np.atleast_2d(strike_point_targets)


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

    def _interp_profile_onto_psi(self, data_tree, var_name, time, profile_type='linterp'):
        r'''! Interpolate a single TORAX profile snapshot onto self._psi_N.

        This is the inner workhorse — no averaging, just one timeslice.
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
            data_on_psi = self._interp_profile_onto_psi(data_tree, var_name, time, profile_type)
        else:
            # Collect all TORAX timesteps inside the window
            mask = (tx_times >= t_start) & (tx_times <= t_end)
            win_times = tx_times[mask]
            if len(win_times) == 0:
                data_on_psi = self._interp_profile_onto_psi(data_tree, var_name, time, profile_type)
            else:
                stack = np.stack([
                    self._interp_profile_onto_psi(data_tree, var_name, wt, profile_type)
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


    def _get_tx_config(self):
        r'''! Generate config object for Torax simulation.

        Build order
        -----------
        1. Deep-copy BASE_CONFIG (from baseconfig.py).
        2. Deep-merge the loaded config on top (if load_config() was called).
           Every key in the loaded config overwrites the matching base key;
           keys only in BASE_CONFIG are kept as-is.
        3. Override geometry (always set by TokTox / TokaMaker equilibria).
        4. Override t_initial / t_final / fixed_dt from __init__ params.
        5. Use psi profile from loop 0 (if available) from profile_conditions.
        6. Apply any explicit set_*() overrides (only when the attribute is
           not None, i.e. the user called the setter after load_config).

        @return Torax config object.
        '''

        # ── 1. Start from base config ──────────────────────────────────────
        myconfig = copy.deepcopy(BASE_CONFIG)

        # ── 2. Deep-merge loaded config ────────────────────────────────────
        if self._loaded_config is not None:
            self._config_merge(myconfig, self._loaded_config)

        # ── 3. Geometry (always set by TokTox) ─────────────────────────────
        myconfig['geometry'] = {
            'geometry_type': 'eqdsk',
            'geometry_directory': os.getcwd(),
            'last_surface_factor': self._last_surface_factor,
            'n_surfaces': 50,
            'Ip_from_parameters': True, # True tells TX to pull from config, not from eqdsk, in case eqdsks fail TX retains correct Ip targets
        }
        if self._current_loop == 1:
            eq_safe = []
            t_safe = []
            t_skipped = []
            for i, t in enumerate(self._eqtimes):
                eq = self._init_files[i]
                if self._test_eqdsk(eq):
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
            for i, t in enumerate(self._tm_times):
                eqdsk = os.path.join(self._eqdsk_dir, f'{self._current_loop - 1:03d}.{i:03d}.eqdsk')
                tm_ok = (eqdsk not in self._eqdsk_skip) and self._test_eqdsk(eqdsk)
                if tm_ok:
                    full_eqdsk_map[t] = eqdsk
                    n_tm += 1
            # If i=0 TM failed, always fall back to the seed EQDSK so TORAX
            # has a valid initial geometry. Other failed timesteps are left out of the
            # map and TORAX interpolates from neighboring solved entries.
            t0 = self._tm_times[0]
            if t0 not in full_eqdsk_map:
                seed_eqdsk = self._init_files[0]
                if self._test_eqdsk(seed_eqdsk):
                    full_eqdsk_map[t0] = seed_eqdsk
                    self._log(f'Loop {self._current_loop}: TM failed at t=0, falling back to seed EQDSK for t=0.')
                else:
                    self._log(f'Warning: Loop {self._current_loop}: TM failed at t=0 and seed EQDSK is also invalid.')

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

        # ── 5. Psi profile from loop 0  ──────────
        myconfig.setdefault('profile_conditions', {})
        if self._psi_init is not None:
            myconfig['profile_conditions']['psi'] = self._psi_init
            myconfig['profile_conditions']['initial_psi_mode'] = 'profile_conditions'
            myconfig['profile_conditions']['initial_psi_from_j'] = False
        else:
            myconfig['profile_conditions']['initial_psi_mode'] = 'geometry' # if loop 0 wasn't run, uses psi from initial eqdsk, not ideal

        # ── 6. Explicit set_*() overrides ──────────────────────────────────
        #    Only applied when the attribute is not None (i.e. the user called the setter
        #    explicitly after load_config(); None means fall through to the loaded/base config).
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

            if self._use_nbi_current:
                myconfig['sources'].setdefault('generic_current', {})
                myconfig['sources']['generic_current']['use_absolute_current'] = True
                myconfig['sources']['generic_current']['I_generic'] = (nbi_times, _NBI_W_TO_MA * np.array(nbi_pow))
                myconfig['sources']['generic_current']['gaussian_location'] = self._generic_heat_loc

        if self._pedestal_config is not None:
            # Full pedestal dict replacement requested via load_pedestal_config().
            myconfig['pedestal'] = copy.deepcopy(self._pedestal_config)
        else:
            # Accept bool or time-varying dict for set_pedestal.
            ped_enabled = (
                self._set_pedestal is not None
                and not (
                    isinstance(self._set_pedestal, (bool, np.bool_))
                    and (not bool(self._set_pedestal))
                )
            )

            if ped_enabled:
                # Build in required key order:
                # model_name -> set_pedestal -> rho_norm_ped_top -> T_i/T_e/n_e.
                ped_cfg = {}
                ped_cfg['model_name'] = 'set_T_ped_n_ped'
                ped_cfg['set_pedestal'] = self._set_pedestal
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

        if self._output_mode in ('normal', 'debug') and self._out_dir is not None:
            cfg_name = f'tx_config_loop{self._current_loop:03d}.py'
            if self._output_file_tag is not None:
                cfg_name = f'{self._output_file_tag}_{cfg_name}'
            config_filename = os.path.join(self._out_dir, cfg_name)
            with open(config_filename, 'w') as f:
                f.write('# Torax configuration\n')
                f.write(f'# Loop {self._current_loop}\n\n')
                f.write('tx_config = ')
                f.write(pprint.pformat(self._to_plain_python(myconfig), width=100, sort_dicts=False))

        tx_config = torax.ToraxConfig.from_dict(myconfig)
        return tx_config

    def _test_eqdsk(self, eqdsk):
            myconfig = copy.deepcopy(BASE_CONFIG)
            if self._loaded_config is not None:
                self._config_merge(myconfig, self._loaded_config)
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
                self._log(f"TEST EQDSK FAILED: {eqdsk} — {repr(e)}")
                self._print(f'    EQDSK rejected by TORAX: {os.path.basename(eqdsk)}')
                return False

    def _run_tx_init(self):
        r'''! Loop 0: Run a short TORAX simulation with eqdsk geometry to equilibrate initial inputs.

        Run TORAX for 1 second with steady-state (time-flattened) inputs lets TORAX
        evolve them to a more physical state.  The relaxed values are then
        injected into the config used by the main simulation (loop 1 onwards).
        '''
        INIT_RUNTIME = 0.5
        INIT_DT = 0.01
        self._log('Transport init: building steady-state init config...')

        init_config = copy.deepcopy(BASE_CONFIG)
        if self._loaded_config is not None:
            self._config_merge(init_config, self._loaded_config)

        # eqdsk geometry from the first seed file — same geometry as loop 1, i=0.
        # This ensures the psi evolved here satisfies the same GS metric coefficients
        # as the main sim, so injected psi produces smooth j at t=0 of loop 1.
        init_eqdsk = self._init_files[0]
        if not self._test_eqdsk(init_eqdsk):
            raise ValueError(f'Transport init: first seed eqdsk is not valid: {init_eqdsk}')
        init_config['geometry'] = {
            'geometry_type': 'eqdsk',
            'geometry_directory': os.getcwd(),
            'last_surface_factor': self._last_surface_factor,
            'n_surfaces': 50,
            'Ip_from_parameters': True,
            'geometry_configs': {self._t_init: {'geometry_file': init_eqdsk, 'cocos': self._cocos}},
        }

        if self._tx_grid_type == 'n_rho':
            init_config['geometry']['n_rho'] = self._tx_grid
        elif self._tx_grid_type == 'face_centers':
            init_config['geometry']['face_centers'] = self._tx_grid

        # Numerics: 1-second steady-state init sim
        init_config.setdefault('numerics', {})
        init_config['numerics']['t_initial'] = self._t_init
        init_config['numerics']['t_final'] = self._t_init + INIT_RUNTIME
        init_config['numerics']['fixed_dt'] = INIT_DT
        
        init_config['numerics']['evolve_current'] = True # Let current evolve to relax to psi profile
        init_config['numerics']['evolve_density'] = False # Fix ne and Te/Ti profiles
        init_config['numerics']['evolve_ion_heat'] = False
        init_config['numerics']['evolve_electron_heat'] = False

        # Explicitly use J mode (nu formula) to initialize psi in the loop-0 sim.
        # This is the same mode the main sim (loop 1+) will use, so T_e/T_i
        # injected from here will be self-consistent with the initial current profile.
        init_config.setdefault('profile_conditions', {})
        init_config['profile_conditions']['initial_psi_mode'] = 'geometry'

        # Propagate Ip from set_Ip() into the Loop 0 config.
        # Without this, TORAX runs without an Ip constraint and may produce
        # psi values with the wrong magnitude or sign (observed: ~10× too large
        # for ITER). The loaded_config Ip (if any) is already merged above;
        # set_Ip() values are stored separately in self._Ip and must be applied here.
        if self._Ip is not None:
            init_config['profile_conditions']['Ip'] = copy.deepcopy(self._Ip)

        # Flatten all time-dependent values to their initial value (steady-state inputs)
        self._flatten_time_dependent(init_config)

        if self._output_mode in ('normal', 'debug') and self._out_dir is not None:
            cfg_name = 'tx_config_loop000.py'
            if self._output_file_tag is not None:
                cfg_name = f'{self._output_file_tag}_{cfg_name}'
            config_filename = os.path.join(self._out_dir, cfg_name)
            with open(config_filename, 'w') as f:
                f.write('# Torax configuration\n# Loop 0 (transport init)\n\n')
                f.write('tx_config = ')
                f.write(pprint.pformat(self._to_plain_python(init_config), width=100, sort_dicts=False))

        self._log(f'Transport init: running ~{INIT_RUNTIME}s steady-state TORAX simulation...')
        tx_config = torax.ToraxConfig.from_dict(init_config)
        data_tree, hist = torax.run_simulation(tx_config, log_timestep_info=False)

        if hist.sim_error != torax.SimError.NO_ERROR:
            raise ValueError(f'Transport init simulation failed: {hist.sim_error}')

        t_final_init = self._t_init + INIT_RUNTIME

        # Propagate psi from init into main simulation config
        main_config = self._loaded_config if self._loaded_config is not None else BASE_CONFIG
        main_config.setdefault('profile_conditions', {})

        # Extract psi directly on its own grid
        psi_xr = data_tree.profiles.psi.sel(time=t_final_init, method='nearest')
        rho_psi_arr = psi_xr.coords['rho_norm'].to_numpy()
        psi_arr = psi_xr.to_numpy()
        self._psi_init = ([self._t_init], rho_psi_arr.tolist(), [psi_arr.tolist()])
        


    def _run_tx(self):
        r'''! Run the TORAX transport simulation.
        @return Tuple (consumed_flux, consumed_flux_integral).
        '''
        self._print('  TORAX: running simulation...')
        myconfig = self._get_tx_config()
        try:
            data_tree, hist = torax.run_simulation(myconfig, log_timestep_info=False)
        except Exception as e:
            self._print(f'  TORAX: config/init FAILED — {e}')
            raise

        if hist.sim_error != torax.SimError.NO_ERROR:
            self._print(f'  TORAX: sim FAILED ({hist.sim_error})')
            raise ValueError(f'TORAX failed to run the simulation: {hist.sim_error}')
        
        self._data_tree = data_tree  # store for visualization at full TORAX resolution

        v_loops = np.zeros(len(self._tm_times))
        for i, t in enumerate(self._tm_times):
            self._tx_update(i, data_tree)
            v_loops[i] = data_tree.scalars.v_loop_lcfs.sel(time=t, method='nearest')

        if self._save_outputs:
            self._res_update(data_tree)

        # Flux consumption: positive Ip drives psi_lcfs down in TM-native.
        # Convention: consumed_flux > 0 means plasma has consumed flux from CS.
        consumed_flux = -2.0 * np.pi * (self._state['psi_lcfs_tx'][-1] - self._state['psi_lcfs_tx'][0])
        consumed_flux_integral = np.trapezoid(v_loops[0:], self._tm_times[0:])
        self._log(f"Loop {self._current_loop} TORAX: cflux={consumed_flux:.4f} Wb")
        self._print(f'  TORAX: done (cflux={consumed_flux:.4f} Wb)')
        return consumed_flux, consumed_flux_integral

    def _tx_update(self, i, data_tree):
        r'''! Update the simulation state from TORAX results at timestep i.

        All profile and scalar extractions use time-averaged methods to smooth
        sawtooth oscillations when averaging is enabled.

        @param i Timestep index.
        @param data_tree Result object from Torax.
        '''
        t = self._tm_times[i]

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

        ffp_ni = self._calc_ffp_ni(i)
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
        # Transform: psi_TM = -psi_COCOS11 / (2π)
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

    def _calc_ffp_ni(self, i):
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
        from tqdm import tqdm
        self._print(f'  TokaMaker: solving {len(self._tm_times)} equilibria...')
        self._log(f"Loop {self._current_loop} TokaMaker:")

        self._eqdsk_skip = []
        _loop_level_log = []

        # ── Per-loop initialization (before timestep sweep) ──────────────────
        self._state['psi_grid_prev_tm'] = {}

        # Seed coil regularization targets
        # Loop 1: use zero targets (None) so the solver freely finds the correct
        # coil configuration without being biased toward the initial equilibrium.
        # Loop 2+: seed from the last solve of the previous loop for warm-starting.
        if getattr(self, '_coil_reg_config', None):
            init_targets = None
            if self._current_loop > 1:
                try:
                    init_targets, _ = self._tm.get_coil_currents()
                except Exception as e:
                    self._log(f'TM: could not read initial equilibrium coil currents: {e}')
            self._apply_coil_reg(targets=init_targets)

        # Debug: log coil bounds at the start of each loop
        if self._debug_mode and hasattr(self, '_coil_bounds') and self._coil_bounds:
            self._log('  TM coil bounds [A/turn] (turns * A/turn = total A-turns):')
            for cname, (lo, hi) in self._coil_bounds.items():
                n_turns = self._tm.coil_sets.get(cname, {}).get('net_turns', 1.0)
                self._log(f'    {cname}: [{lo:.3g}, {hi:.3g}] A/turn  x {n_turns:.0f} turns = [{lo*n_turns:.3g}, {hi*n_turns:.3g}] A-turns')


        if 0 in self._psi_warm_start and self._psi_warm_start[0] is not None:
            self._tm.set_psi_dt(psi0=self._psi_warm_start[0], dt=1.0e10)

        _pbar = tqdm(enumerate(self._tm_times), total=len(self._tm_times),
                    desc=f'  TM loop {self._current_loop}', unit='solve',
                    bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_inv_fmt}]{postfix}')
        for i, t in _pbar:
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
                perc_limit = 0.60       # LCFS points above percentage limit* max(abs(Z)) are removed from isoflux targets
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

            if i > 0:
                if self._state['psi_grid_prev_tm'][i-1] is not None:
                    self._tm.set_psi_dt(psi0=self._state['psi_grid_prev_tm'][i-1], dt=self._tm_times[i]-self._tm_times[i-1])
            
            skip_coil_update = False
            eq_name = os.path.join(self._eqdsk_dir, f'{self._current_loop:03d}.{i:03d}.eqdsk')

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
            ffp_1, pp_1 = self._level1_raw(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
            level_profiles.append({'ffp': ffp_1, 'pp': pp_1, 'name': 'lv1: raw'})
            
            # Level 2: sign flip
            ffp_2, pp_2 = self._level2_sign_flip(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
            level_profiles.append({'ffp': ffp_2, 'pp': pp_2, 'name': 'lv2: sign_flip'})
            
            # Level 3: pedestal smoothing (takes p_profile as input) # TODO: read in actual n_rho_ped_top, have to add to state first
            ffp_3, pp_3 = self._level3_pedestal_smoothing(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw), copy.deepcopy(self._state['p_prof_tx'][i])) 
            level_profiles.append({'ffp': ffp_3, 'pp': pp_3, 'name': 'lv3: pedestal_smoothing'})
            
            # Level 4: power flux
            ffp_4, pp_4 = self._level4_power_flux(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
            level_profiles.append({'ffp': ffp_4, 'pp': pp_4, 'name': 'lv4: power_flux'})

            # Level 5: power flux + pax from initial eqdsk
            ffp_5, pp_5 = self._level4_power_flux(copy.deepcopy(ffp_prof_raw), copy.deepcopy(pp_prof_raw))
            level_profiles.append({'ffp': ffp_5, 'pp': pp_5, 'name': 'lv5: power_flux + pax'})

            # Try each level
            for level_idx, level_prof in enumerate(level_profiles):
                level_name = level_prof['name']
                ffp_level = level_prof['ffp']
                pp_level = level_prof['pp']

                # Initialize psi from geometry parameters
                self._tm.init_psi(self._state['R0_mag'][i], self._state['Z'][i], self._state['a'][i], self._state['kappa'][i], self._state['delta'][i])

                if i in self._psi_warm_start and self._psi_warm_start[i] is not None:
                    self._tm.set_psi(self._psi_warm_start[i], update_bounds=True)
                elif i > 0 and (i-1) in self._state.get('psi_grid_prev_tm', {}) and self._state['psi_grid_prev_tm'][i-1] is not None:
                    self._tm.set_psi(self._state['psi_grid_prev_tm'][i-1], update_bounds=True)

                try:
                    with self._quiet_tm():
                        self._tm.set_profiles(ffp_prof=ffp_level, pp_prof=pp_level,
                                              ffp_NI_prof=self._state['ffp_ni_prof'][i])
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

            if not solve_succeeded:
                self._eqdsk_skip.append(eq_name)
                skip_coil_update = True
                self._log(f'\tTM: Solve failed at t={t} (all levels attempted).')
                self._state['psi_grid_prev_tm'][i] = None  # if solve failed, set psi grid to None
            
            if solve_succeeded:
                with self._quiet_tm():
                    self._state['equil'][i].save_eqdsk(eq_name,
                        lcfs_pad=1-self._last_surface_factor, run_info='TokaMaker EQDSK',
                        cocos=self._cocos, nr=300, nz=300, truncate_eq=self._truncate_eq)
                    self._tm_update(i)

                # Store diverted/limited flag for this timestep
                if not hasattr(self, '_diverted_flags'):
                    self._diverted_flags = {}
                self._diverted_flags[i] = self._state['equil'][i].diverted

                # Store psi on nodes for later movie generation
                self._tm_psi_on_nodes.setdefault(self._current_loop, {})[i] = self._state['equil'][i].get_psi(normalized=False)

            _winning = next((a for a in level_attempts if a['succeeded']), None)
            _last_attempt = level_attempts[-1] if level_attempts else {}
            _loop_level_log.append({
                'i': i, 't': t,
                'succeeded': solve_succeeded,
                'level': _winning['level'] if _winning else None,
                'level_name': _winning['name'] if _winning else None,
                'error': _last_attempt.get('error') if not solve_succeeded else None,
            })

            if self._output_mode in ('debug', 'normal') and self._out_dir is not None:
                if self._output_mode == 'debug':
                    diag_name = f'tm_diag_loop{self._current_loop:03d}_tidx{i:03d}.png'
                    if self._output_file_tag is not None:
                        diag_name = f'{self._output_file_tag}_{diag_name}'
                    _diag_path = os.path.join(self._out_dir, diag_name)
                    try:
                        tm_diagnostic_plot(self, i, t, level_attempts, solve_succeeded,
                                           save_path=_diag_path, display=False)
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

            # Update progress bar postfix; print FAIL messages above the bar
            if solve_succeeded:
                lvl = _winning['level']
                _pbar.set_postfix_str(f't={t:.2f}s OK(L{lvl})', refresh=False)
            else:
                err_short = (_last_attempt.get('error') or 'unknown')[:60]
                tqdm.write(f'    WARNING: TM FAIL at t={t:.2f}s — {err_short}')
                self._log(f'    TM FAIL at t={t:.2f}s — {err_short}')
                _pbar.set_postfix_str(f't={t:.2f}s FAIL', refresh=False)

            if not skip_coil_update and getattr(self, '_coil_reg_config', None):
                prev_coil_targets, _ = self._state['equil'][i].get_coil_currents()
                self._apply_coil_reg(targets=prev_coil_targets)

        # Flux consumption: positive Ip drives psi_lcfs down in TM convention
        # Convention: consumed_flux > 0
        consumed_flux = -(self._state['psi_lcfs_tm'][-1] - self._state['psi_lcfs_tm'][0]) * 2.0 * np.pi
        consumed_flux_integral = np.trapezoid(self._state['vloop_tm'][0:], self._tm_times[0:])

        n_ok = sum(1 for e in _loop_level_log if e['succeeded'])
        self._print(f'  TokaMaker: {n_ok}/{len(self._tm_times)} solved (cflux={consumed_flux:.4f} Wb)')

        # Compact level-usage summary for log
        from collections import Counter
        _lvl_counts = Counter(e['level_name'] for e in _loop_level_log if e['succeeded'])
        _lvl_summary = ', '.join(f'{name}: {cnt}' for name, cnt in sorted(_lvl_counts.items()))
        n_fail = len(self._tm_times) - n_ok
        self._log(f'\tTM summary: {n_ok}/{len(self._tm_times)} solved. Levels: {_lvl_summary}.'
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

        return consumed_flux, consumed_flux_integral
        
    # ── Profile level functions ──────────────────────────────────────────
    # Each level takes (self, ffp_prof, pp_prof, i) and returns (ffp_prof, pp_prof).
    # All levels receive deep copies of the raw TORAX profiles (not cumulative).

    def _level1_raw(self, ffp_prof, pp_prof):
        r'''! Raw TORAX profiles passed through unchanged.'''
        return ffp_prof, pp_prof

    def _level2_sign_flip(self, ffp_prof, pp_prof):
        r'''! Sign-flip clipping: clip each profile to its dominant sign.'''
        def _clip(prof):
            y = prof['y']
            sign = 1 if np.sum(y > 0) >= np.sum(y < 0) else -1
            y_new = np.clip(y, 0, None) if sign > 0 else np.clip(y, None, 0)
            return {**prof, 'y': y_new}
        return _clip(ffp_prof), _clip(pp_prof)

    def _level3_pedestal_smoothing(self, ffp_prof, pp_prof, p_prof, transition_psi_N = 0.6, gauss_sigma=8, blend_width=0.02, sav_window=41, sav_order=3):
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

    def _level4_power_flux(self, ffp_prof, pp_prof):
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
        psi_geo, q_tm, geo, _, _, _ = self._state['equil'][i].get_q(npsi=N_PSI, psi_pad=0.02)

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
    #  fly — main simulation loop
    # =========================================================================


    # ─── Main Simulation Loop ───────────────────────────────────────────────────

    def fly(self, run_name='tmp', convergence_threshold=-1.0, max_loop=3,
            output_mode=False, skip_bad_init_eqdsks=False, run_tx_init=True,
            t_ave_toggle='off', t_ave_window=0.5, t_ave_causal=True, t_ave_ignore_start=0.25):
        r'''! Run TokaMaker-TORAX coupled simulation loop.

        @param convergence_threshold Max fractional change in consumed flux between loops for convergence.
        @param max_loop Maximum number of loops.
        @param run_name Name for this run (used in output directory and log file).
        @param output_mode Output level selector: False (or None), 'minimal', 'normal', or 'debug'.
        @param skip_bad_init_eqdsks If True, skip broken initial gEQDSK files instead of raising.
        @param run_tx_init If True (default), run Loop 0 transport initialization before the main
               coupling loop. If False, skip initialization and start at loop 1.
        @param t_ave_toggle Time-averaging mode: 'off' (no averaging), 'flattop' (average only
               during flat-top), or 'pulse' (average over the whole pulse).
        @param t_ave_window Averaging window size in seconds. Default 0.5 s.
        @param t_ave_causal If True, window is entirely behind the timepoint (backward-looking).
               If False, window is centred on the timepoint.
        @param t_ave_ignore_start Ignore the first N seconds of the pulse when building the
               averaging window (avoids numerical transients). Default 0.25 s.
        '''
        import tempfile

        # Disable JAX's persistent XLA compilation cache before any TORAX/JAX JIT
        # compilation occurs.  Since JAX 0.4.x the cache stores serialized XLA
        # executables keyed by a hash that does NOT include the XLA/JAX version.
        # After a JAX upgrade the old entries are loaded anyway (triggering the
        # "Assume version compatibility. PjRt-IFRT does not track XLA executable
        # versions" warnings) and can produce silently wrong numerical results or
        # semaphore leaks.  Disabling the persistent cache here means JAX still
        # JIT-compiles in memory for the duration of the session (fast after the
        # first call) but never reads or writes stale on-disk entries.
        try:
            import jax
            jax.config.update('jax_enable_compilation_cache', False)
        except Exception:
            pass  # non-fatal: older JAX versions may not have this config key

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

        # Time-averaging settings for sawtooth smoothing
        self._t_ave_toggle = t_ave_toggle
        self._t_ave_window = t_ave_window
        self._t_ave_causal = t_ave_causal
        self._t_ave_ignore_start = t_ave_ignore_start

        dt_str = datetime.now().strftime('%Y-%m-%d_%H%M%S')
        _sim_start_time = time.time()

        self._run_timestamp = None if run_name == 'tmp' else dt_str
        self._output_file_tag = None if run_name == 'tmp' else f'{run_name}_{dt_str}'

        # ── Log file: same directory as toktox_outputs (i.e. cwd / './') ──
        if run_name == 'tmp':
            self._log_file = os.path.abspath('toktox_log_tmp.log')
        else:
            self._log_file = os.path.abspath(f'toktox_log_{run_name}_{dt_str}.log')
        with open(self._log_file, 'w'):
            pass
        print(f'  Log file: {self._log_file}', flush=True)
        self._log(f'Log file: {self._log_file}')

        # In debug mode, attach file handler to Python logging so library
        # messages (TORAX, JAX, etc.) are captured in the log file.
        if self._debug_mode:
            self._logging_configured = False
            self.configure_redirect_to_log()

        # ── Output directory ──
        if self._output_mode is not False:
            if run_name == 'tmp':
                self._out_dir = os.path.join('./toktox_outputs', 'tmp')
                if os.path.exists(self._out_dir):
                    shutil.rmtree(self._out_dir)
            else:
                self._out_dir = os.path.join('./toktox_outputs', self._output_file_tag)
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

        # ── EQDSK directory: persisted only for debug mode ──
        if self._output_mode == 'debug':
            self._eqdsk_dir = self._out_dir
            self._eqdsk_dir_is_temp = False
        else:
            self._eqdsk_dir = tempfile.mkdtemp(prefix='toktox_equil_')
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
        self._print(f'\n{"="*60}\n TokaMaker + TORAX (TokTox) \n run_name = {run_name} | t=[{self._t_init:.1f}, {self._t_final:.1f}] s '
                      f'| {len(self._tm_times)} timepoints | dt={self._tx_dt} s | max_loop={max_loop}')

        err = convergence_threshold + 1.0
        cflux_tx_prev = 0.0
        tm_cflux_psi = []
        tm_cflux_vloop = []
        tx_cflux_psi = []
        tx_cflux_vloop = []

        try:
            # ── Loop 0: Transport initialization (optional) ──
            if run_tx_init:
                self._print(f'\n{"="*60}\n  Loop 0: Transport Initialization\n{"="*60}')
                self._run_tx_init()

            self._current_loop = 1

            # ── Main coupling loop ──
            while err > convergence_threshold and self._current_loop <= max_loop:
                self._print(f'\n{"="*60}\n  Loop {self._current_loop}\n{"="*60}')

                cflux_tx, cflux_tx_vloop = self._run_tx()

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
        if self._output_mode is not False:
            if self._output_mode in ('minimal', 'debug') and self._out_dir is not None:
                profile_evo_name = 'profile_evolution.png'
                lcfs_evo_name = 'lcfs_evolution.png'
                if self._output_file_tag is not None:
                    profile_evo_name = f'{self._output_file_tag}_{profile_evo_name}'
                    lcfs_evo_name = f'{self._output_file_tag}_{lcfs_evo_name}'
                _profile_evo_path = os.path.join(self._out_dir, profile_evo_name)
                _lcfs_evo_path = os.path.join(self._out_dir, lcfs_evo_name)
                try:
                    plot_profile_evolution(self, save_path=_profile_evo_path, display=False)
                except Exception as _e:
                    self._log(f'plot_profile_evolution failed: {_e}')
                try:
                    plot_lcfs_evolution(self, save_path=_lcfs_evo_path, display=False)
                except Exception as _e:
                    self._log(f'plot_lcfs_evolution failed: {_e}')

            if self._output_mode in ('normal', 'minimal', 'debug') and self._out_dir is not None:
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
        n_loops = self._current_loop
        converged = err <= convergence_threshold
        self._print(f'\n{"="*60}')
        if converged:
            self._print(f'  CONVERGED in {n_loops} loops (err={err*100:.3f}%)')
        else:
            self._print(f'  Max loops ({max_loop}) reached (err={err*100:.3f}%)')

        # Print convergence history
        self._print(f'\n  {"Loop":<6} {"cflux TX [Wb]":<16} {"cflux TM [Wb]":<16} {"TX-TM diff %":<14}')
        self._print(f'  {"-"*52}')
        for s in range(len(tx_cflux_psi)):
            diff_pct = np.abs(tx_cflux_psi[s] - tm_cflux_psi[s]) / tm_cflux_psi[s] * 100 if tm_cflux_psi[s] != 0 else np.inf
            self._print(f'  {s+1:<6} {tx_cflux_psi[s]:<16.4f} {tm_cflux_psi[s]:<16.4f} {diff_pct:<14.4f}')
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

    # ─── Visualization ──────────────────────────────────────────────────────────

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


# ── Helpers ──────────────────────────────────────────────────────────────────

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
    for ax in (axes if hasattr(axes, '__iter__') else [axes]): # TODO remove __iter__
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


# ── Profile plot (per-timestep diagnostic) ───────────────────────────────────

def profile_plot(tt, i, t, save_path=None, display=True):
    """Detailed profile comparison at a single timestep."""
    s = tt._state
    psi_N = tt._psi_N

    tm_psi, tm_f_prof, tm_fp_prof, tm_p_prof, tm_pp_prof = tt._tm.get_profiles(npsi=len(tt._psi_N))

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


# ── TokaMaker diagnostic plot ─────────────────────────────────────────────────

def tm_diagnostic_plot(tt, i, t, level_attempts, solve_succeeded, save_path=None, display=True):
    """TokaMaker input/output diagnostic plot for a single timestep."""
    s = tt._state

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


# ── TM loop summary plot ──────────────────────────────────────────────────────

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


# ── Profile evolution plot ────────────────────────────────────────────────────

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

        fig, axes = plt.subplots(2, 4, figsize=(18, 10))
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


# ── Scalar time-series plot ───────────────────────────────────────────────────

def plot_scalars(tt, save_path=None, display=True):
    """Plot 4x3 grid of time-series scalars."""
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
    t_psi_lcfs, y_psi_lcfs = _tx_profile_at_rho(tt, 'psi', 1.0, scale=-1.0/(2.0*np.pi)) # TODO: just pull from state
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

    plt.suptitle('Scalars', fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    _save_or_display(fig, save_path, display)


# ── Coil current plot ─────────────────────────────────────────────────────────

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


# ── LCFS evolution plot ───────────────────────────────────────────────────────

def plot_lcfs_evolution(tt, save_path=None, display=True, one_plot=False):
    """Plot time evolution of the last closed flux surface for each phase.

    Produces phase-split figures by default (rampup, flattop, rampdown),
    or one combined figure when one_plot=True.
    """
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


# ── Movie generation ──────────────────────────────────────────────────────────

def make_movie(tt, save_path=None, display=True, speed_factor=1.0, loop=None, notebook_mode=None):
    """Create pulse movie from simulation data.

    Renders equilibrium plots from stored equilibrium snapshots, generates
    composite frames in a temp directory, encodes to MP4, then cleans up.

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


# ── Interactive widgets (Jupyter) ─────────────────────────────────────────────

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


# ── Physics summary ───────────────────────────────────────────────────────────

def summary(tt):
    """Print and return key physics figures of merit from the simulation."""
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

    out['Ip_max_MA'] = float(np.max(np.abs(s['Ip_tm'])) / 1e6)

    _, y_bN = _tx_scalar(tt, 'beta_N')
    if y_bN is not None:
        out['beta_N_max'] = float(np.nanmax(y_bN))
    out['beta_N_tm_max'] = float(np.nanmax(s['beta_N_tm']))

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
        out['vloop_tm_flattop_avg_V'] = float(np.nanmean(np.array(s['vloop_tm'])[ft_mask]))
        out['vloop_tx_flattop_avg_V'] = float(np.nanmean(np.array(s['vloop_tx'])[ft_mask]))

    print(f"\n{'=' * 55}")
    print("  TokTox Physics Summary")
    print(f"{'=' * 55}")
    for key, val in out.items():
        if val is not None:
            if isinstance(val, float):
                print(f"  {key:30s}  {val:.4g}")
            else:
                print(f"  {key:30s}  {val}")
    print(f"{'=' * 55}\n")
    return out

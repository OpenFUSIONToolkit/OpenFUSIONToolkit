import numpy as np
import pprint
from scipy.interpolate import interp1d, CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import copy
import json
import os
import shutil
from datetime import datetime
import time

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, create_power_flux_fun
from OpenFUSIONToolkit.TokaMaker import toktox_visualization

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
    # 'transport': { # recommended in torax documentation
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
        'apply_outer_patch': True,
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
import logging
import sys
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
import torax


from contextlib import contextmanager
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

    def set_heating(self, generic_heat=None, generic_heat_loc=None, nbi_current=False, ecrh=None, ecrh_loc=None, ohmic=None):
        r'''! Set heating sources for Torax.

        Ohmic heating is always enabled (it is on by default in BASE_CONFIG).
        @param generic_heat Generic heating (dictionary of {time: power_in_watts}).
        @param generic_heat_loc Generic heating deposition location (normalized rho).
        @param ecrh ECRH heating (dictionary of {time: power_in_watts}).
        @param ecrh_loc ECRH deposition location (normalized rho).
        @param nbi_current Whether to include NBI current drive, uses _NBI_W_TO_MA = 1/16e6 to convert heating to current driven.
        '''
        if generic_heat is not None and generic_heat_loc is not None:
            self._generic_heat = generic_heat
            self._generic_heat_loc = generic_heat_loc
        if ecrh is not None and ecrh_loc is not None:
            self._ecrh_heating = ecrh
            self._ecrh_loc = ecrh_loc
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

        if self._generic_heat is not None:
            nbi_times, nbi_pow = zip(*self._generic_heat.items())    
            myconfig.setdefault('sources', {})
            myconfig['sources'].setdefault('generic_heat', {})
            myconfig['sources']['generic_heat']['P_total'] = (nbi_times, nbi_pow)
            myconfig['sources']['generic_heat']['gaussian_location'] = self._generic_heat_loc

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
        self._print(f'  TORAX: running simulation...')
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
                from .toktox_visualization import tm_diagnostic_plot, profile_plot
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
            from .toktox_visualization import tm_loop_summary_plot
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
        import contextlib, os, sys
        @contextlib.contextmanager
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
            output_mode=False, skip_bad_init_eqdsks=False,
            t_ave_toggle='off', t_ave_window=0.5, t_ave_causal=True, t_ave_ignore_start=0.25):
        r'''! Run TokaMaker-TORAX coupled simulation loop.

        @param convergence_threshold Max fractional change in consumed flux between loops for convergence.
        @param max_loop Maximum number of loops.
        @param run_name Name for this run (used in output directory and log file).
        @param output_mode Output level selector: False (or None), 'minimal', 'normal', or 'debug'.
        @param skip_bad_init_eqdsks If True, skip broken initial gEQDSK files instead of raising.
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
            # ── Loop 0: Transport initialization ──
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
                    from .toktox_visualization import plot_scalars
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
                from .toktox_visualization import plot_profile_evolution, plot_lcfs_evolution
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

    # =========================================================================
    #  Visualization wrapper methods (lazy-import from toktox_visualization)
    # =========================================================================

    def make_movie(self, save_path=None, **kwargs):
        r'''! Generate pulse movie from stored psi snapshots.
        @param save_path Path to save MP4 file. If None, does not save.
        '''
        from .toktox_visualization import make_movie
        return make_movie(self, save_path=save_path, **kwargs)

    def plot_scalars(self, save_path=None, display=True, **kwargs):
        r'''! Plot scalar time traces (Ip, Q, Te, ne, power channels, etc.).
        @param save_path Path to save figure. If None, does not save.
        @param display Whether to show the plot.
        '''
        from .toktox_visualization import plot_scalars
        return plot_scalars(self, save_path=save_path, display=display, **kwargs)

    def plot_profiles(self, **kwargs):
        r'''! Interactive profile viewer (ipywidgets slider in Jupyter, static otherwise).'''
        from .toktox_visualization import plot_profiles_interactive
        return plot_profiles_interactive(self, **kwargs)

    def plot_profile_evolution(self, save_path=None, display=True, one_plot=False, **kwargs):
        r'''! Plot profile evolution over time.
        @param save_path Path to save figure. If None, does not save.
        @param display Whether to show the plot.
        @param one_plot If True, combine all pulse phases into one figure.
        '''
        from .toktox_visualization import plot_profile_evolution
        return plot_profile_evolution(self, save_path=save_path, display=display, one_plot=one_plot, **kwargs)


    def plot_coils(self, save_path=None, display=True, **kwargs):
        r'''! Plot coil current traces over the pulse.
        @param save_path Path to save figure. If None, does not save.
        @param display Whether to show the plot.
        '''
        from .toktox_visualization import plot_coils
        return plot_coils(self, save_path=save_path, display=display, **kwargs)

    def plot_lcfs_evolution(self, save_path=None, display=True, one_plot=False, **kwargs):
        r'''! Plot time evolution of the LCFS for each pulse phase (rampup, flattop, rampdown).
        @param save_path Path prefix to save figures. If None, does not save.
        @param display Whether to show the plots.
        @param one_plot If True, combine all pulse phases into one figure.
        '''
        from .toktox_visualization import plot_lcfs_evolution
        return plot_lcfs_evolution(self, save_path=save_path, display=display, one_plot=one_plot, **kwargs)

    def summary(self, **kwargs):
        r'''! Print/display a physics summary of the simulation.'''
        from .toktox_visualization import summary
        return summary(self, **kwargs)


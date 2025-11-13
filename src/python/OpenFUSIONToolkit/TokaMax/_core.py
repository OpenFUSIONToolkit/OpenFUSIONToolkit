import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.interpolate import make_smoothing_spline, interp1d
import torax
import copy
import json
import os
import shutil

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk, create_isoflux

from .baseconfig import BASE_CONFIG

LCFS_WEIGHT = 100.0
N_PSI = 100
_NBI_W_TO_MA = 1/16e6

class MyEncoder(json.JSONEncoder):
    '''! JSON Encoder Object to store simulation results.'''
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

class TokaMax:
    '''! Coupled Grad-Shafranov/Transport Solver Object.'''

    def __init__(self, t_init, t_final, times, g_eqdsk_arr, dt=1, t_res=None, prescribed_currents=False):
        r'''! Initialize the Coupled Grad-Shafranov/Transport Solver Object.
        @param t_init Start time (s).
        @param t_final End time (s).
        @param times Time points of each gEQDSK file.
        @param g_eqdsk_arr Filenames of each gEQDSK file.
        @param dt Time step (s).
        @param t_res Time points to sample output at.
        '''
        self._oftenv = OFT_env(nthreads=2)
        self._gs = TokaMaker(self._oftenv)
        self._state = {}
        self._times = times
        self._boundary = {}
        self._results = {}
        self._init_files = g_eqdsk_arr
        self._init_ip = None
        self._t_init = t_init
        self._t_final = t_final
        self._dt = dt
        self._prescribed_currents = prescribed_currents

        if t_res is None:
            self._t_res = times
        else:
            self._t_res = t_res

        self._state['R'] = np.zeros(len(times))
        self._state['Z'] = np.zeros(len(times))
        self._state['a'] = np.zeros(len(times))
        self._state['kappa'] = np.zeros(len(times))
        self._state['delta'] = np.zeros(len(times))    
        self._state['deltaU'] = np.zeros(len(times))    
        self._state['deltaL'] = np.zeros(len(times))    
        self._state['B0'] = np.zeros(len(times))
        self._state['V0'] = np.zeros(len(times))
        self._state['Ip'] = np.zeros(len(times))
        self._state['pax'] = np.zeros(len(times))
        self._state['beta_pol'] = np.zeros(len(times))
        self._state['vloop'] = np.zeros(len(times))
        self._state['q95'] = np.zeros(len(times))
        self._state['psi_lcfs'] = np.zeros(len(times))
        self._state['psi_axis'] = np.zeros(len(times))

        self._state['ffp_prof'] = {}
        self._state['pp_prof'] = {}
        self._state['ffp_prof_save'] = {}
        self._state['pp_prof_save'] = {}
        self._state['eta_prof'] = {}
        self._state['T_e'] = {}
        self._state['T_i'] = {}
        self._state['n_e'] = {}
        self._state['n_i'] = {}
        self._state['Ptot'] = {}
        # self._state['psi'] = {}

        self._results['lcfs'] = {}
        self._results['dpsi_lcfs_dt'] = {}
        self._results['vloop_tmaker'] = np.zeros([20, len(times)])
        self._results['vloop_torax'] = np.zeros([20, len(times)])
        self._results['q'] = {}
        self._results['jtot'] = {}
        self._results['n_e'] = {}
        self._results['T_e'] = {}
        self._results['T_i'] = {}

        for i, _ in enumerate(times):
            # Calculate geometry
            g = read_eqdsk(g_eqdsk_arr[i])

            self._boundary[i] = g['rzout'].copy()
            self._results['lcfs'][i] = g['rzout'].copy()
            zmax = np.max(self._boundary[i][:,1])
            zmin = np.min(self._boundary[i][:,1])
            rmax = np.max(self._boundary[i][:,0])
            rmin = np.min(self._boundary[i][:,0])
            minor_radius = (rmax - rmin) / 2.0
            rgeo = (rmax + rmin) / 2.0
            highest_pt_idx = np.argmax(self._boundary[i][:,1])
            lowest_pt_idx = np.argmin(self._boundary[i][:,1])
            rupper = self._boundary[i][highest_pt_idx][0]
            rlower = self._boundary[i][lowest_pt_idx][0]
            delta_upper = (rgeo - rupper) / minor_radius
            delta_lower = (rgeo - rlower) / minor_radius

            # Default Scalars
            self._state['R'][i] = g['rcentr']
            self._state['Z'][i] = g['zmid']
            self._state['a'][i] = minor_radius
            self._state['kappa'][i] = (zmax - zmin) / (2.0 * minor_radius)
            self._state['delta'][i] = (delta_upper + delta_lower) / 2.0
            self._state['deltaU'][i] = delta_upper
            self._state['deltaL'][i] = delta_lower
            self._state['B0'][i] = g['bcentr']
            self._state['V0'][i] = g['zaxis']
            self._state['pax'][i] = g['pres'][0]
            self._state['q95'][i] = np.percentile(g['qpsi'], 95)
            self._state['Ip'][i] = abs(g['ip'])
            self._state['psi_axis'][i] = abs(g['psimag'])
            self._state['psi_lcfs'][i] = abs(g['psibry'])

            # Default Profiles
            psi_sample = np.linspace(0.0, 1.0, N_PSI)
            psi_eqdsk = np.linspace(0.0, 1.0, g['nr'])
            ffp_prof = np.interp(psi_sample, psi_eqdsk, g['ffprim'])
            pp_prof = np.interp(psi_sample, psi_eqdsk, g['pprime'])
            self._state['ffp_prof'][i] = {'x': psi_sample, 'y': ffp_prof, 'type': 'linterp'}
            self._state['pp_prof'][i] = {'x': psi_sample, 'y': pp_prof, 'type': 'linterp'}
            # self._state['psi'][i] = np.linspace(g['psimag'], g['psibry'], N_PSI)

            # Normalize profiles
            self._state['ffp_prof'][i]['y'] -= self._state['ffp_prof'][i]['y'][-1]
            self._state['pp_prof'][i]['y'] -= self._state['pp_prof'][i]['y'][-1]
            self._state['ffp_prof'][i]['y'] /= self._state['ffp_prof'][i]['y'][0]
            self._state['pp_prof'][i]['y'] /= self._state['pp_prof'][i]['y'][0]

            self._state['eta_prof'][i]= {
                'x': np.linspace(0.0, 1.0, N_PSI),
                'y': np.zeros(N_PSI),
                'type': 'linterp',
            }
        
        self._nbi_heating = {t_init: 0, t_final: 0}
        self._eccd_heating = {t_init: 0, t_final: 0}
        self._eccd_loc = 0.1
        self._nbi_loc = 0.25
        self._ohmic_power = None

        self._z_eff = None

        self._evolve_density = True
        self._evolve_current = True
        self._evolve_Ti = True
        self._evolve_Te = True

        self._nbar = None
        self._n_e = None
        self._T_i = None
        self._T_e = None

        self._T_i_ped = None
        self._T_e_ped = None
        self._n_e_ped = None

        self._Te_right_bc = None
        self._Ti_right_bc = None
        self._ne_right_bc = None

        self._ohmic = None

        self._gp_s = None
        self._gp_dl = None

        self._chi_min = 0.05
        self._chi_max = 100.0
        self._De_min = 0.05
        self._De_max = 50.0
        self._Ve_min = -10.0
        self._Ve_max = 10.0

        self._targets = None
        self._flux = None
        
    def initialize_gs(self, mesh, weights=None, vsc=None):
        r'''! Initialize GS Solver Object.
        @param mesh Filename of reactor mesh.
        @param vsc Vertical Stability Coil.
        '''
        mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh(mesh)
        self._gs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
        self._gs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
        self._gs.setup(order = 2, F0 = self._state['R'][0]*self._state['B0'][0])

        self._gs.settings.maxits = 500

        if vsc is not None:
            self._gs.set_coil_vsc({vsc: 1.0})
        # self.set_coil_reg(targets, weights=weights, weight_mult=0.1)

    def set_Ip(self, ip):
        r'''! Set plasma current (Amps).
        @param ip Plasma current.
        '''
        self._init_ip = ip
    
    def set_density(self, n_e):
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

    def set_Zeff(self, z_eff):
        r'''! Set plasma effective charge.
        @param z_eff Effective charge.
        '''
        self._z_eff = z_eff
    
    def set_nbar(self, nbar):
        r'''! Set line averaged density over time.
        @param nbar Density (m^-3).
        '''
        self._nbar = nbar
    
    def set_right_bc(self, ne_right_bc=None, Te_right_bc=None, Ti_right_bc=None):
        if ne_right_bc:
            self._ne_right_bc = ne_right_bc
        if Te_right_bc:
            self._Te_right_bc = Te_right_bc
        if Ti_right_bc:
            self._Ti_right_bc = Ti_right_bc

    def set_heating(self, nbi=None, nbi_loc=None, eccd=None, eccd_loc=None):
        r'''! Set heating sources for Torax.
        @param nbi NBI heating (dictionary of heating at times).
        @param nbi_loc Location of NBI heating (rho norm).
        @param eccd ECCD heating (dictionary of heating at times).
        @param eccd_loc Location of ECCD heating (rho norm).
        '''
        if nbi and nbi_loc:
            self._nbi_heating = nbi
            self._nbi_loc = nbi_loc
        if eccd and eccd_loc:
            self._eccd_heating = eccd
            self._eccd_loc = eccd_loc

    def set_pedestal(self, T_i_ped=None, T_e_ped=None, n_e_ped=None):
        r'''! Set pedestals for ion and electron temperatures.
        @pararm T_i_ped Ion temperature pedestal (dictionary of temperature at times).
        @pararm T_e_ped Electron temperature pedestal (dictionary of temperature at times).
        '''
        if T_i_ped:
            self._T_i_ped = T_i_ped
        if T_e_ped:
            self._T_e_ped = T_e_ped
        if n_e_ped:
            self._n_e_ped = n_e_ped

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
    
    def set_Bp(self, Bp):
        for i in range(len(self._times)):
            self._state['beta_pol'][i] = Bp[i]

    def set_Vloop(self, vloop):
        for i in range(len(self._times)):
            self._state['vloop'][i] = vloop[i]
    
    def set_gaspuff(self, s=None, decay_length=None):
        r'''! Set gas puff particle source.
        @param s Particle source (particles/s).
        @param decay_length Decay length from edge (normalized rho coordinates).
        '''
        self._gp_s = s
        self._gp_dl = decay_length
    
    def set_flux(self, flux):
        self._flux = flux
    
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
            
    def set_coil_reg(self, targets=None, i=0, updownsym=False, weights=None, strict_limit=50.0E6, disable_virtual_vsc=True, weight_mult=1.0):
        r'''! Set coil regularization terms.
        @param targets Target values for each coil.
        @param weights Default weight for each coil.
        @param strict_limit Strict limit for coil currents.
        @param disable_virtual_vsc Disable VSC virtual coil. 
        @param weight_mult Factor by which to multiply target weights (reduce to allow for more flexibility).
        '''
        # Set hard constraints
        coil_bounds = {key: [-strict_limit, strict_limit] for key in self._gs.coil_sets}
        # for key in [x for x in self._gs.coil_sets if 'DIV' in x]:   
        #     coil_bounds[key] = [0, 0] # turn off div coils, for now
        self._gs.set_coil_bounds(coil_bounds)

        if self._prescribed_currents and targets:
            self._targets = targets
        
        # Set soft constraints
        regularization_terms = []
        if self._prescribed_currents:
            for name, currents in self._targets.items():
                if name == 'time':
                    continue
                t_current = np.interp(self._times[i], self._targets['time'], currents)
                regularization_terms.append(self._gs.coil_reg_term({name: 1.0},target=t_current,weight=1.0E-3))
        else:
            for name, target_current in targets.items():
                if name == 'time':
                    continue
                regularization_terms.append(self._gs.coil_reg_term({name: 1.0},target=target_current,weight=1.0E-3))

        # Pass regularization terms to TokaMaker
        self._gs.set_coil_reg(reg_terms=regularization_terms)

    def _run_gs(self, step, graph=False):
        r'''! Run the GS solve across n timesteps using TokaMaker.
        @param step Iteration number of the Torax-Tokamaker simulation loop.
        @param graph Whether to display psi graphs at each iteration (for testing).
        @return Consumed flux.
        '''
        for i, t in enumerate(self._times):
            self._gs.set_isoflux(None)
            self._gs.set_flux(None,None)

            Ip_target = abs(self._state['Ip'][i])
            # P0_target = abs(self._state['pax'][i])
            # V0_target = self._state['V0'][i]
            # Ip_ratio = 0.05
            if self._state['beta_pol'][i] != 0:
                Ip_ratio= 1.0/self._state['beta_pol'][i] - 1.0
            self._gs.set_targets(Ip=Ip_target, Ip_ratio = Ip_ratio)
            # self._gs.set_targets(Ip=Ip_target, Ip_ratio=P0_target)

            ffp_prof = self._state['ffp_prof'][i]
            pp_prof = self._state['pp_prof'][i]
            ffp_prof_save = self._state['ffp_prof_save'][i]
            pp_prof_save = self._state['pp_prof_save'][i]

            def mix_profiles(tmp, curr):
                my_prof = {'x': np.zeros(len(curr['x'])), 'y': np.zeros(len(curr['x'])), 'type': 'linterp'}
                for i, x in enumerate(curr['x']):
                    my_prof['x'][i] = x
                    my_prof['y'][i] = 0.1 * tmp['y'][i] + 0.9 * curr['y'][i]

            self._gs.set_profiles(ffp_prof=mix_profiles(ffp_prof_save, ffp_prof), pp_prof=mix_profiles(pp_prof_save, pp_prof))

            self._gs.set_resistivity(eta_prof=self._state['eta_prof'][i])

            lcfs = self._boundary[i]
            isoflux_weights = LCFS_WEIGHT * np.ones(len(lcfs))
            lcfs_psi_target = self._state['psi_lcfs'][i]
            # t_flux = sorted(self._flux.keys())
            # flux_list = [self._flux[t][0] for t in t_flux]
            # lcfs_psi_target = np.interp(t, t_flux, flux_list)

            self._gs.set_flux(lcfs, targets=lcfs_psi_target*np.ones_like(lcfs[:,0]), weights=isoflux_weights)
            # self._gs.set_isoflux(lcfs, weights=1.E3*np.ones_like(lcfs[:,0]))

            # x_points = np.zeros((2,2))
            # x_points[0,:] = lcfs[np.argmin(lcfs[:,1]),:]
            # x_points[1,:] = lcfs[np.argmax(lcfs[:,1]),:]
            # x_weights = 1.E3*np.ones(2)
            # self._gs.set_saddles(x_points, x_weights)

            err_flag = self._gs.init_psi(self._state['R'][i],
                                         self._state['Z'][i],
                                         self._state['a'][i],
                                         self._state['kappa'][i],
                                         self._state['delta'][i])
            
            if err_flag:
                print("Error initializing psi.")

            if graph:
                fig, ax = plt.subplots(1,1)
                self._gs.plot_machine(fig,ax,coil_colormap='seismic',coil_symmap=True,coil_scale=1.E-6,coil_clabel=r'$I_C$ [MA]')
                self._gs.plot_psi(fig,ax,xpoint_color='r',vacuum_nlevels=4)
                ax.plot(self._boundary[i][:, 0], self._boundary[i][:, 1], color='r')
                # ax.set_title(f'i={i}')
                ax.set_title(f't={t}')
                plt.savefig(f'nsf/t={t}.png')
                # plt.show()

            self._gs.update_settings()
            err_flag = self._gs.solve()

            self._gs.print_info()
            self._gs_update(i)

            # if i:
            #     # Compute loop voltage
            #     dt = self._times[i] - self._times[i-1]
            #     dpsi_lcfs_dt = (self._state['psi_lcfs'][i] - self._state['psi_lcfs'][i-1]) / dt
            #     self._results['dpsi_lcfs_dt'][i] = dpsi_lcfs_dt
            self._gs.save_eqdsk('tmp/{:03}.{:03}.eqdsk'.format(step, i),lcfs_pad=0.001,run_info='TokaMaker EQDSK', cocos=2)

            if self._prescribed_currents:
                if i < len(self._times) - 1:
                    self.set_coil_reg(i=i+1)
            else:
                coil_targets, _ = self._gs.get_coil_currents()
                self.set_coil_reg(targets=coil_targets)

        # consumed_flux = self._state['psi_lcfs'][-1] - self._state['psi_lcfs'][0]
        consumed_flux = np.trapezoid(self._times, self._state['vloop'])
        return consumed_flux
        
    def _gs_update(self, i):
        r'''! Update internal state and coil current results based on results of GS solver.
        @param i Timestep of the solve.
        '''
        eq_stats = self._gs.get_stats()
        self._state['Ip'][i] = eq_stats['Ip']

        self._state['psi_lcfs'][i] = self._gs.psi_bounds[0]
        self._state['psi_axis'][i] = self._gs.psi_bounds[1]

        if 'psi_lcfs_tmaker' not in self._results:
            self._results['psi_lcfs_tmaker'] = {'x': np.zeros(len(self._times)), 'y': np.zeros(len(self._times))}
        self._results['psi_lcfs_tmaker']['x'][i] = self._times[i]
        self._results['psi_lcfs_tmaker']['y'][i] = self._state['psi_lcfs'][i]

        self._state['vloop'][i] = self._gs.calc_loopvoltage()
        
        # Update Results
        coils, _ = self._gs.get_coil_currents()
        if i == 0:
            self._results['COIL'] = {coil: {} for coil in coils}
        for coil, current in coils.items():
            self._results['COIL'][coil][self._times[i]] = current * 1.0 # TODO: handle nturns > 1

    def _get_torax_config(self, step):
        r'''! Generate config object for Torax simulation. Modifies BASE_CONFIG based on current simulation state.
        @param step Iteration number of the Torax-Tokamaker simulation loop.
        @return Torax config object.
        '''
        myconfig = copy.deepcopy(BASE_CONFIG)

        myconfig['numerics'] = {
            't_initial': self._t_init,
            't_final': self._t_final,  # length of simulation time in seconds
            'fixed_dt': self._dt, # fixed timestep
            'evolve_ion_heat': self._evolve_Ti, # solve ion heat equation
            'evolve_electron_heat': self._evolve_Te, # solve electron heat equation
            'evolve_current': self._evolve_current, # solve current equation
            'evolve_density': self._evolve_density, # solve density equation
        }

        myconfig['geometry'] = {
            'geometry_type': 'eqdsk',
            'geometry_directory': '/Users/johnl/Desktop/discharge-model', 
            'last_surface_factor': 0.95,  # TODO: tweak
            # 'n_surfaces': 10,
            'Ip_from_parameters': True,
            'geometry_configs': {
                t: {'geometry_file': self._init_files[i]} for i, t in enumerate(self._times)
            },
        }
        if step:
            myconfig['geometry']['geometry_configs'] = {
                t: {'geometry_file': 'tmp/{:03}.{:03}.eqdsk'.format(step, i)} for i, t in enumerate(self._times)
            }

        myconfig['profile_conditions']['Ip'] = {
            t: abs(self._state['Ip'][i]) for i, t in enumerate(self._times)
        }
        myconfig['profile_conditions']['psi'] = {
            t: {0.0: self._state['psi_axis'][i], 1.0: self._state['psi_lcfs'][i]} for i, t in enumerate(self._times)
        }
        if self._flux:
            myconfig['profile_conditions']['psi'] = self._flux

        if self._init_ip:
             myconfig['profile_conditions']['Ip'] = self._init_ip
        
        if self._n_e:
            myconfig['profile_conditions']['n_e'] = self._n_e
        
        if self._T_e:
            myconfig['profile_conditions']['T_e'] = self._T_e
        
        if self._T_i:
            myconfig['profile_conditions']['T_i'] = self._T_i
        
        if self._z_eff:
            myconfig['plasma_composition']['Z_eff'] = self._z_eff
        
        myconfig['sources']['ecrh']['P_total'] = self._eccd_heating
        myconfig['sources']['ecrh']['gaussian_location'] = self._eccd_loc

        if self._ohmic_power:
            myconfig['sources']['ohmic']['mode'] = 'PRESCRIBED'
            myconfig['sources']['ohmic']['prescribed_values'] = self._ohmic_power

        nbi_times, nbi_pow = zip(*self._nbi_heating.items())
        myconfig['sources']['generic_heat']['P_total'] = (nbi_times, nbi_pow)
        myconfig['sources']['generic_heat']['gaussian_location'] = self._nbi_loc
        myconfig['sources']['generic_current']['I_generic'] = (nbi_times, _NBI_W_TO_MA * np.array(nbi_pow))
        myconfig['sources']['generic_current']['gaussian_location'] = self._nbi_loc

        # if self._flux:
        #     myconfig['profile_conditions']['psi'] = self._flux

        if self._T_i_ped:
            myconfig['pedestal']['T_i_ped'] = self._T_i_ped
        if self._T_e_ped:
            myconfig['pedestal']['T_e_ped'] = self._T_e_ped
        
        if self._n_e_ped:
            myconfig['pedestal']['n_e_ped_is_fGW'] = False
            myconfig['pedestal']['n_e_ped'] = self._n_e_ped
        
        if self._nbar:
            myconfig['profile_conditions']['nbar'] = self._nbar

        if self._ne_right_bc:
            myconfig['profile_conditions']['n_e_right_bc_is_fGW'] = False
            myconfig['profile_conditions']['n_e_right_bc'] = self._ne_right_bc

        if self._Te_right_bc:
            myconfig['profile_conditions']['T_e_right_bc'] = self._Te_right_bc
        if self._Ti_right_bc:
            myconfig['profile_conditions']['T_i_right_bc'] = self._Ti_right_bc
        
        # if self._ohmic:
        #     myconfig['sources']['ohmic'] = {
        #         'mode': 'PRESCRIBED',
        #         'prescribed_values': self._ohmic,
        #     }
        
        if self._gp_s and self._gp_dl:
            myconfig['sources']['gas_puff'] = {
                'S_total': self._gp_s,
                'puff_decay_length': self._gp_dl,
            }

        myconfig['transport']['chi_min'] = self._chi_min
        myconfig['transport']['chi_max'] = self._chi_max
 
        myconfig['transport']['D_e_min'] = self._De_min
        myconfig['transport']['D_e_max'] = self._De_max

        myconfig['transport']['V_e_min'] = self._Ve_min
        myconfig['transport']['V_e_max'] = self._Ve_max

        # print(myconfig)
        with open('torax_config.json', 'w') as json_file:
            json.dump(myconfig, json_file, indent=4, cls=MyEncoder)
        torax_config = torax.ToraxConfig.from_dict(myconfig)
        return torax_config

    def _run_transport(self, step, graph=False):
        r'''! Run the Torax simulation.
        @param step Iteration number of the Torax-Tokamaker simulation loop.
        @param graph Whether to display profiles at each iteration (for testing).
        @return Consumed flux.
        '''
        myconfig = self._get_torax_config(step)
        data_tree, hist = torax.run_simulation(myconfig, log_timestep_info=False)

        if hist.sim_error != torax.SimError.NO_ERROR:
            print(hist.sim_error)
            raise ValueError(f'TORAX failed to run the simulation.')
        
        v_loops = np.zeros(len(self._times))
        for i, t in enumerate(self._times):
            self._transport_update(step, i, data_tree)
            v_loops[i] = data_tree.scalars.v_loop_lcfs.sel(time=t, method='nearest') / (2.0 * np.pi)
            self._state['vloop'][i] = v_loops[i]
        
        self._res_update(data_tree)

        consumed_flux = self._state['psi_lcfs'][-1] - self._state['psi_lcfs'][0]
        return consumed_flux
    
    def _transport_update(self, step, i, data_tree, smooth=False):
        r'''! Update the simulation state and simulation results based on results of the Torax simulation.
        @param i Timestep of the solve.
        @param data_tree Result object from Torax.
        @smooth Whether to smooth profiles generated by Torax.
        '''
        t = self._times[i]
        
        self._state['Ip'][i] = data_tree.scalars.Ip.sel(time=t, method='nearest')
        self._state['pax'][i] = data_tree.profiles.pressure_thermal_total.sel(time=t, rho_norm=0.0, method='nearest')
        self._state['beta_pol'][i] = data_tree.scalars.beta_pol.sel(time=t, method='nearest')
        self._state['q95'][i] = data_tree.scalars.q95.sel(time=t, method='nearest')

        self._state['ffp_prof_save'][i] = self._state['ffp_prof'][i]
        self._state['pp_prof_save'][i] = self._state['pp_prof'][i]

        ffprime = data_tree.profiles.FFprime.sel(time=t, method='nearest')
        pprime = data_tree.profiles.pprime.sel(time=t, method='nearest')

        psi_norm = data_tree.profiles.psi_norm.sel(time=t, method='nearest')

        self._state['ffp_prof'][i] = {
            'x': [psi_norm.sel(rho_face_norm=rfn, method='nearest') for rfn in ffprime.coords['rho_face_norm'].values],
            'y': ffprime.to_numpy(),
            'type': 'linterp',
        }

        psi_sample = np.linspace(0.0, 1.0, N_PSI)
        ffp_sample = np.interp(psi_sample, self._state['ffp_prof'][i]['x'], self._state['ffp_prof'][i]['y'])
        self._state['ffp_prof'][i]['x'] = psi_sample
        self._state['ffp_prof'][i]['y'] = ffp_sample

        self._state['pp_prof'][i] = {
            'x': [psi_norm.sel(rho_face_norm=rfn, method='nearest') for rfn in pprime.coords['rho_face_norm'].values],
            'y': pprime.to_numpy(),
            'type': 'linterp',
        }

        pp_sample = np.interp(psi_sample, self._state['pp_prof'][i]['x'], self._state['pp_prof'][i]['y'])
        self._state['pp_prof'][i]['x'] = psi_sample
        self._state['pp_prof'][i]['y'] = pp_sample

        # Normalize profiles
        self._state['ffp_prof'][i]['y'] -= self._state['ffp_prof'][i]['y'][-1]
        self._state['pp_prof'][i]['y'] -= self._state['pp_prof'][i]['y'][-1]
        self._state['ffp_prof'][i]['y'] /= self._state['ffp_prof'][i]['y'][0]
        self._state['pp_prof'][i]['y'] /= self._state['pp_prof'][i]['y'][0]

        # Smooth Profiles
        def make_smooth(x, y):
            spline = make_smoothing_spline(x, y, lam=0.1)
            smoothed = spline(x)
            return smoothed

        if smooth:
            self._state['ffp_prof'][i]['y'] = make_smooth(self._state['ffp_prof'][i]['x'], self._state['ffp_prof'][i]['y'])
            self._state['pp_prof'][i]['y'] = make_smooth(self._state['pp_prof'][i]['x'], self._state['pp_prof'][i]['y'])

        t_i = data_tree.profiles.T_i.sel(time=t, method='nearest')
        t_e = data_tree.profiles.T_e.sel(time=t, method='nearest')
        n_i = data_tree.profiles.n_i.sel(time=t, method='nearest')
        n_e = data_tree.profiles.n_e.sel(time=t, method='nearest')

        self._state['T_i'][i] = {
            'x': np.pow(t_i.coords['rho_norm'].values, 2),
            'y': t_i.to_numpy(),
        }
        self._state['T_e'][i] = {
            'x': np.pow(t_e.coords['rho_norm'].values, 2),
            'y': t_e.to_numpy(),
        }
        self._state['n_i'][i] = {
            'x': np.pow(n_i.coords['rho_norm'].values, 2),
            'y': n_i.to_numpy(),
        }
        self._state['n_e'][i] = {
            'x': np.pow(n_e.coords['rho_norm'].values, 2),
            'y': n_e.to_numpy(),
        }

        # def spitzer_resistivity(n,T):
        #     def log_lambda(n,T):
        #         return 24.0-np.log(np.sqrt(n/1.E6)/T)
        #     return 5.253E-5*log_lambda(n,T)/np.power(T,1.5)
        # self._state['eta_prof'][i] = spitzer_resistivity(n_e.to_numpy(), t_e.to_numpy())

        ptot = data_tree.profiles.pressure_thermal_total.sel(time=t, method='nearest')
        self._state['Ptot'][i] = {
            'x': np.pow(ptot.coords['rho_norm'].values, 2),
            'y': ptot.to_numpy(),
        }

        conductivity = data_tree.profiles.sigma_parallel.sel(time=t, method='nearest')
        self._state['eta_prof'][i] = {
            'x': np.pow(conductivity.coords['rho_norm'].values, 2),
            'y': 1.0 / conductivity.to_numpy(),
            'type': 'linterp',
        }
        psi_sample = np.linspace(0.0, 1.0, N_PSI)
        eta_sample = np.interp(psi_sample, self._state['eta_prof'][i]['x'], self._state['eta_prof'][i]['y'])
        self._state['eta_prof'][i]['x'] = psi_sample
        self._state['eta_prof'][i]['y'] = eta_sample

        self._state['psi_lcfs'][i] = data_tree.profiles.psi.sel(time=t, rho_norm=1.0, method='nearest') / (2.0 * np.pi)
        self._state['psi_axis'][i] = data_tree.profiles.psi.sel(time=t, rho_norm=0.0, method='nearest') / (2.0 * np.pi)

    def _res_update(self, data_tree):

        self._results['t_res'] = self._t_res
        self._results['eqtimes'] = self._times

        for t in self._t_res:
            self._results['T_e'][t] = {
                'x': np.pow(list(data_tree.profiles.T_e.coords['rho_norm'].values), 2),
                'y': data_tree.profiles.T_e.sel(time=t, method='nearest').to_numpy()
            }
            self._results['T_i'][t] = {
                'x': np.pow(list(data_tree.profiles.T_i.coords['rho_norm'].values), 2),
                'y': data_tree.profiles.T_i.sel(time=t, method='nearest').to_numpy()
            }
            self._results['n_e'][t] = {
                'x': np.pow(list(data_tree.profiles.n_e.coords['rho_norm'].values), 2),
                'y': data_tree.profiles.n_e.sel(time=t, method='nearest').to_numpy()
            }

        self._results['E_fusion'] = {
            'x': list(data_tree.scalars.E_fusion.coords['time'].values),
            'y': data_tree.scalars.E_fusion.to_numpy()
        }

        self._results['Q'] = {
            'x': list(data_tree.scalars.Q_fusion.coords['time'].values),
            'y': data_tree.scalars.Q_fusion.to_numpy(),
        }

        self._results['Ip'] = {
            'x': list(data_tree.scalars.Ip.coords['time'].values),
            'y': data_tree.scalars.Ip.to_numpy(),
        }

        self._results['B0'] = {
            'x': list(data_tree.scalars.B_0.coords['time'].values),
            'y': data_tree.scalars.B_0.to_numpy(),
        }

        # self._results['q'][i] = {
        #     'x': list(data_tree.profiles.q.coords['rho_face_norm'].values),
        #     'y': data_tree.profiles.q.sel(time=t, method='nearest').to_numpy()
        # }

        # self._results['jtot'][i] = {
        #     'x': list(data_tree.profiles.j_total.coords['rho_norm'].values),
        #     'y': data_tree.profiles.j_total.sel(time=t, method='nearest').to_numpy()
        # }

        self._results['n_e_line_avg'] = {
            'x': list(data_tree.scalars.n_e_line_avg.coords['time'].values),
            'y': data_tree.scalars.n_e_line_avg.to_numpy(),
        }

        self._results['n_i_line_avg'] = {
            'x': list(data_tree.scalars.n_i_line_avg.coords['time'].values),
            'y': data_tree.scalars.n_i_line_avg.to_numpy(),
        }

        my_times = list(data_tree.profiles.T_e.coords['time'].values)
        T_e_line_avg = np.array([])
        for my_t in my_times:
            T_e_line_avg = np.append(T_e_line_avg, data_tree.profiles.T_e.sel(time=my_t).mean(dim='rho_norm'))
        self._results['T_e_line_avg'] = {
            'x': my_times,
            'y': T_e_line_avg,
        }

        my_times = list(data_tree.profiles.T_i.coords['time'].values)
        T_i_line_avg = np.array([])
        for my_t in my_times:
            T_i_line_avg = np.append(T_i_line_avg, data_tree.profiles.T_i.sel(time=my_t).mean(dim='rho_norm'))
        self._results['T_i_line_avg'] = {
            'x': my_times,
            'y': T_i_line_avg,
        }
        
        n_e_core = data_tree.profiles.n_e.sel(rho_norm=0.0)
        self._results['n_e_core'] = {
            'x': list(n_e_core.coords['time'].values),
            'y': n_e_core.to_numpy(),
        }

        n_i_core = data_tree.profiles.n_i.sel(rho_norm=0.0)
        self._results['n_i_core'] = {
            'x': list(n_i_core.coords['time'].values),
            'y': n_i_core.to_numpy(),
        }

        T_e_core = data_tree.profiles.T_e.sel(rho_norm=0.0)
        self._results['T_e_core'] = {
            'x': list(T_e_core.coords['time'].values),
            'y': T_e_core.to_numpy(),
        }

        T_i_core = data_tree.profiles.T_i.sel(rho_norm=0.0)
        self._results['T_i_core'] = {
            'x': list(T_i_core.coords['time'].values),
            'y': T_i_core.to_numpy(),
        }

        self._results['beta_N'] = {
            'x': list(data_tree.scalars.beta_N.coords['time'].values),
            'y': data_tree.scalars.beta_N.to_numpy(),
        }

        self._results['q95'] = {
            'x': list(data_tree.scalars.q95.coords['time'].values),
            'y': data_tree.scalars.q95.to_numpy(),
        }

        self._results['H98'] = {
            'x': list(data_tree.scalars.H98.coords['time'].values),
            'y': data_tree.scalars.H98.to_numpy(),
        }

        self._results['v_loop_lcfs'] = {
            'x': list(data_tree.scalars.v_loop_lcfs.coords['time'].values),
            'y': data_tree.scalars.v_loop_lcfs.to_numpy(),
        }

        psi_lcfs = data_tree.profiles.psi.sel(rho_norm = 1.0)
        self._results['psi_lcfs_torax'] = {
            'x': list(psi_lcfs.coords['time'].values),
            'y': psi_lcfs.to_numpy(),
        }

        self._results['li3'] = {
            'x': list(data_tree.scalars.li3.coords['time'].values),
            'y': data_tree.scalars.li3.to_numpy(),
        }

        self._results['P_alpha_total'] = {
            'x': list(data_tree.scalars.P_alpha_total.coords['time'].values),
            'y': data_tree.scalars.P_alpha_total.to_numpy(),
        }

        self._results['P_aux_total'] = {
            'x': list(data_tree.scalars.P_aux_total.coords['time'].values),
            'y': data_tree.scalars.P_aux_total.to_numpy(),
        }

        self._results['P_ohmic_e'] = {
            'x': list(data_tree.scalars.P_ohmic_e.coords['time'].values),
            'y': data_tree.scalars.P_ohmic_e.to_numpy(),
        }

        self._results['P_radiation_e'] = {
            'x': list(data_tree.scalars.P_radiation_e.coords['time'].values),
            'y': -1.0 * data_tree.scalars.P_radiation_e.to_numpy(),
        }

        self._results['P_SOL_total'] = {
            'x': list(data_tree.scalars.P_SOL_total.coords['time'].values),
            'y': -1.0 * data_tree.scalars.P_SOL_total.to_numpy(),
        }

    def save_state(self, fname):
        r'''! Save intermediate simulation state to JSON.
        @param fname Filename to save to.
        '''
        with open(fname, 'w') as f:
            json.dump(self._state, f, cls=MyEncoder)
    
    def save_res(self):
        r'''! Save simulation results to JSON.'''
        with open(self._fname_out, 'w') as f:
            json.dump(self._results, f, cls=MyEncoder)

    def fly(self, convergence_threshold=-1.0, save_states=False, graph=False, max_step=50, out='res.json', remove_tmp=False):
        r'''! Run Tokamaker-Torax simulation loop until convergence or max_step reached. Saves results to JSON object.
        @pararm convergence_threshold Maximum percent difference between iterations allowed for convergence.
        @param save_states Save intermediate simulation states (for testing).
        @param graph Whether to display psi and profile graphs at each iteration (for testing).
        @param max_step Maximum number of simulation iterations allowed.
        '''
        if not remove_tmp:
            del_tmp = input('Delete temporary storage? [y/n] ')
            if del_tmp != 'y':
                quit()
        with open('convergence_history.txt', 'w'):
            pass
        shutil.rmtree('./tmp')
        os.mkdir('./tmp')

        self._fname_out = out

        err = convergence_threshold + 1.0
        step = 0
        cflux_prev = 0.0

        while err > convergence_threshold and step < max_step:
            cflux = self._run_transport(step, graph=graph)
            if save_states:
                self.save_state('tmp/ts_state{}.json'.format(step))
            step += 1

            cflux_gs = self._run_gs(step, graph=graph)
            if save_states:
                self.save_state('tmp/gs_state{}.json'.format(step-1))

            self.save_res()

            with open('convergence_history.txt', 'a') as f:
                print("GS CF = {}".format(cflux_gs), file=f)
                print("TS CF = {}".format(cflux), file=f)

            err = np.abs(cflux - cflux_prev) / cflux_prev
            cflux_prev = cflux
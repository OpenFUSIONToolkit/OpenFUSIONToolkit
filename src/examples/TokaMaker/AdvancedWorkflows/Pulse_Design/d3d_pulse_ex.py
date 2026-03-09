from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import read_eqdsk

import os
import json

import numpy as np
import matplotlib.pyplot as plt

# To access DIII-D eqdsk files, please contact Columbia University Fusion Research Center
eqdsk_names = sorted(os.listdir('163303/eqs'))
eqdsks = []
eqtimes = []

ffp = np.array([])
pp = np.array([])

for fname in eqdsk_names:
    if 'OMFIT' in fname or 'DS_Store' in fname:
        continue
    tag, _ = fname.split('.')
    _, t = tag.split('-')
    t = float(t) / 1e3
    eqtimes.append(t)
    eqdsks.append(f'163303/eqs/{fname}')

    # Create and set profiles
    g = read_eqdsk(eqdsks[-1])
    ffprim = g['ffprim']
    pprime = g['pprime']
    ffp = np.append(ffp, ffprim)
    pp = np.append(pp, pprime)

def read_pfile(path):
    data = {}
    key = ''
    with open(path) as f:
        for line in f:
            if '3 N Z A' in line:
                break
            if line.startswith('201'):
                key = line.split()[2]
                data[key] = np.array([])
            else:
                psi, dat, _ = line.split()
                psi = float(psi)
                dat = float(dat)
                data[key] = np.append(data[key], dat)
    return data

# To access DIII-D profile data, please contact Columbia University Fusion Research Center
prof_names = sorted(os.listdir('163303/profs'))
prof_times = []
ne = []
Te = []
Ti = []
for fname in prof_names:
    if 'OMFIT' in fname or 'DS_Store' in fname:
        continue
    _, t = fname.split('.')
    t = float(t) / 1e3
    prof_times.append(t)

    path = f'163303/profs/{fname}'
    data = read_pfile(path)

    Te.append(data['te(KeV)'] * 1e3)
    Ti.append(data['ti(KeV)'] * 1e3)
    ne.append(data['ne(10^20/m^3)'] * 1e20)

# To access DIII-D current/bp data, please contact Columbia University Fusion Research Center
ip_f = open('163303/ip.json')
ip_json = json.load(ip_f)
ip_t = np.array(ip_json['time']) / 1e3
ip = ip_json['data']
bpol_f = open('163303/betap.json')
bpol_json = json.load(bpol_f)
bpol_t = np.array(bpol_json['time']) / 1e3
bpol = bpol_json['data']

myOFT = OFT_env(nthreads=2)
mygs = TokaMaker(myOFT)
mygs.settings.maxits = 200

mesh_pts,mesh_lc,mesh_reg,coil_dict,cond_dict = load_gs_mesh('163303/DIIID_mesh.h5')
mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
mygs.setup_regions(cond_dict=cond_dict,coil_dict=coil_dict)
g1 = read_eqdsk(eqdsks[0])
mygs.setup(order=2,F0=g1['rcentr']*g1['bcentr'])

coil_bounds = {key: [-5.0E6, 5.0E6] for key in mygs.coil_sets}
mygs.set_coil_bounds(coil_bounds)

target_currents = {
    'ECOILA': 0.0,
    'ECOILB': 0.0,
    'F1A': 0.0,
    'F2A': 0.0,
    'F3A': 0.0,
    'F4A': 0.0,
    'F5A': 0.0,
    'F6A': 0.0,
    'F7A': 0.0,
    'F8A': 0.0,
    'F9A': 0.0,
    'F1B': 0.0,
    'F2B': 0.0,
    'F3B': 0.0,
    'F4B': 0.0,
    'F5B': 0.0,
    'F6B': 0.0,
    'F7B': 0.0,
    'F8B': 0.0,
    'F9B': 0.0,
}
 
# Set regularization weights
def set_coil_reg(gs, target_currents):
    regularization_terms = []
    for name, target_current in target_currents.items():
        # Set specific target currents from input equilibrium and different weights depending on the coil set
        if name.startswith('ECOIL'):
            regularization_terms.append(gs.coil_reg_term({name: 1.0},target=target_current,weight=61.0))
        elif name.startswith('F'):
            if name.startswith('F5'):
                regularization_terms.append(gs.coil_reg_term({name: 1.0},target=target_current,weight=1.E2))
            else:
                regularization_terms.append(gs.coil_reg_term({name: 1.0},target=target_current,weight=1.E0))
    # Set zero target current and small weight on virtual VSC to allow up-down adjustment
    regularization_terms.append(gs.coil_reg_term({'#VSC': 1.0},target=0.0,weight=1.E-2))
    # Pass regularization terms to TokaMaker
    gs.set_coil_reg(reg_terms=regularization_terms)

set_coil_reg(mygs, target_currents)

# Helper functions
def spitzer_resistivity(n,T): # NOTE: requires T in eV and n in m^-3
    def log_lambda(n,T):
        return 24.0-np.log(np.sqrt(n/1.E6)/T)
    return 5.253E-5*log_lambda(n,T)/np.power(T,1.5)

def coil_vec2dict(self,currents):
    current_array = np.zeros((self.ncoils,), dtype=np.float64)
    for coil_key, coil_current in currents.items():
        current_array[self.coil_sets[coil_key]['id']] = coil_current
    return current_array

# Calculate plasma shape parameters
R = []
Z = []
a = []
kappa = []
delta = []
lcfs = []

# To access DIII-D data, please contact Columbia University Fusion Research Center
for i, eqdsk in enumerate(eqdsks):
    t = eqtimes[i]
    g = read_eqdsk(eqdsks[i])
    mylcfs = g['rzout']
    zmax = np.max(mylcfs[:,1])
    zmin = np.min(mylcfs[:,1])
    rmax = np.max(mylcfs[:,0])
    rmin = np.min(mylcfs[:,0])
    minor_radius = (rmax - rmin) / 2.0
    rgeo = (rmax + rmin) / 2.0
    highest_pt_idx = np.argmax(mylcfs[:,1])
    lowest_pt_idx = np.argmin(mylcfs[:,1])
    rupper = mylcfs[highest_pt_idx][0]
    rlower = mylcfs[lowest_pt_idx][0]
    delta_upper = (rgeo - rupper) / minor_radius
    delta_lower = (rgeo - rlower) / minor_radius
    mykappa = (zmax - zmin) / (2.0 * minor_radius)
    mydelta = (delta_upper + delta_lower) / 2.0

    R.append(g['rcentr'])
    Z.append(g['zmid'])
    a.append(minor_radius)
    kappa.append(mykappa)
    delta.append(mydelta)

    lcfs.append(g['rzout'])

frames = []
coil_hist = []
flux_hist = []
volt_hist = []
Lcoils = None
prev_coils = None

psi_target = 0.0
psi = 0.0
loop_voltage = 0.0
dt = 0.1
timesteps = np.arange(0.1, 5.0, dt)
graph_shape = False
for i, t in enumerate(timesteps):
    def interp_prof(t, times, profs):
        for i in range(1, len(times)):
            if t <= times[i-1]:
                return profs[i-1]
            if t > times[i-1] and t < times[i]:
                dt = times[i] - times[i-1]
                alpha = (times[i] - t) / dt
                return (1.0 - alpha) * profs[i-1] + alpha * profs[i]
            
    ffp_prof = interp_prof(t, eqtimes, ffp)
    pp_prof = interp_prof(t, eqtimes, pp)
    psi_eqdsk = np.linspace(0.0,1.0,np.size(ffprim))
    psi_sample = np.linspace(0.0, 1.0, 50)

    psi_prof = np.copy(psi_sample)
    ffp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_eqdsk,ffprim)))).copy()
    pp_prof = np.transpose(np.vstack((psi_prof,np.interp(psi_sample,psi_eqdsk,pprime)))).copy()

    mygs.set_profiles(ffp_prof={'type': 'linterp', 'y': ffp_prof[:,1], 'x': psi_sample},pp_prof={'type': 'linterp', 'y': pp_prof[:,1], 'x': psi_sample})

    # Set current
    Ip = np.interp(t, ip_t, ip)
    bp = np.interp(t, bpol_t, bpol)
    mygs.set_targets(Ip=Ip,Ip_ratio=(1.0/bp - 1.0))

    # Set plasma shape
    R0 = np.interp(t, eqtimes, R)
    Z0 = np.interp(t, eqtimes, Z)
    a0 = np.interp(t, eqtimes, a)
    kappa0 = np.interp(t, eqtimes, kappa)
    delta0 = np.interp(t, eqtimes, delta)

    # lcfs = create_isoflux(20, R0, Z0, a0, kappa0, delta0)
    mylcfs = interp_prof(t, eqtimes, lcfs)
    if i == 0:
        mygs.set_isoflux(mylcfs, 1.0E3 * np.ones_like(mylcfs[:, 0]))    # Target shape
    if i > 0:
        delta_psi = (loop_voltage * dt) / (2 * np.pi)
        psi_target -= delta_psi
        mygs.set_flux(mylcfs, psi_target*np.ones_like(mylcfs[:,0]), weights=1.0E3 * np.ones_like(mylcfs[:, 0])) # Target specific flux value
    # TODO: Set saddles?

    mygs.init_psi(R0, Z0, a0, kappa0, delta0)
    if i > 0:
        mygs.set_psi_dt(psi0,dt)
    mygs.solve()

    psi0 = mygs.get_psi(False)

    if graph_shape:
        xpts, _ = mygs.get_xpoints()
        fig, ax = plt.subplots(1,1)
        mygs.plot_machine(fig,ax,coil_colormap='seismic',coil_symmap=True,coil_scale=1.0E-3,coil_clabel=r'$I_{coil}$ [kA]')
        mygs.plot_psi(fig,ax,xpoint_color='k',vacuum_nlevels=6,plasma_nlevels=6)
        mygs.plot_constraints(fig,ax,isoflux_color='tab:red',isoflux_marker='.')
        ax.set_ylabel('Z [m]')
        _ =ax.set_xlabel('R [m]')
        ax.scatter(xpts[:, 0], xpts[:, 1], color='lime')
        plt.show()

    # Save coil currents and boundary flux
    coil_currents, _ = mygs.get_coil_currents()
    if i > 0:
        volt_hist.append(np.dot(Lcoils[:-1,:-1],coil_vec2dict(mygs,coil_currents)-coil_vec2dict(mygs,prev_coils))/dt)

    prev_coils = coil_currents.copy()
    coil_hist.append(coil_vec2dict(mygs,coil_currents))
    flux_hist.append([psi_target,mygs.psi_bounds[0]])

    # Update coil reg terms, loop voltage, and psi target
    if i == 0:
        psi_target = mygs.psi_bounds[0]
        Lcoils = mygs.get_coil_Lmat()

    density = interp_prof(t, prof_times, ne)
    temp = interp_prof(t, prof_times, Te)
    eta_prof = spitzer_resistivity(density, temp)
    psi_prof = np.linspace(0.0, 1.0, len(eta_prof))
    mygs.set_resistivity({'type': 'linterp', 'x': psi_prof, 'y': eta_prof})
    loop_voltage = mygs.calc_loopvoltage()

    set_coil_reg(mygs, coil_currents)

flux_hist = np.array(flux_hist)
fig, ax = plt.subplots(1,1)
ax.plot(timesteps,flux_hist[:,0],'k',label=r'Target')
ax.plot(timesteps,flux_hist[:,1],'--',color='tab:red',label=r'Actual')
ax.grid(True)
ax.set_xlim(left=0.0)
ax.set_xlabel('Time [ms]')
ax.set_ylabel(r'$\psi$ [mWb/rad]')
_ = ax.legend()
plt.show()

coil_hist = np.array(coil_hist)
fig, ax = plt.subplots(2,1,sharex=True)
for name, coil in mygs.coil_sets.items():
    if name.startswith('E'):
        ax[0].plot(timesteps,coil_hist[:,coil['id']],label=name)
    if name.startswith('F'):
        ax[1].plot(timesteps,coil_hist[:,coil['id']],label=name)
for ax_tmp in ax:
    ax_tmp.grid(True)
    ax_tmp.set_xlim(left=0.0)
    ax_tmp.set_ylim(-2.5E6, 2.5E6)
    ax_tmp.set_ylabel('Coil Current [A]')
    ax_tmp.legend()
_ = ax[-1].set_xlabel('Time [s]')
plt.show()

volt_hist = np.array(volt_hist)
fig, ax = plt.subplots(2,1,sharex=True)
for name, coil in mygs.coil_sets.items():
    if name.startswith('E'):
        ax[0].plot(timesteps[1:],volt_hist[:,coil['id']],label=name)
    if name.startswith('F'):
        ax[1].plot(timesteps[1:],volt_hist[:,coil['id']],label=name)
for ax_tmp in ax:
    ax_tmp.grid(True)
    ax_tmp.set_xlim(left=0.0)
    ax_tmp.set_ylim(-1.E3,1.E3)
    ax_tmp.set_ylabel('Coil Voltage [V]')
    ax_tmp.legend(ncols=2,loc='center left',bbox_to_anchor=(1.05, 0.5))
_ = ax[-1].set_xlabel('Time [ms]')
plt.show()

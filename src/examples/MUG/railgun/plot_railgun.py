#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Post-process the MUG plasma railgun example (MUG-railgun) and compare the current sheet speed with two
# estimates for a voltage-driven accelerator, where E is the applied bias field and E = v*B at the sheet:
#  - snowplow (all swept-up mass collected in the sheet): B^2/(2 mu0) = rho0 v^2, so v_sp = (E^2/(2 mu0 rho0))^(1/4)
#  - magnetic piston driving a gas-dynamic shock: B^2/(2 mu0) = p0 + rho0 u_s v, where the shock speed u_s
#    follows from the Rankine-Hugoniot relations for a piston moving at v into gas at rest
#
# Usage: run the parallel-plate and coaxial cases, then a plotting run (plot_run=T) in each directory, then
#   python plot_railgun.py --slab SLAB_DIR --coax COAX_DIR [--time T] [--outdir DIR]
#
# Each directory must contain the input file (oft_slab.in or oft_coax.in) and the plot files
# (oft_xdmf.XXXX.h5) written by the plotting run. Either case may be omitted. The script prints the
# measured and predicted sheet speeds and writes the figures MUG_railgun_ex-J.png, MUG_railgun_ex-V.png
# and MUG_railgun_ex-Sheet.png.
import argparse
import glob
import os
import re
import numpy as np
from scipy.optimize import brentq
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable

MU0 = 4.E-7*np.pi
PROTON_MASS = 1.67262192E-27
ELEM_CHARGE = 1.602176634E-19
GAMMA = 5./3.
SERIES_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
TEXT_COLOR = '#3d3d3a'
X_TRACK = (-0.2, 0.35)  # Range of sheet positions used to fit the sheet speed [m]
X_SEARCH = 0.05  # Exclude this distance from the breech (surface current) when locating the sheet [m]


def read_inputs(case_dir, case):
    '''Read the parameters needed for the snowplow estimate from the case input file'''
    with open(os.path.join(case_dir, 'oft_{0}.in'.format(case)), 'r') as fid:
        contents = fid.read()
    params = {'te_factor': 1.}
    for group, keys in (('railgun_options', ('n0', 't0', 'voltage')), ('xmhd_options', ('mu_ion', 'te_factor'))):
        group_text = re.search(r'&{0}(.*?)^\s*/'.format(group), contents, re.MULTILINE | re.DOTALL).group(1)
        for key in keys:
            match = re.search(r'^\s*{0}\s*=\s*([-+0-9.EeDd]+)'.format(key), group_text, re.MULTILINE)
            if match is not None:
                params[key] = float(match.group(1).replace('D', 'E').replace('d', 'E'))
    params['rho0'] = params['n0']*params['mu_ion']*PROTON_MASS
    params['p0'] = params['n0']*params['t0']*(1.+params['te_factor'])*ELEM_CHARGE
    return params


def read_plot_files(case_dir):
    '''Read vertex locations and fields from all per-rank plot files'''
    files = sorted(glob.glob(os.path.join(case_dir, 'oft_xdmf.*.h5')))
    if len(files) == 0:
        raise FileNotFoundError('No plot files found in "{0}" (run with plot_run=T first)'.format(case_dir))
    r = []
    fields = {}
    times = None
    for filename in files:
        with h5py.File(filename, 'r') as h5_file:
            mesh = h5_file['mug/vmesh']
            r.append(np.asarray(mesh['R']))
            steps = sorted(key for key in mesh.keys() if key.isdigit() and int(key) > 0)
            file_times = np.array([mesh[step]['TIME'][0] for step in steps])
            if times is None:
                times = file_times
            for name in ('B', 'J', 'V', 'N', 'T'):
                fields.setdefault(name, []).append(np.array([np.asarray(mesh[step][name]) for step in steps]))
    r = np.vstack(r)
    for name in fields:
        fields[name] = np.concatenate(fields[name], axis=1)
    return times, r, fields


class RailgunCase:
    '''Fields from one case in local coordinates (x, transverse coordinate s)'''
    def __init__(self, case_dir, case):
        self.case = case
        self.params = read_inputs(case_dir, case)
        self.times, r, fields = read_plot_files(case_dir)
        self.x = r[:, 0]
        self.coax = (case == 'coax')
        if self.coax:
            # Cylindrical coordinates about the x-axis: s = radius
            self.s = np.sqrt(r[:, 1]**2 + r[:, 2]**2)
            rhat = np.stack([np.zeros_like(self.s), r[:, 1]/self.s, r[:, 2]/self.s], axis=1)
            that = np.stack([np.zeros_like(self.s), -r[:, 2]/self.s, r[:, 1]/self.s], axis=1)
            self.j_s = np.einsum('tnk,nk->tn', fields['J'], rhat)/MU0  # Radial current density [A/m^2]
            self.b_perp = np.einsum('tnk,nk->tn', fields['B'], that)  # Azimuthal field
            self.in_plane = np.abs(r[:, 2]) < 1.E-6  # Points in the theta = 0 plane
            self.gap = (self.s.min(), self.s.max())
        else:
            # Parallel plates: s = y
            self.s = r[:, 1]
            self.j_s = fields['J'][:, :, 1]/MU0  # Plot files store curl(B) = mu0*J
            self.b_perp = fields['B'][:, :, 2]
            self.in_plane = np.abs(r[:, 2]) < 1.E-6  # Points in the z = 0 plane
            self.gap = (self.s.min(), self.s.max())
        self.vx = fields['V'][:, :, 0]
        self.n = fields['N']
        self.temp = fields['T']
        # Orient so that the driven current is positive
        sign = np.sign(np.sum(self.j_s[-1, self.x < 0.]))
        self.j_s *= sign
        self.b_perp *= sign

    def efield(self, s):
        '''Applied bias field at transverse location s'''
        if self.coax:
            return self.params['voltage']/(s*np.log(self.gap[1]/self.gap[0]))
        return self.params['voltage']/(self.gap[1]-self.gap[0])

    def snowplow_speed(self, s):
        '''Voltage-driven snowplow speed at transverse location s'''
        return (self.efield(s)**2/(2.*MU0*self.params['rho0']))**0.25

    def piston_speed(self, s):
        '''Speed of a magnetic piston driving a shock at transverse location s'''
        rho0 = self.params['rho0']
        p0 = self.params['p0']
        c0 = np.sqrt(GAMMA*p0/rho0)
        efield = self.efield(s)

        def shock_speed(v):
            return (GAMMA+1.)/4.*v + np.sqrt(((GAMMA+1.)/4.*v)**2 + c0**2)

        def balance(v):
            return (efield/v)**2/(2.*MU0) - p0 - rho0*shock_speed(v)*v
        return brentq(balance, 1.E-3*c0, 1.E3*c0)

    def flat_sheet_speed(self):
        '''Snowplow speed for a flat sheet spanning the whole gap'''
        if self.coax:
            ri, ro = self.gap
            return (self.params['voltage']**2/(MU0*self.params['rho0']*(ro**2-ri**2)*np.log(ro/ri)))**0.25
        return self.snowplow_speed(0.)

    def profile(self, field, itime, s0, width=8.E-3, nbins=200):
        '''Average a field over points near transverse location s0, binned in x'''
        mask = np.abs(self.s - s0) < width
        edges = np.linspace(self.x.min(), self.x.max(), nbins+1)
        index = np.clip(np.digitize(self.x[mask], edges)-1, 0, nbins-1)
        counts = np.bincount(index, minlength=nbins)
        values = np.bincount(index, weights=field[itime, mask], minlength=nbins)
        centers = (edges[1:]+edges[:-1])/2.
        valid = counts > 0
        return centers[valid], values[valid]/counts[valid]

    def sheet_position(self, s0):
        '''Location of the current density peak at transverse location s0 for each time'''
        positions = np.zeros(self.times.shape)
        for itime in range(self.times.shape[0]):
            xbins, jbins = self.profile(self.j_s, itime, s0)
            search = xbins > self.x.min() + X_SEARCH
            xbins, jbins = xbins[search], jbins[search]
            imax = np.argmax(jbins)
            # Refine with a parabola through the peak and its neighbors
            if 0 < imax < jbins.shape[0]-1:
                y0, y1, y2 = jbins[imax-1:imax+2]
                denom = y0 - 2.*y1 + y2
                shift = 0.5*(y0-y2)/denom if denom != 0. else 0.
                positions[itime] = xbins[imax] + shift*(xbins[1]-xbins[0])
            else:
                positions[itime] = xbins[imax]
        return positions


def fit_speed(times, positions):
    '''Fit a constant speed to the sheet positions inside the tracking window'''
    mask = (positions > X_TRACK[0]) & (positions < X_TRACK[1])
    if np.count_nonzero(mask) < 3:
        return None, None, mask
    speed, offset = np.polyfit(times[mask], positions[mask], 1)
    return speed, offset, mask


def plot_maps(cases, field_name, itimes, filename, label):
    '''Plot a field in the (x, s) plane for each case'''
    fig, axes = plt.subplots(len(cases), 1, figsize=(8., 1.25*len(cases)+0.4), squeeze=False)
    for ax, case, itime in zip(axes[:, 0], cases, itimes):
        field = case.j_s if field_name == 'J' else case.vx
        scale = 1.E-6 if field_name == 'J' else 1.E-3
        mask = case.in_plane
        tri = mtri.Triangulation(case.x[mask], case.s[mask])
        values = field[itime, mask]*scale
        vmax = np.percentile(np.abs(values), 99.5)
        cont = ax.tricontourf(tri, values, levels=np.linspace(0., vmax, 41), cmap='Blues', extend='both')
        ax.set_aspect('equal')
        ax.set_xlim(case.x.min(), case.x.max())
        ax.set_ylabel('r [m]' if case.coax else 'y [m]')
        title = '{0}, t = {1:.1f} $\\mu$s'.format('Coaxial' if case.coax else 'Parallel plate', case.times[itime]*1.E6)
        ax.set_title(title, fontsize=10, loc='left', color=TEXT_COLOR)
        cax = make_axes_locatable(ax).append_axes('right', size='1.5%', pad=0.08)
        ticks = [tick for tick in matplotlib.ticker.MaxNLocator(3).tick_values(0., vmax) if tick <= vmax]
        cbar = fig.colorbar(cont, cax=cax, ticks=ticks)
        cbar.set_label(label, fontsize=9)
    axes[-1, 0].set_xlabel('x [m]')
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Post-process the MUG plasma railgun example')
    parser.add_argument('--slab', type=str, default=None, help='Directory with the parallel-plate run')
    parser.add_argument('--coax', type=str, default=None, help='Directory with the coaxial run')
    parser.add_argument('--time', type=float, default=None, help='Time for field maps [s] (default: last plot time)')
    parser.add_argument('--outdir', type=str, default='.', help='Directory for figures')
    options = parser.parse_args()

    cases = []
    if options.coax is not None:
        cases.append(RailgunCase(options.coax, 'coax'))
    if options.slab is not None:
        cases.append(RailgunCase(options.slab, 'slab'))
    if len(cases) == 0:
        parser.error('At least one of --slab or --coax is required')

    # Sheet trajectories and speeds
    fig, axes = plt.subplots(1, 2, figsize=(10., 3.8))
    series = []
    for case in cases:
        if case.coax:
            ri, ro = case.gap
            for frac in (0.1, 0.5, 0.8):
                s0 = ri + frac*(ro-ri)
                series.append((case, s0, 'Coax, r = {0:.2f} m'.format(s0)))
        else:
            series.append((case, 0., 'Parallel plate'))
    print('Sheet speed [km/s]')
    print('{0:<24} {1:>10} {2:>10} {3:>10}'.format('Location', 'Measured', 'Piston', 'Snowplow'))
    for (case, s0, name), color in zip(series, SERIES_COLORS):
        positions = case.sheet_position(s0)
        speed, offset, mask = fit_speed(case.times, positions)
        v_sp = case.snowplow_speed(s0)
        v_piston = case.piston_speed(s0)
        if speed is None:
            print('{0:<24} {1:>10} {2:>10.1f} {3:>10.1f}'.format(name, 'n/a', v_piston/1.E3, v_sp/1.E3))
            continue
        print('{0:<24} {1:>10.1f} {2:>10.1f} {3:>10.1f}'.format(name, speed/1.E3, v_piston/1.E3, v_sp/1.E3))
        times_us = case.times*1.E6
        axes[0].plot(times_us, positions, 'o', color=color, markersize=4,
                     label='{0}: {1:.0f} km/s ({2:.0f})'.format(name, speed/1.E3, v_piston/1.E3))
        t_fit = case.times[mask]
        axes[0].plot(t_fit*1.E6, offset + v_piston*(t_fit - t_fit[0]) + speed*t_fit[0], '--', color=color, linewidth=1.5)
    for case in cases:
        if case.coax:
            print('Coax, flat sheet: snowplow {0:.1f} km/s'.format(case.flat_sheet_speed()/1.E3))
    axes[0].axhspan(X_TRACK[0], X_TRACK[1], color='#e8e8e3', zorder=0)
    axes[0].set_xlabel(r't [$\mu$s]')
    axes[0].set_ylabel('Current sheet position [m]')
    axes[0].set_title('Sheet position and speed (dashed/parentheses: piston estimate)', fontsize=10, loc='left', color=TEXT_COLOR)
    axes[0].legend(fontsize=8, frameon=False, loc='upper left')

    # Centerline profiles for the first parallel-plate (or coaxial) case
    case = cases[-1]
    s0 = 0. if not case.coax else 0.5*(case.gap[0]+case.gap[1])
    itime = case.times.shape[0]//2 if options.time is None else np.argmin(np.abs(case.times-options.time))
    xbins, bbins = case.profile(case.b_perp, itime, s0)
    _, nbins = case.profile(case.n, itime, s0)
    _, vbins = case.profile(case.vx, itime, s0)
    b_sp = case.efield(s0)/case.snowplow_speed(s0)
    axes[1].plot(xbins, bbins/b_sp, color=SERIES_COLORS[0], linewidth=2., label=r'$B / B_{sp}$')
    axes[1].plot(xbins, nbins/case.params['n0'], color=SERIES_COLORS[1], linewidth=2., label=r'$n / n_0$')
    axes[1].plot(xbins, vbins/case.snowplow_speed(s0), color=SERIES_COLORS[2], linewidth=2., label=r'$V_x / v_{sp}$')
    axes[1].set_xlabel('x [m]')
    axes[1].set_title('{0} centerline, t = {1:.1f} $\\mu$s'.format('Coaxial' if case.coax else 'Parallel plate',
                      case.times[itime]*1.E6), fontsize=10, loc='left', color=TEXT_COLOR)
    axes[1].legend(fontsize=8, frameon=False)
    for ax in axes:
        ax.grid(color='#e8e8e3', linewidth=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(options.outdir, 'MUG_railgun_ex-Sheet.png'), dpi=150)
    plt.close(fig)

    # Field maps
    if options.time is None:
        itimes = [c.times.shape[0]-1 for c in cases]
    else:
        itimes = [np.argmin(np.abs(c.times-options.time)) for c in cases]
    plot_maps(cases, 'J', itimes, os.path.join(options.outdir, 'MUG_railgun_ex-J.png'),
              r'$J_{\perp}$ [MA/m$^2$]')
    plot_maps(cases, 'V', itimes, os.path.join(options.outdir, 'MUG_railgun_ex-V.png'),
              r'$V_x$ [km/s]')


if __name__ == '__main__':
    main()

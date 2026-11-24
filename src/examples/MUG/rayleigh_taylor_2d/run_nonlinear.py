#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Run the 2D Rayleigh-Taylor example (MUG_rt_2d) with the tanh transition layer, linearized and nonlinear, for
# one or more in-plane field strengths, and compare the growth of the seeded mode.
#
# The nonlinear runs start from the projected equilibrium, which is not in exact discrete balance, so the total
# kinetic energy also contains a small, x-independent adjustment of the layer. The comparison therefore uses the
# amplitude of the cos(kx) component of v_z (k = 2 pi/L_x), computed from the plot files of each run. Particle
# diffusion and thermal conduction are turned off (D_diff = chi = 0), since they also spread the layer in the
# nonlinear runs, which the linearized model does not include; the other options come from oft.in.
#
# Usage (from the directory containing MUG_rt_2d, oft.in, oft_in.xml, and oft_in_mpi.xml for --np > 1):
#   python run_nonlinear.py [--B_0 0. 2.E-4] [--nsteps_lin 400] [--np 1] [--plot] [--skip_run]
import argparse
import os
import shutil
import subprocess
import numpy as np
from OpenFUSIONToolkit.io import XDMF_plot_file
from run_scan import read_option, set_option, parallel_settings, run_command
from plot_rt_2d import plot_files


def mode_amplitude(run_dir, k):
    '''Return (t, a): the amplitude of the cos(kx) component of v_z over the whole domain at each plot time'''
    times, amps = None, []
    filenames = plot_files(run_dir)
    if len(filenames) == 0:
        raise FileNotFoundError('No plot files in {0} (see {0}/run.log and plot.log)'.format(run_dir))
    meshes = [XDMF_plot_file(f)['xmhd_2d']['smesh'] for f in filenames]
    times = meshes[0].times
    for istep in range(times.shape[0]):
        num, den = 0.0, 0.0
        for mesh in meshes:
            c = np.cos(k*mesh.r[:, 0])
            # The plot vector is (vx, vz, vy): the second component is the vertical velocity
            num += np.sum(mesh.get_field('V', timestep=istep)[:, 1]*c)
            den += np.sum(c**2)
        amps.append(num/den)
    return np.asarray(times), np.abs(np.asarray(amps))


def run(exe, text, run_dir, nprocs):
    '''Run a case and convert its restart files to plot files'''
    os.makedirs(run_dir, exist_ok=True)
    text, xml_file = parallel_settings(text, nprocs)
    shutil.copy(xml_file, run_dir)
    cmd = run_command(exe, xml_file, nprocs)
    for plot in ('F', 'T'):
        with open(os.path.join(run_dir, 'oft.in'), 'w') as fid:
            fid.write(set_option(text, 'plot_run', plot))
        with open(os.path.join(run_dir, 'run.log' if plot == 'F' else 'plot.log'), 'w') as log:
            result = subprocess.run(cmd, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print('  WARNING: {0} stopped early (see {1})'.format('plotting' if plot == 'T' else 'run', run_dir))


def fit_rate(t, a, a_max):
    '''Growth rate from the slope of ln(a) over the part of the run where 0.03 a_max < a < 0.3 a_max'''
    mask = (a > 0.03*a_max) & (a < 0.3*a_max)
    if np.count_nonzero(mask) < 3:
        print('  WARNING: fewer than 3 plot times in the fit window (run longer, or change --a_max)')
        return np.nan, np.nan, np.nan
    return np.polyfit(t[mask], np.log(a[mask]), 1)[0], t[mask].min(), t[mask].max()


def main():
    parser = argparse.ArgumentParser(description='Linear vs nonlinear Rayleigh-Taylor growth (tanh layer)')
    parser.add_argument('--B_0', type=float, nargs='+', default=[0.0, 2.E-4], help='In-plane fields [T]')
    parser.add_argument('--nsteps_lin', type=int, default=400, help='Time steps for the linear runs')
    parser.add_argument('--a_max', type=float, default=100.0,
                        help='Mode amplitude [m/s] that sets the fit window (0.03-0.3 a_max) in both runs')
    parser.add_argument('--np', type=int, default=1, help='Number of MPI processes per run')
    parser.add_argument('--exe', default='./MUG_rt_2d', help='Path to the MUG_rt_2d executable')
    parser.add_argument('--plot', action='store_true', help='Plot the mode amplitudes (rt_2d_nonlinear.png)')
    parser.add_argument('--skip_run', action='store_true',
                        help='Analyze existing runs (and plot files) instead of running the cases')
    options = parser.parse_args()
    with open('oft.in', 'r') as fid:
        template = fid.read()
    if read_option(template, 'V1')[0] > 2.0:
        print('Note: a seed V1 above ~1 m/s leaves little time before the mode becomes nonlinear')
    # The linearized model holds the equilibrium fixed, so compare without particle diffusion and thermal conduction,
    # which spread the layer in the nonlinear run only
    template = set_option(template, 'D_diff', '0.')
    template = set_option(template, 'chi', '0.')
    k = 2.0*np.pi/read_option(template, 'rscale')[0]
    exe = os.path.abspath(options.exe)
    results = []
    for B_0 in options.B_0:
        text = set_option(template, 'B_0', '{0:.6E}'.format(B_0))
        curves = {}
        for linear in (True, False):
            run_dir = 'tanh_B_{0:.2E}_{1}'.format(B_0, 'lin' if linear else 'nl')
            case = set_option(text, 'linear', 'T' if linear else 'F')
            if linear:
                case = set_option(case, 'nsteps', str(options.nsteps_lin))
            if not options.skip_run:
                run(exe, case, run_dir, options.np)
            curves[linear] = mode_amplitude(run_dir, k)
        # Fit both runs over the same time window, set by the nonlinear run's amplitude
        gamma_nl, t0, t1 = fit_rate(*curves[False], options.a_max)
        t_lin, a_lin = curves[True]
        mask = (t_lin >= t0) & (t_lin <= t1)
        gamma_lin = np.polyfit(t_lin[mask], np.log(a_lin[mask]), 1)[0] if np.count_nonzero(mask) >= 2 else np.nan
        results.append((B_0, gamma_lin, gamma_nl, t0, t1, curves))
    print('  B_0 [T]    window [us]      gamma linear   gamma nonlinear   ratio')
    for B_0, gamma_lin, gamma_nl, t0, t1, _ in results:
        print('  {0:8.2E}  {1:6.0f}-{2:<6.0f}  {3:13.4E}  {4:16.4E}  {5:7.4f}'.format(
            B_0, t0*1.E6, t1*1.E6, gamma_lin, gamma_nl, gamma_nl/gamma_lin))
    if options.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        for i, (B_0, _, _, t0, t1, curves) in enumerate(results):
            color = 'C{0}'.format(i)
            ax.semilogy(curves[True][0]*1.E6, curves[True][1], '--', color=color,
                        label=r'Linear, $B_0$ = {0:.1E} T'.format(B_0))
            ax.semilogy(curves[False][0]*1.E6, curves[False][1], '-', color=color,
                        label=r'Nonlinear, $B_0$ = {0:.1E} T'.format(B_0))
        ax.axhspan(0.03*options.a_max, 0.3*options.a_max, color='0.9')
        ax.set_xlabel(r'$t$ [$\mu$s]')
        ax.set_ylabel(r'$|\hat{v}_z|$, $\cos(kx)$ component [m/s]')
        ax.grid(True)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig('rt_2d_nonlinear.png', dpi=100)


if __name__ == '__main__':
    main()

#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Plot the density from one or more runs of the 2D Rayleigh-Taylor example (MUG_rt_2d) at a set of times, and
# report the minimum density and temperature over each run.
#
# Each run directory must contain the XDMF/HDF5 plot files written with `plot_run=T` (oft_xdmf.*.h5). One row is
# drawn per run and one column per time; the nearest plotted time is used. Times are given either as one list for
# all runs, or as one comma-separated group per run (each group with the same number of times).
#
# Usage:
#   python plot_rt_2d.py --runs b0 b05 --labels "B_0 = 0" "B_0 = 0.5 B_c" --times 2.E-4 3.E-4 4.E-4 --out rt_2d_nl.png
#   python plot_rt_2d.py --runs b0 b05 --times 3.E-4,4.E-4 4.E-4,5.E-4
import argparse
import glob
import os
import numpy as np
from OpenFUSIONToolkit.io import XDMF_plot_file


def plot_files(run_dir):
    '''Plot files of a run (one per MPI process), in `run_dir` or, from older versions, in `run_dir/all`'''
    for sub_dir in ('', 'all'):
        files = sorted(glob.glob(os.path.join(run_dir, sub_dir, 'oft_xdmf.*.h5')))
        if len(files) > 0:
            return files
    return []


def load_run(run_dir):
    '''Load the plot meshes of a run (one per MPI process)'''
    files = plot_files(run_dir)
    if len(files) == 0:
        raise FileNotFoundError('No plot files in {0} (run with plot_run=T first)'.format(run_dir))
    return [XDMF_plot_file(filename)['xmhd_2d']['smesh'] for filename in files]


def field_at(meshes, name, istep):
    '''Return (r, lc, values) for a field at plot step `istep`, merged over all processes'''
    r_list, lc_list, val_list = [], [], []
    offset = 0
    for mesh in meshes:
        r_list.append(mesh.r)
        lc_list.append(mesh.lc + offset)
        val_list.append(mesh.get_field(name, timestep=istep))
        offset += mesh.np
    return np.vstack(r_list), np.vstack(lc_list), np.concatenate(val_list)


def main():
    parser = argparse.ArgumentParser(description='Density snapshots from 2D Rayleigh-Taylor runs')
    parser.add_argument('--runs', nargs='+', default=['.'], help='Run directories (one row each)')
    parser.add_argument('--labels', nargs='+', default=None, help='Row labels')
    parser.add_argument('--times', nargs='+', default=None,
                        help='Times to plot [s]: one list for all runs, or one comma-separated group per run '
                             '(default: last time)')
    parser.add_argument('--out', default='rt_2d_density.png', help='Output image')
    parser.add_argument('--dpi', type=int, default=100, help='Image resolution')
    options = parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    labels = options.labels if options.labels is not None else options.runs
    runs = [load_run(run_dir) for run_dir in options.runs]
    if options.times is None:
        row_times = [[meshes[0].times[-1]] for meshes in runs]
    elif all(',' in group for group in options.times):
        row_times = [[float(val) for val in group.split(',')] for group in options.times]
        if len(row_times) != len(runs) or len(set(len(group) for group in row_times)) != 1:
            raise ValueError('Give one group of times per run, all with the same number of times')
    else:
        row_times = [[float(val) for val in options.times] for _ in runs]
    ncols = len(row_times[0])
    # Minimum density and temperature over each run
    for label, meshes in zip(labels, runs):
        n_min, T_min = np.inf, np.inf
        for istep in range(meshes[0].times.shape[0]):
            n_min = min(n_min, field_at(meshes, 'n', istep)[2].min())
            T_min = min(T_min, field_at(meshes, 'T', istep)[2].min())
        print('{0}: t_max = {1:.3E} s, min(n) = {2:.4E} m^-3, min(T) = {3:.4E} eV'.format(
            label, meshes[0].times[-1], n_min, T_min))
    fig, axes = plt.subplots(len(runs), ncols, figsize=(2.2*ncols+0.8, 2.3*len(runs)), squeeze=False,
                             sharex=True, sharey=True)
    n_range = None
    for i, (label, meshes) in enumerate(zip(labels, runs)):
        for j, time in enumerate(row_times[i]):
            istep = int(np.argmin(np.abs(meshes[0].times - time)))
            r, lc, n = field_at(meshes, 'n', istep)
            if n_range is None:
                n_range = (n.min(), n.max())
            ax = axes[i, j]
            tri = mtri.Triangulation(r[:, 0], r[:, 1], lc)
            image = ax.tripcolor(tri, n, shading='gouraud', cmap='viridis', vmin=n_range[0], vmax=n_range[1],
                                 rasterized=True)
            ax.set_xlim(r[:, 0].min(), r[:, 0].max())
            ax.set_ylim(r[:, 1].min(), r[:, 1].max())
            ax.set_aspect('equal')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(r'$t$ = {0:.0f} $\mu$s'.format(meshes[0].times[istep]*1.E6), fontsize=9)
            if j == 0:
                ax.set_ylabel(label, fontsize=10)
    cbar = fig.colorbar(image, ax=axes, shrink=0.9, pad=0.02)
    cbar.set_label(r'$n$ [m$^{-3}$]')
    fig.savefig(options.out, dpi=options.dpi, bbox_inches='tight')


if __name__ == '__main__':
    main()

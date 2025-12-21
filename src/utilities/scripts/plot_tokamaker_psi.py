import sys
import numpy as np
from scipy.interpolate import griddata
import h5py
import argparse
#
parser = argparse.ArgumentParser()
parser.description = "Flux contour plotting script for TokaMaker"
parser.add_argument("--rstfile", default='tokamaker_gs.rst', type=str, help='Name of TokaMaker data file (default: "psi_gs.rst")')
parser.add_argument("--nplevels", default=10, type=int, help="Number of levels for plasma region contour plot")
parser.add_argument("--nvlevels", default=20, type=int, help="Number of levels for vacuum region contour plot")
parser.add_argument("--title", default=None, type=str, help="Title for contour plot")
parser.add_argument("--filename", default=None, type=str, help="Filename to save plot (without extension, saved as PNG)")
parser.add_argument("--press_color", action="store_true", default=False, help="Use pressure to shade interior")
parser.add_argument("--plot_b_perp", action="store_true", default=False, help="Plot B_perp on mid-plane")
parser.add_argument("--no_color", action="store_true", default=False, help="Do not shade interior")
parser.add_argument("--no_saddle", action="store_true", default=False, help="Do not mark saddle points (O-, X-points)")
parser.add_argument("--no_isoflux", action="store_true", default=False, help="Do not mark isoflux points (if present)")
parser.add_argument("--pdf_fig", action="store_true", default=False, help="Create PDF figure instead of PNG")
options = parser.parse_args()
if options.filename is not None:
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Load data from TokaMaker equilibrium restart file
with h5py.File(options.rstfile, 'r') as fid:
    lc = np.asarray(fid['mesh/lc'])-1
    r = np.asarray(fid['mesh/r'])
    regions = np.asarray(fid['mesh/regions'])
    lc_plot = np.asarray(fid['mesh/lc_plot'])
    r_plot = np.asarray(fid['mesh/r_plot'])
    psi = np.asarray(fid['gs/psi'])
    psi_bounds = np.asarray(fid['gs/bounds'])
    o_point = np.asarray(fid['gs/o_point'])
    p_profile = np.asarray(fid['/gs/p/sample'])
    if 'gs/nx_points' in fid:
        nx_points = int(fid['gs/nx_points'][0])
        x_points = np.asarray(fid['gs/x_points'])
    else:
        nx_points = 0
    if ('gs/isoflux_ntargets' in fid) and (not options.no_isoflux):
        isoflux_ntargets = int(fid['gs/isoflux_ntargets'][0])
        isoflux_targets = np.asarray(fid['gs/isoflux_targets'])
    else:
        isoflux_ntargets = 0
    if 'gs/vac_flag' in fid:
        vac_flag = bool(fid['gs/vac_flag'][0])
    else:
        vac_flag = False
    if 'gs/B' in fid:
        B = np.asarray(fid['gs/B'])
    else:
        B = None
if options.press_color:
    psi_sample = np.linspace(psi_bounds[0], psi_bounds[1], p_profile.shape[0])
    shading = np.interp(psi, psi_sample, p_profile[:,1], left=0.0, right=p_profile[-1,1])
else:
    shading = psi
# Plot shading (flux or pressure)
fig, ax = plt.subplots(tight_layout=True)
if not options.no_color:
    if options.press_color:
        psi_sample = np.linspace(psi_bounds[0], psi_bounds[1], p_profile.shape[0])
        shading = np.interp(psi, psi_sample, p_profile[:,1], left=0.0, right=p_profile[-1,1])
        label = 'Pressure [Pa]'
    else:
        shading = psi
        label = 'Poloidal flux [Wb]'
    pcm = ax.tricontourf(r_plot[:,0], r_plot[:,1], lc_plot, shading, min(100,(options.nplevels+options.nvlevels)*2))
    colorbar = fig.colorbar(pcm, ax=ax)
    colorbar.set_label(label)
# Plot regions
if regions.max() > 1.:
    mask = np.where(regions > 1.5)[0]
    mask_vals = np.ones(r.shape[0])
    ax.tricontourf(r[:,0], r[:,1], lc[mask,:], mask_vals, colors='k', alpha=0.3)
# Plot contours
if vac_flag:
    ax.tricontour(r_plot[:,0], r_plot[:,1], lc_plot, psi,
        np.linspace(psi.min(),psi.max(),options.nvlevels), colors='k')
else:
    ax.tricontour(r_plot[:,0], r_plot[:,1], lc_plot, psi-psi_bounds[0],
        np.linspace(0.0,psi_bounds[1]-psi_bounds[0],options.nplevels), colors='k')
    ax.tricontour(r_plot[:,0], r_plot[:,1], lc_plot, psi-psi_bounds[0],
        np.linspace(psi.min()-psi_bounds[0],0.0,options.nvlevels), colors='k')
    # Plot O-point and X-points
    if not options.no_saddle:
        ax.plot(o_point[0], o_point[1], 'r+')
        for i in range(nx_points):
            if i == nx_points-1:
                ax.plot(x_points[i,0], x_points[i,1], 'rx')
            else:
                ax.plot(x_points[i,0], x_points[i,1], 'bx')
    # Plot isoflux constraints
    if isoflux_ntargets > 0:
        ax.plot(isoflux_targets[:,0], isoflux_targets[:,1], 'yo', fillstyle='none')
# Style plot
x1 = r_plot[:,0].min(); x2 = r_plot[:,0].max()
x1 = x1 - (x2-x1)*0.1; x2 = x2 + (x2-x1)*0.1
ax.set_xlim(max(0.,x1), x2)
y1 = r_plot[:,1].min(); y2 = r_plot[:,1].max()
y1 = y1 - (y2-y1)*0.1; y2 = y2 + (y2-y1)*0.1
ax.set_ylim(y1,y2)
ax.set_ylabel('Cylindrical Height [m]')
ax.set_aspect('equal')
#
if options.plot_b_perp:
    ax.xaxis.set_tick_params(labelbottom=False)
    divider = make_axes_locatable(ax)
    ax_perp = divider.append_axes("bottom", 1.2, pad=0.1, sharex=ax)
    r_sample = np.linspace(r_plot[:,0].min(), r_plot[:,0].max(), 1000)
    if B is None:
        dr = r_sample[1] - r_sample[0]
        psi_sample = griddata(r_plot, psi, (r_sample, np.zeros(r_sample.shape)))
        r_mid = (r_sample[1:]+r_sample[:-1])/2.0
        ax_perp.semilogy(r_mid, np.abs(np.diff(psi_sample)/dr/r_mid))
    else:
        ax_perp.semilogy(r_sample, griddata(r_plot, abs(B[:,2]), (r_sample, np.zeros(r_sample.shape))))
    ax_perp.grid(True)
    ax_perp.set_ylabel(r'$|B_{\perp}| [T]$')
    ax_perp.set_xlabel('Major Radius [m]')
else:
    ax.set_xlabel('Major Radius [m]')
#
if options.title is not None:
    fig.set_title(options.title)
# Handle saving if necessary
if options.filename is None:
    plt.show()
else:
    if options.pdf_fig:
        fig.savefig('{0}.pdf'.format(options.filename.strip()))
    else:
        fig.savefig('{0}.png'.format(options.filename.strip()), dpi=400)

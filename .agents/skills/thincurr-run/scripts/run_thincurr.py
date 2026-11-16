#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-only
"""Run a ThinCurr calculation on an existing model from a JSON case definition.

Usage:
  python run_thincurr.py case.json [--output-dir DIR] [--nthreads N] [--no-plot]

Supported analyses ("analysis" in the case file):
  eig   L/R eigenmodes (current decay times and mode structures)

See `../reference.md` for the case-file schema and `../examples/` for complete
examples. Build and verify the model first (thincurr-model skill).
"""
import argparse
import json
import os
import sys
import time

import numpy

# Above this many unknowns HODLR is used when "hodlr" is "auto"
HODLR_AUTO_NELEMS = 10000
# At or below this many unknowns a dense model uses the direct (LAPACK) eigensolver when "direct" is "auto"
DIRECT_AUTO_NELEMS = 4000


def _load_thincurr():
    """Import the OpenFUSIONToolkit python package, honoring OFT_ROOTPATH.

    Returns:
        Tuple of (OFT_env, ThinCurr) classes.
    """
    root_path = os.getenv('OFT_ROOTPATH')
    if root_path is not None:
        sys.path.append(os.path.join(root_path, 'python'))
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.ThinCurr import ThinCurr
    return OFT_env, ThinCurr


def resolve_path(path, base_dir):
    """Resolve a path from the case file relative to the case file's directory.

    Args:
        path: Path as given in the case file (or None).
        base_dir: Directory containing the case file.

    Returns:
        Absolute path, or None.
    """
    if path is None:
        return None
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base_dir, path))


def _insulated_junctions(spec):
    """Convert JSON insulated-junction entries (ints or 2-lists) to setup_model() form.

    Args:
        spec: List from the case file, or None.

    Returns:
        List of region ints and region-pair tuples, or None.
    """
    if not spec:
        return None
    return [tuple(item) if isinstance(item, (list, tuple)) else int(item) for item in spec]


def setup_model(case, base_dir, out_dir, prefix, nthreads):
    """Create the ThinCurr object, load the model, and apply resistivity overrides.

    Args:
        case: Full case dictionary.
        base_dir: Directory containing the case file.
        out_dir: Output directory.
        prefix: Output file prefix.
        nthreads: OpenMP thread count.

    Returns:
        Tuple of (ThinCurr instance, geometry dictionary with "r", "lc", "reg" or None).
    """
    OFT_env, ThinCurr = _load_thincurr()
    model = case['model']
    myOFT = OFT_env(nthreads=nthreads)
    tw_model = ThinCurr(myOFT)
    mesh_file = resolve_path(model['mesh_file'], base_dir)
    tw_model.setup_model(mesh_file=mesh_file,
                         xml_filename=resolve_path(model.get('xml_file'), base_dir),
                         jumper_start=int(model.get('jumper_start', 0)),
                         insulated_junctions=_insulated_junctions(model.get('insulated_junctions')))
    tw_model.setup_io(basepath=os.path.join(out_dir, prefix+'_plot'))
    # Optional resistivity override (eg. for parameter scans without editing the model)
    eta_spec = model.get('eta')
    if eta_spec:
        kwargs = {key: numpy.asarray(eta_spec[key], dtype=numpy.float64)
                  for key in ('eta_surf', 'eta_vol', 'thickness') if key in eta_spec}
        tw_model.set_eta_values(**kwargs)
    eta_surf, _ = tw_model.get_eta_values(include_eta_vol=True)
    if numpy.any(eta_surf <= 0.0):
        raise ValueError('Surface resistivity is missing or non-positive for regions {0}; fix the model '
                         '(thincurr-model skill) or set "model.eta"'.format(
                             (numpy.where(eta_surf <= 0.0)[0]+1).tolist()))
    # Mesh geometry for diagnostics and plots (the ThinCurr object does not expose it for file models)
    from OpenFUSIONToolkit.ThinCurr.meshing import read_ThinCurr_mesh
    mesh_info = read_ThinCurr_mesh(mesh_file)
    geom = None
    if (mesh_info['r'].shape[0] == tw_model.np) and (mesh_info['lc'].shape[0] == tw_model.nc):
        geom = {'r': mesh_info['r'], 'lc': mesh_info['lc'], 'reg': numpy.asarray(mesh_info['reg'], dtype=numpy.int32)}
    else:
        print('WARNING: mesh file and model sizes differ; skipping per-region diagnostics and plots')
    return tw_model, geom


def build_operators(tw_model, case, out_dir, prefix):
    """Compute the inductance (dense or HODLR) and resistance operators.

    Args:
        tw_model: ThinCurr instance after setup_model().
        case: Full case dictionary.
        out_dir: Output directory.
        prefix: Output file prefix.

    Returns:
        Dictionary describing the operators (HODLR flag, cache file, timings).
    """
    hodlr_spec = case.get('hodlr', 'auto')
    if hodlr_spec == 'auto':
        use_hodlr = bool(tw_model.nelems > HODLR_AUTO_NELEMS)
    else:
        use_hodlr = bool(hodlr_spec)
    if use_hodlr and case.get('hodlr_options'):
        tw_model.set_hodlr_options(**case['hodlr_options'])
    cache_file = None
    if case.get('cache', True):
        # Cache files are validated against mesh/coil/HODLR-setting hashes, so stale caches are recomputed
        cache_file = os.path.join(out_dir, '{0}_L{1}.save'.format(prefix, '_HODLR' if use_hodlr else ''))
    print('\nComputing L ({0}, {1} unknowns)'.format('HODLR' if use_hodlr else 'dense', tw_model.nelems), flush=True)
    t0 = time.perf_counter()
    tw_model.compute_Lmat(use_hodlr=use_hodlr, cache_file=cache_file)
    t1 = time.perf_counter()
    tw_model.compute_Rmat()
    t2 = time.perf_counter()
    return {'hodlr': use_hodlr, 'L_cache': cache_file, 'time_L': t1-t0, 'time_R': t2-t1}


def region_fractions(tw_model, geom, vec):
    """Fraction of a mode's resistive dissipation in each region.

    Args:
        tw_model: ThinCurr instance.
        geom: Geometry dictionary from setup_model().
        vec: Mode (current potential) [nelems].

    Returns:
        List of fractions, one per region.
    """
    jcell = tw_model.reconstruct_current(vec, centering='cell')
    r, lc, reg = geom['r'], geom['lc'], geom['reg']
    v1 = r[lc[:, 1], :] - r[lc[:, 0], :]
    v2 = r[lc[:, 2], :] - r[lc[:, 0], :]
    areas = 0.5*numpy.linalg.norm(numpy.cross(v1, v2), axis=1)
    eta_surf, _ = tw_model.get_eta_values(include_eta_vol=True)
    power = numpy.sum(jcell**2, axis=1)*areas*eta_surf[reg-1]
    total = power.sum()
    fracs = [float(power[reg == i+1].sum()/total) if total > 0.0 else 0.0 for i in range(tw_model.nregs)]
    return fracs


def plot_modes(tw_model, geom, vecs, taus, filename):
    """Plot the current magnitude and direction of the leading modes.

    Args:
        tw_model: ThinCurr instance.
        geom: Geometry dictionary from setup_model().
        vecs: Mode vectors [nmodes, nelems].
        taus: Decay times [nmodes].
        filename: Output image path.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    r, lc = geom['r'], geom['lc']
    centers = r[lc, :].mean(axis=1)
    nmodes = vecs.shape[0]
    ncols = min(nmodes, 2)
    nrows = (nmodes + ncols - 1)//ncols
    fig = plt.figure(figsize=(6*ncols, 5*nrows))
    extent = numpy.ptp(r, axis=0)
    box = numpy.maximum(extent, 1.E-2*extent.max())
    narrows = min(lc.shape[0], 400)
    arrow_ids = numpy.linspace(0, lc.shape[0]-1, narrows).astype(int)
    for i in range(nmodes):
        jcell = tw_model.reconstruct_current(vecs[i, :], centering='cell')
        jmag = numpy.linalg.norm(jcell, axis=1)
        jmax = jmag.max() if jmag.max() > 0.0 else 1.0
        ax = fig.add_subplot(nrows, ncols, i+1, projection='3d')
        colors = plt.get_cmap('viridis')(jmag/jmax)
        ax.add_collection3d(Poly3DCollection(r[lc, :], facecolors=colors, edgecolors='none', alpha=0.9))
        scale = 0.1*extent.max()/jmax
        sel = centers[arrow_ids, :]
        jsel = jcell[arrow_ids, :]*scale
        ax.quiver(sel[:, 0], sel[:, 1], sel[:, 2], jsel[:, 0], jsel[:, 1], jsel[:, 2], color='k', lw=0.8, arrow_length_ratio=0.4)
        mid = 0.5*(r.max(axis=0) + r.min(axis=0))
        ax.set_xlim(mid[0]-box[0]/2.0, mid[0]+box[0]/2.0)
        ax.set_ylim(mid[1]-box[1]/2.0, mid[1]+box[1]/2.0)
        ax.set_zlim(mid[2]-box[2]/2.0, mid[2]+box[2]/2.0)
        ax.set_box_aspect(box.copy())  # set_box_aspect() normalizes its argument in place
        ax.set_title('Mode {0}: tau = {1:.4G} ms'.format(i+1, 1.E3*taus[i]))
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
    fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.94, wspace=0.15, hspace=0.15)
    fig.savefig(filename, dpi=130)
    plt.close(fig)


def run_eig(ctx, make_plot):
    """L/R eigenmode analysis: L I = tau R I.

    Args:
        ctx: Run context (keys "tw", "geom", "case", "ops", "out_dir", "prefix").
        make_plot: Write a PNG of the leading modes?

    Returns:
        Tuple of (results dictionary, list of written files).
    """
    tw_model, geom, case, ops = ctx['tw'], ctx['geom'], ctx['case'], ctx['ops']
    out_dir, prefix = ctx['out_dir'], ctx['prefix']
    spec = case.get('eig', {})
    neigs = int(spec.get('neigs', 5))
    if neigs >= tw_model.nelems:
        raise ValueError('"neigs" ({0}) must be smaller than the number of unknowns ({1})'.format(neigs, tw_model.nelems))
    direct = spec.get('direct', 'auto')
    if direct == 'auto':
        direct = (not ops['hodlr']) and bool(tw_model.nelems <= DIRECT_AUTO_NELEMS)
    direct = bool(direct)
    if direct and ops['hodlr']:
        raise ValueError('"direct": true is not supported with HODLR; use the iterative (ARPACK) solver')
    print('\nComputing {0} eigenmodes ({1} solver)'.format(neigs, 'direct' if direct else 'ARPACK'), flush=True)
    t0 = time.perf_counter()
    eig_vals, eig_vecs = tw_model.get_eigs(neigs, direct)
    t_eig = time.perf_counter() - t0
    # Slowest-decaying (longest-lived) modes first
    order = numpy.argsort(eig_vals)[::-1]
    eig_vals = eig_vals[order]
    eig_vecs = eig_vecs[order, :]

    written = []
    nsave = min(int(spec.get('save_modes', 4)), neigs)
    modes = []
    for i in range(neigs):
        mode = {'index': i+1, 'tau_s': float(eig_vals[i]), 'tau_ms': float(1.E3*eig_vals[i])}
        if (tw_model.nregs > 1) and (geom is not None):
            mode['region_dissipation_fraction'] = region_fractions(tw_model, geom, eig_vecs[i, :])
        modes.append(mode)
    print('\nL/R decay times:')
    for mode in modes:
        line = '  {0:4d}  tau = {1:12.5E} s'.format(mode['index'], mode['tau_s'])
        if 'region_dissipation_fraction' in mode:
            line += '  region fractions = {0}'.format(
                ', '.join('{0:.2f}'.format(frac) for frac in mode['region_dissipation_fraction']))
        print(line)

    if case.get('output', {}).get('xdmf', True) and nsave > 0:
        for i in range(nsave):
            tw_model.save_current(eig_vecs[i, :], 'J_{0:02d}'.format(i+1))
        tw_model.build_XDMF()
        written.append(tw_model._io_basepath)
    if case.get('output', {}).get('vectors', False):
        vec_file = os.path.join(out_dir, prefix+'_eigvecs.npy')
        numpy.save(vec_file, eig_vecs)
        written.append(vec_file)
    if make_plot and (nsave > 0) and (geom is not None):
        plot_file = os.path.join(out_dir, prefix+'_modes.png')
        plot_modes(tw_model, geom, eig_vecs[:nsave, :], eig_vals[:nsave], plot_file)
        written.append(plot_file)
    results = {'neigs': neigs, 'solver': 'direct' if direct else 'arpack', 'time_eig': t_eig, 'modes': modes}
    return results, written


# Analyses available through "analysis"; time-domain ("td") and frequency-response ("fr")
# runs belong here and reuse setup_model() and build_operators()
ANALYSES = {
    'eig': run_eig,
}


def main():
    parser = argparse.ArgumentParser(description='Run a ThinCurr calculation from a JSON case file')
    parser.add_argument('case_file', help='JSON case definition')
    parser.add_argument('--output-dir', default=None, help='Directory for output files (default: case "output.dir", else the working directory)')
    parser.add_argument('--nthreads', type=int, default=None, help='OpenMP thread count (default: case file or 2)')
    parser.add_argument('--no-plot', action='store_true', help='Skip writing PNG plots')
    args = parser.parse_args()

    with open(args.case_file, 'r') as fid:
        case = json.load(fid)
    base_dir = os.path.dirname(os.path.abspath(args.case_file))
    out_spec = case.get('output', {})
    # Outputs go to --output-dir, else "output.dir" (relative to the case file), else the working directory
    if args.output_dir is not None:
        out_dir = os.path.abspath(args.output_dir)
    elif 'dir' in out_spec:
        out_dir = resolve_path(out_spec['dir'], base_dir)
    else:
        out_dir = os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_spec.get('prefix', os.path.splitext(os.path.basename(args.case_file))[0])
    analysis = case.get('analysis', 'eig')
    if analysis not in ANALYSES:
        raise ValueError('Unknown analysis "{0}" (supported: {1})'.format(analysis, ', '.join(ANALYSES)))
    if 'model' not in case or 'mesh_file' not in case['model']:
        raise ValueError('"model.mesh_file" is required')

    nthreads = args.nthreads or int(case.get('nthreads', 2))
    tw_model, geom = setup_model(case, base_dir, out_dir, prefix, nthreads)
    ops = build_operators(tw_model, case, out_dir, prefix)
    ctx = {'tw': tw_model, 'geom': geom, 'case': case, 'ops': ops, 'out_dir': out_dir, 'prefix': prefix}
    results, written = ANALYSES[analysis](ctx, out_spec.get('plot', True) and (not args.no_plot))

    summary = {
        'analysis': analysis,
        'model': {
            'mesh_file': resolve_path(case['model']['mesh_file'], base_dir),
            'xml_file': resolve_path(case['model'].get('xml_file'), base_dir),
            'nelems': int(tw_model.nelems),
            'nholes': int(tw_model.nholes),
            'nregs': int(tw_model.nregs),
            'n_vcoils': int(tw_model.n_vcoils),
            'n_icoils': int(tw_model.n_icoils),
            'eta_surf': [float(val) for val in tw_model.get_eta_values(include_eta_vol=True)[0]],
        },
        'operators': ops,
        'results': results,
    }
    summary_file = os.path.join(out_dir, '{0}_{1}.json'.format(prefix, analysis))
    with open(summary_file, 'w') as fid:
        json.dump(summary, fid, indent=2)
    written.insert(0, summary_file)
    print('\nWrote:')
    for path in written:
        print('  '+path)
    return 0


if __name__ == '__main__':
    sys.exit(main())

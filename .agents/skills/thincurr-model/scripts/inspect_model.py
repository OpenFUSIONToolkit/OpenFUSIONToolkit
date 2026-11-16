#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-only
"""Inspect and verify a ThinCurr model (native mesh file plus optional XML).

Usage:
  python inspect_model.py MESH.h5 [--xml oft_in.xml] [--jumper-start N]
                          [--insulated-junctions 3 1,2] [--no-setup]
                          [--plot FILE.png] [--json FILE.json]

Reports mesh size, regions, topology (boundary loops, junction edges, connected
pieces), how holes/closures are defined, and where resistivity and coils come
from. Unless --no-setup is given, the model is also set up with ThinCurr, which
verifies the holes and closures and reports the number of unknowns, the resolved
per-region resistivity, and the coil sets.

Exits non-zero if the model cannot be used as-is (eg. missing resistivity or
failed hole verification).
"""
import argparse
import json
import os
import sys

import h5py
import numpy

# Above this many unknowns a dense L matrix is large (8*N^2 bytes); use HODLR
HODLR_SUGGEST_NELEMS = 10000


def _load_oft():
    """Make the OpenFUSIONToolkit package importable, honoring OFT_ROOTPATH."""
    root_path = os.getenv('OFT_ROOTPATH')
    if root_path is not None:
        sys.path.append(os.path.join(root_path, 'python'))


def parse_insulated(values):
    """Parse --insulated-junctions entries ("3" or "1,2") into setup_model() form.

    Args:
        values: List of strings from the command line (or None).

    Returns:
        List of region ints and region-pair tuples, or None.
    """
    if not values:
        return None
    out = []
    for value in values:
        parts = [int(part) for part in value.split(',')]
        if len(parts) == 1:
            out.append(parts[0])
        elif len(parts) == 2:
            out.append(tuple(parts))
        else:
            raise ValueError('Invalid insulated junction "{0}"'.format(value))
    return out


def read_file_info(mesh_file):
    """Read mesh arrays and ThinCurr metadata directly from the HDF5 file.

    Args:
        mesh_file: Native mesh file.

    Returns:
        Dictionary of raw mesh data and metadata.
    """
    info = {}
    with h5py.File(mesh_file, 'r') as fid:
        if 'mesh' not in fid:
            raise ValueError('"{0}" is not a native OFT mesh file (no "mesh" group)'.format(mesh_file))
        r = fid['mesh/R'][()]
        if r.shape[1] == 2:
            r = numpy.hstack((r, numpy.zeros((r.shape[0], 1))))
        info['r'] = r
        info['lc'] = fid['mesh/LC'][()] - 1
        info['reg'] = fid['mesh/REG'][()]
        if info['lc'].shape[1] != 3:
            raise ValueError('ThinCurr needs a 3-node triangular surface mesh, cells have {0} nodes'.format(info['lc'].shape[1]))
        # Legacy (Cubit) nodesets: used as holes when no "thincurr/holes" group exists
        nodesets = []
        if 'mesh/NUM_NODESETS' in fid:
            i = 1
            while 'mesh/NODESET{0:04d}'.format(i) in fid:
                nodesets.append(fid['mesh/NODESET{0:04d}'.format(i)][()] - 1)
                i += 1
        info['nodesets'] = nodesets
        tc = fid['thincurr'] if 'thincurr' in fid else {}
        info['holes'] = []
        if 'holes' in tc:
            nholes = tc['holes/NHOLES'][0]
            info['holes'] = [tc['holes/HOLE{0:04d}'.format(i+1)][()] - 1 for i in range(nholes)]
        info['holes_version'] = int(tc['HOLES_VERSION'][0]) if 'HOLES_VERSION' in tc else None
        info['closures'] = (tc['CLOSURES'][()] - 1).tolist() if 'CLOSURES' in tc else []
        info['junctions'] = 'junctions' in tc
        info['junction_holes'] = int(tc['junctions/NHOLES'][0]) if info['junctions'] else 0
        info['njumpers'] = int(tc['jumpers/NJUMPERS'][0]) if 'jumpers' in tc else 0
        info['eta_surf'] = tc['ETA_SURF'][()].tolist() if 'ETA_SURF' in tc else None
        info['eta_vol'] = tc['ETA_VOL'][()].tolist() if 'ETA_VOL' in tc else None
        info['thickness'] = tc['THICKNESS'][()].tolist() if 'THICKNESS' in tc else None
        info['file_coil_sets'] = int(tc['coils/NCOILSETS'][0]) if 'coils' in tc else 0
        info['nfp'] = int(tc['periodicity/NFP'][0]) if 'periodicity/NFP' in tc else None
        info['has_pmap'] = 'periodicity/PMAP' in tc
    return info


def mesh_statistics(info):
    """Compute region, edge-length, and topology statistics.

    Args:
        info: Dictionary from read_file_info().

    Returns:
        Dictionary of statistics.
    """
    r, lc, reg = info['r'], info['lc'], info['reg']
    # Cell areas
    v1 = r[lc[:, 1], :] - r[lc[:, 0], :]
    v2 = r[lc[:, 2], :] - r[lc[:, 0], :]
    areas = 0.5*numpy.linalg.norm(numpy.cross(v1, v2), axis=1)
    regions = {}
    for reg_id in numpy.unique(reg):
        mask = (reg == reg_id)
        regions[int(reg_id)] = {'cells': int(mask.sum()), 'area': float(areas[mask].sum())}
    # Unique edges and how many cells share each
    edges = numpy.sort(numpy.vstack((lc[:, [0, 1]], lc[:, [1, 2]], lc[:, [2, 0]])), axis=1)
    uedges, counts = numpy.unique(edges, axis=0, return_counts=True)
    lengths = numpy.linalg.norm(r[uedges[:, 1], :] - r[uedges[:, 0], :], axis=1)
    # Connected pieces (through shared vertices)
    import scipy.sparse
    from scipy.sparse.csgraph import connected_components
    used = numpy.unique(lc)
    graph = scipy.sparse.coo_matrix((numpy.ones(uedges.shape[0]), (uedges[:, 0], uedges[:, 1])),
                                    shape=(r.shape[0], r.shape[0]))
    _, labels = connected_components(graph, directed=False)
    pieces = numpy.unique(labels[used])
    npieces = pieces.shape[0]
    # Boundary loops: connected components of the boundary-edge graph
    bedges = uedges[counts == 1, :]
    nloops = 0
    loop_piece = numpy.zeros((0,), dtype=numpy.int32)
    if bedges.shape[0] > 0:
        bverts = numpy.unique(bedges)
        bgraph = scipy.sparse.coo_matrix((numpy.ones(bedges.shape[0]), (bedges[:, 0], bedges[:, 1])),
                                         shape=(r.shape[0], r.shape[0]))
        _, blabels = connected_components(bgraph, directed=False)
        loop_ids, first_vert = numpy.unique(blabels[bverts], return_index=True)
        nloops = loop_ids.shape[0]
        loop_piece = labels[bverts[first_vert]]
    # Expected holes/closures from the Euler characteristic of each piece (manifold meshes only):
    # chi = V - E + F = 2 - 2*genus - nboundary, holes = 2*genus + max(nboundary-1, 0), closure if closed
    expected = None
    if (counts > 2).sum() == 0:
        expected = {'holes': 0, 'closures': 0}
        cell_piece = labels[lc[:, 0]]
        edge_piece = labels[uedges[:, 0]]
        for piece in pieces:
            chi = (labels[used] == piece).sum() - (edge_piece == piece).sum() + (cell_piece == piece).sum()
            nbound = (loop_piece == piece).sum()
            if nbound > 0:
                expected['holes'] += int(1 - chi)
            else:
                expected['holes'] += int(2 - chi)
                expected['closures'] += 1
    return {
        'np': int(r.shape[0]),
        'nc': int(lc.shape[0]),
        'nregions': len(regions),
        'regions': regions,
        'bbox_min': r[used, :].min(axis=0).tolist(),
        'bbox_max': r[used, :].max(axis=0).tolist(),
        'edge_length': {'min': float(lengths.min()), 'mean': float(lengths.mean()), 'max': float(lengths.max())},
        'connected_pieces': int(npieces),
        'boundary_loops': int(nloops),
        'junction_edges': int((counts > 2).sum()),
        'unused_points': int(r.shape[0] - used.shape[0]),
        'expected_topology': expected,
    }


def setup_and_query(mesh_file, xml_file, jumper_start, insulated, nthreads):
    """Set up the model with ThinCurr (verifies holes) and query the result.

    Args:
        mesh_file: Native mesh file.
        xml_file: XML file or None.
        jumper_start: Index of the first legacy nodeset that is a jumper (0 for none).
        insulated: insulated_junctions argument or None.
        nthreads: OpenMP thread count.

    Returns:
        Dictionary of model properties.
    """
    _load_oft()
    from OpenFUSIONToolkit import OFT_env
    from OpenFUSIONToolkit.ThinCurr import ThinCurr
    myOFT = OFT_env(nthreads=nthreads)
    tw_model = ThinCurr(myOFT)
    tw_model.setup_model(mesh_file=mesh_file, xml_filename=xml_file, jumper_start=jumper_start,
                         insulated_junctions=insulated)
    sys.stdout.flush()
    eta_surf, eta_vol = tw_model.get_eta_values(include_eta_vol=True)
    thickness = tw_model.get_thickness()
    return {
        'nelems': int(tw_model.nelems),
        'np_active': int(tw_model.np_active),
        'nholes': int(tw_model.nholes),
        'nregs': int(tw_model.nregs),
        'eta_surf': [float(val) for val in eta_surf],
        'eta_vol': None if eta_vol is None else [float(val) for val in eta_vol],
        'thickness': [float(val) for val in thickness],
        'coils': {name.strip(): coil['type'] for name, coil in tw_model.coils.items()},
    }


def write_plot(info, filename):
    """Plot the mesh colored by region, with node-loop holes and nodesets overlaid.

    Args:
        info: Dictionary from read_file_info().
        filename: Output image path.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r, lc, reg = info['r'], info['lc'], info['reg']
    fig = plt.figure(figsize=(8, 7), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    cmap = plt.get_cmap('tab10')
    for i, reg_id in enumerate(numpy.unique(reg)):
        mask = (reg == reg_id)
        ax.plot_trisurf(r[:, 0], r[:, 1], r[:, 2], triangles=lc[mask, :], color=cmap(i % 10),
                        alpha=0.35, linewidth=0.1, edgecolor='k', label='region {0}'.format(reg_id))
    loops = [(hole, 'tab:red', '-') for hole in info['holes']]
    if not info['holes']:
        loops += [(nodeset, 'tab:blue', '--') for nodeset in info['nodesets'] if len(nodeset) > 1]
    for nodes, color, style in loops:
        nodes = numpy.append(nodes, nodes[0])
        ax.plot(r[nodes, 0], r[nodes, 1], r[nodes, 2], color=color, ls=style, lw=2)
    ax.set_box_aspect(numpy.maximum(numpy.ptp(r, axis=0), 1.E-3*numpy.ptp(r, axis=0).max()))
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    ax.legend(loc='upper left')
    ax.set_title('Regions (shaded), holes (red) / legacy nodesets (blue dashed)')
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Inspect and verify a ThinCurr model')
    parser.add_argument('mesh_file', help='Native ThinCurr mesh file (.h5)')
    parser.add_argument('--xml', default=None, help='XML file with <oft><thincurr> settings (resistivity, coils)')
    parser.add_argument('--jumper-start', type=int, default=0,
                        help='Index of the first legacy mesh nodeset that is a jumper, not a hole '
                             '(1-based, or negative to count from the end; default: 0, no jumpers)')
    parser.add_argument('--insulated-junctions', nargs='+', default=None,
                        help='Insulated junctions for split meshes: region ("3") or region pair ("1,2")')
    parser.add_argument('--no-setup', action='store_true', help='Only read the file; skip ThinCurr setup/verification')
    parser.add_argument('--plot', default=None, help='Write a PNG of the mesh, regions, and holes')
    parser.add_argument('--json', default=None, help='Write the report as JSON')
    parser.add_argument('--nthreads', type=int, default=2, help='OpenMP threads for setup (default: 2)')
    args = parser.parse_args()

    errors = []
    warnings = []
    info = read_file_info(args.mesh_file)
    stats = mesh_statistics(info)
    report = {'mesh_file': os.path.abspath(args.mesh_file),
              'xml_file': None if args.xml is None else os.path.abspath(args.xml),
              'mesh': stats}

    # Hole/closure definition source
    if info['junctions']:
        hole_source = 'junction-split mesh ({0} junction holes)'.format(info['junction_holes'])
    elif info['holes'] or info['closures'] or (info['holes_version'] is not None):
        hole_source = 'thincurr/holes ({0} holes, {1} closures{2})'.format(
            len(info['holes']), len(info['closures']),
            '' if info['holes_version'] is None else ', OFT_ThinCurr_holes v{0}'.format(info['holes_version']))
    elif info['nodesets']:
        hole_source = 'legacy mesh nodesets ({0})'.format(len(info['nodesets']))
    else:
        hole_source = 'none'
    report['hole_source'] = hole_source
    skip_setup = False
    expected = stats['expected_topology']
    if expected is not None and not info['junctions']:
        if hole_source == 'none':
            if expected['holes'] > 0 or expected['closures'] > 0:
                errors.append('Mesh topology needs {holes} hole(s) and {closures} closure(s) but none are defined; '
                              'run OFT_ThinCurr_holes (compute_homology) on it'.format(**expected))
                skip_setup = True
        elif hole_source.startswith('legacy'):
            nsets = len(info['nodesets'])
            if args.jumper_start == 0 and nsets > expected['holes']:
                errors.append('{0} legacy nodesets but the topology needs only {1} hole(s); the extra nodesets are '
                              'probably jumpers (try --jumper-start {2}, and "jumper_start" when running), or '
                              'recreate the holes with OFT_ThinCurr_holes'.format(nsets, expected['holes'],
                                                                                  expected['holes'] - nsets))
                skip_setup = True
        elif len(info['holes']) != expected['holes'] or len(info['closures']) != expected['closures']:
            warnings.append('File defines {0} hole(s) and {1} closure(s); the topology suggests {2} and {3}'.format(
                len(info['holes']), len(info['closures']), expected['holes'], expected['closures']))
    if stats['junction_edges'] > 0 and not info['junctions']:
        errors.append('{0} edges are shared by more than two cells (junctions) but the mesh is not split; '
                      'run OFT_ThinCurr_holes on it'.format(stats['junction_edges']))
    if stats['unused_points'] > 0:
        warnings.append('{0} mesh points are not used by any cell'.format(stats['unused_points']))

    # Resistivity and coil sources (mesh file takes precedence over XML)
    xml_text = ''
    if args.xml is not None:
        with open(args.xml, 'r') as fid:
            xml_text = fid.read()
    xml_has_eta = any(tag in xml_text for tag in ('<eta>', '<eta_surf>', '<eta_vol>'))
    xml_has_coils = any(tag in xml_text for tag in ('<icoils', '<vcoils'))
    file_has_eta = (info['eta_surf'] is not None) or (info['eta_vol'] is not None)
    if file_has_eta and xml_has_eta:
        warnings.append('Resistivity is set in both the mesh file and the XML; the mesh file values are used')
    if (info['file_coil_sets'] > 0) and xml_has_coils:
        warnings.append('Coils are defined in both the mesh file and the XML; only the mesh file coils are used')
    report['eta_source'] = 'mesh file' if file_has_eta else ('xml' if xml_has_eta else 'none')
    report['coil_source'] = 'mesh file' if info['file_coil_sets'] > 0 else ('xml' if xml_has_coils else 'none')
    if info['nfp'] is not None or info['has_pmap']:
        report['periodic'] = {'nfp': info['nfp'], 'pmap': info['has_pmap']}

    # Print file-level summary
    print('ThinCurr model: {0}'.format(args.mesh_file))
    print('  # of points / cells   = {0} / {1}'.format(stats['np'], stats['nc']))
    print('  Bounding box [m]      = {0} to {1}'.format(
        numpy.round(stats['bbox_min'], 4).tolist(), numpy.round(stats['bbox_max'], 4).tolist()))
    print('  Edge length [m]       = min {min:.3E}, mean {mean:.3E}, max {max:.3E}'.format(**stats['edge_length']))
    print('  Connected pieces      = {0}'.format(stats['connected_pieces']))
    print('  Boundary loops        = {0}'.format(stats['boundary_loops']))
    print('  Junction edges        = {0}{1}'.format(stats['junction_edges'], ' (split mesh)' if info['junctions'] else ''))
    if expected is not None:
        print('  Topology needs        = {holes} hole(s), {closures} closure(s)'.format(**expected))
    print('  Holes/closures from   = {0}'.format(hole_source))
    print('  Resistivity from      = {0}'.format(report['eta_source']))
    print('  Coils from            = {0}'.format(report['coil_source']))
    if 'periodic' in report:
        print('  Periodic              = {0}'.format(report['periodic']))
    print('  Regions:')
    for reg_id, reg_info in stats['regions'].items():
        print('    {0:4d}: {1:8d} cells, area {2:.4E} m^2'.format(reg_id, reg_info['cells'], reg_info['area']))

    if skip_setup:
        print('\nSkipping ThinCurr setup: it would abort on the problem reported below')
    if not (args.no_setup or skip_setup):
        print('\n--- ThinCurr setup (verifies holes and closures) ---', flush=True)
        try:
            model = setup_and_query(args.mesh_file, args.xml, args.jumper_start,
                                    parse_insulated(args.insulated_junctions), args.nthreads)
        except Exception as setup_error:
            errors.append('ThinCurr setup failed: {0}'.format(setup_error))
            model = None
        if model is not None:
            report['model'] = model
            print('\nModel summary:')
            print('  # of unknowns (nelems) = {0}'.format(model['nelems']))
            print('  # of holes             = {0}'.format(model['nholes']))
            print('  Coils                  = {0}'.format(model['coils'] if model['coils'] else 'none'))
            print('  Region  eta_surf [Ohm]  thickness [m]')
            for i, eta in enumerate(model['eta_surf']):
                thick = model['thickness'][i]
                print('    {0:4d}  {1:12.4E}  {2}'.format(i+1, eta, '{0:.4E}'.format(thick) if thick > 0.0 else '-'))
            if any(eta <= 0.0 for eta in model['eta_surf']):
                errors.append('Surface resistivity is missing or non-positive for at least one region; '
                              'set <eta> (one value per region) in the XML or ETA_SURF in the mesh file')
            if model['nelems'] > HODLR_SUGGEST_NELEMS:
                warnings.append('{0} unknowns: a dense L matrix needs {1:.1f} GB; use HODLR when running'.format(
                    model['nelems'], 8.0*model['nelems']**2/1.E9))

    if args.plot is not None:
        write_plot(info, args.plot)
        print('\nWrote plot: {0}'.format(args.plot))

    report['warnings'] = warnings
    report['errors'] = errors
    if args.json is not None:
        with open(args.json, 'w') as fid:
            json.dump(report, fid, indent=2)
        print('Wrote report: {0}'.format(args.json))
    for warning in warnings:
        print('WARNING: {0}'.format(warning))
    for error in errors:
        print('ERROR: {0}'.format(error), file=sys.stderr)
    if errors:
        return 1
    print('\nModel OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())

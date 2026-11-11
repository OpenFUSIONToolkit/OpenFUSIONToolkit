#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python cli for ThinCurr hole construction using homology computations

See @ref thincurr_homology_tool

@authors Chris Hansen
@date July 2026
@ingroup doxy_oft_python
'''
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from .meshing import read_ThinCurr_mesh, write_ThinCurr_mesh, HOLES_VERSION
from . import _tree_homology
from . import _manifold_homology
from .junctions import has_junctions, find_junction_holes, verify_junction_mesh, split_jumpers


def _verify_counts(input_model, nholes, nclosures):
    '''! Check the number of holes and closures in a mesh file against the expected counts

    @param input_model Mesh information from @ref meshing.read_ThinCurr_mesh
    @param nholes Expected number of holes
    @param nclosures Expected number of closures
    '''
    # Check holes
    if 'holes' in input_model['thincurr']:
        if len(input_model['thincurr']['holes']) != nholes:
            raise ValueError("Incorrect number of holes found {1} in file (expected {0})".format(nholes, len(input_model['thincurr']['holes'])))
    else:
        if 'nodesets' in input_model:
            if len(input_model['nodesets']) >= nholes:
                print("""
  Note: Input file appears to be an older ThinCurr mesh file.
    Found {0} nodesets, expected {1} hole(s).""".format(len(input_model['nodesets']), nholes))
                if len(input_model['nodesets']) > nholes:
                    print("More nodesets than holes may indicate the presence of jumpers instead of an issue.")
        else:
            if nholes > 0:
                raise ValueError("Incorrect number of holes found 0 in file (expected {0})".format(nholes))
    # Check closures
    if 'closures' in input_model['thincurr']:
        if len(input_model['thincurr']['closures']) != nclosures:
            raise ValueError("Incorrect number of closures found {1} in file (expected {0})".format(nclosures, len(input_model['thincurr']['closures'])))
    else:
        if 'sidesets' in input_model:
            if len(input_model['sidesets'][0]) == nclosures:
                print("""
  Note: Input file appears to be an older ThinCurr mesh file.
    Found {0} entries in first sideset, expected {1} closure(s).""".format(len(input_model['sidesets'][0]), nclosures))
            else:
                raise ValueError("No closures found and incorrect number of sidesets found {1} in file (expected {0})".format(nclosures, len(input_model['sidesets'][0])))
        else:
            if nclosures > 0:
                raise ValueError("Incorrect number of closures found 0 in file (expected {0})".format(nclosures))


def _write_holes(out_file, input_model, holes, closures):
    '''! Copy a mesh to a new file with the given holes and closures

    @param out_file Output mesh file
    @param input_model Mesh information from @ref meshing.read_ThinCurr_mesh
    @param holes List of hole loops (0-based vertex indices)
    @param closures List of closure cells (0-based)
    '''
    write_ThinCurr_mesh(out_file, input_model['r'], input_model['lc'], input_model['reg'],
                        holes=holes,
                        jumpers=input_model['thincurr'].get('jumpers'), closures=closures,
                        eta_surf=input_model['thincurr'].get('eta_surf'),
                        eta_vol=input_model['thincurr'].get('eta_vol'),
                        thickness=input_model['thincurr'].get('thickness'),
                        coil_sets=input_model['thincurr'].get('coil_sets'),
                        reg_names=input_model.get('reg_names'), reg_attrs=input_model.get('reg_attrs'),
                        holes_version=HOLES_VERSION)


def _compute_junction_holes(input_model, out_file, verify_only, print_info, split_regions):
    '''! Split a mesh with junctions (edges shared by more than two cells) and choose its holes

    See @ref junctions.find_junction_holes. The output file holds the split mesh and, in the
    `thincurr/junctions` group, the junction table, explicit holes, vertex map and P.

    @param input_model Mesh information from @ref meshing.read_ThinCurr_mesh
    @param out_file Output mesh file
    @param verify_only Only report the junction model, do not write output file?
    @param print_info Print information about the junction model?
    @param split_regions Also split the mesh at all region interfaces?
    '''
    out = find_junction_holes(input_model['r'], input_model['lc'], input_model['reg'], split_regions=split_regions)
    if print_info:
        print("Junction mesh:")
        print("  Split at region interfaces = {0}".format(split_regions))
        print("  # of junction lines = {0}".format(len(out['lines'])))
        print("  # of junction edges = {0}".format(len(out['kje'])-1))
        print("  # of vertices after split = {0}".format(out['r'].shape[0]))
        print("  # of holes = {0} (sheet boundaries and cuts)".format(len(out['holes'])))
        print("  # of cuts = {0}".format(len(out['cuts'])))
        print("  # of closures = {0}".format(len(out['info']['closures'])))
        print("  # of unknowns = {0} (independent currents {1})".format(out['n_unknowns'], out['n_loops']))
    if verify_only:
        return
    thincurr = input_model['thincurr']
    jumpers = split_jumpers(input_model['lc'], out['lc'], out['vmap'], out['jedges'], thincurr.get('jumpers', []))
    write_ThinCurr_mesh(out_file, out['r'], out['lc'], input_model['reg'],
                        junctions=out, jumpers=jumpers,
                        eta_surf=thincurr.get('eta_surf'), eta_vol=thincurr.get('eta_vol'),
                        thickness=thincurr.get('thickness'), coil_sets=thincurr.get('coil_sets'),
                        reg_names=input_model.get('reg_names'), reg_attrs=input_model.get('reg_attrs'),
                        holes_version=HOLES_VERSION)


def compute_homology(in_file, out_file=None, plot_final=False, plot_steps=False, show_omitted=False,
                     show_covering=False, debug=False, ref_point=None, optimize_holes=False, verify_only=False,
                     use_legacy=False, verbose=False, split_regions=False):
    '''! Compute holes and closures for ThinCurr meshes

    By default holes are selected using a tree-cotree approach (see @ref _tree_homology). The previous Greedy
    Homotopy approach (see @ref _manifold_homology) is available with `use_legacy=True`.

    @param in_file Input mesh file
    @param out_file Ouput mesh file (default: Input file name with "-homology" appended)
    @param plot_final Show final homology basis?
    @param plot_steps Show intermediate bases for each distinct surface (legacy method only)?
    @param show_omitted Show boundary cycles that are omitted?
    @param show_covering Show covering triangles for boundary cycles (legacy method only, only used if `plot_steps`)?
    @param debug Print additional debug information?
    @param ref_point Reference location for base point [x,y,z] (default: [0,0,0])
    @param optimize_holes Sample additional points to attempt to optimize holes (little effect with the default method)?
    @param verify_only Verify the correct number of holes and closures are found, but do not write output file?
    @param use_legacy Use the legacy Greedy Homotopy method?
    @param verbose Print progress for each search base point (tree-cotree method only)?
    @param split_regions Also split the mesh at all region interfaces (only needed to insulate regions, see @ref thincurr_holes_tool)?
    '''
    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-homology.h5"

    if ref_point is not None:
        ref_point = np.asarray(ref_point)
        if ref_point.shape[0] != 3:
            raise ValueError('`ref_point` must be a 3-vector')
    else:
        ref_point = np.r_[0.0,0.0,0.0]

    if optimize_holes and verify_only:
        raise ValueError('Optimizing holes not necessary or permitted when `verify_only` is True')
    if (plot_steps or show_covering) and (not use_legacy):
        print('Note: "plot_steps" and "show_covering" are only used with the legacy method')
    print_info = debug or (not verify_only)

    # Load mesh
    input_model = read_ThinCurr_mesh(in_file)
    if 'thincurr' in input_model:
        if 'nfp' in input_model['thincurr']:
            if split_regions or has_junctions(input_model['lc'], input_model['reg']):
                raise ValueError('Junctions are not supported for periodic ThinCurr meshes (nfp > 0) at this time')
            if verify_only:
                print("Periodic ThinCurr model detected (nfp > 0), skipping homology calculation")
                return
            else:
                raise ValueError('Homology calculation does not support periodic ThinCurr meshes (nfp > 0) at this time')
    else:
        input_model['thincurr'] = {}
    split_mesh = 'junctions' in input_model['thincurr']
    unsplit_junctions = (not split_mesh) and has_junctions(input_model['lc'], input_model['reg'])
    if use_legacy and (split_mesh or unsplit_junctions or split_regions):
        raise ValueError('The legacy method does not support meshes with junctions (edges shared by more than two cells, '
                         'or sheets touching at a vertex) or "split_regions"')
    if split_mesh:
        if verify_only:
            if split_regions and (not input_model['thincurr']['junctions']['split_regions']):
                print('Note: "split_regions" is ignored when verifying a mesh that is already split; this mesh is not '
                      'split at region interfaces (rerun "OFT_ThinCurr_holes --split_regions" on the unsplit mesh to '
                      'split it there)')
            njedges, nunknowns = verify_junction_mesh(input_model['r'], input_model['lc'], input_model['reg'],
                                                      input_model['thincurr']['junctions'])
            if print_info:
                print("Junction mesh: {0} junction edges, {1} unknowns".format(njedges, nunknowns))
            return
        raise ValueError('Mesh has already been split at junctions')
    if unsplit_junctions or split_regions:
        if unsplit_junctions and verify_only:
            raise ValueError('Mesh has junctions (edges shared by more than two cells, or sheets touching at a vertex) '
                             'but no junction data, run "OFT_ThinCurr_holes" to split it')
        return _compute_junction_holes(input_model, out_file, verify_only, print_info, split_regions)

    # Select (or count) holes and closures
    if use_legacy:
        out = _manifold_homology.find_holes(input_model['r'].copy(), input_model['lc'].copy(), ref_point=ref_point,
                                            optimize_holes=optimize_holes, print_info=print_info, debug=debug,
                                            plot_steps=plot_steps, show_covering=show_covering, show_omitted=show_omitted)
        nholes = len(out['boundary_holes']) + len(out['internal_holes'])
        nclosures = len(out['closures'])
    elif verify_only:
        nholes, nclosures = _tree_homology.count_holes(input_model['r'], input_model['lc'])
    else:
        out = _tree_homology.find_holes(input_model['r'], input_model['lc'], ref_point=ref_point, optimize=optimize_holes, info=print_info,
                                        verbose=verbose)
        nholes = len(out['boundary_holes']) + len(out['internal_holes'])
        nclosures = len(out['closures'])
    if print_info:
        print()
        print("Final model:")
        print("    # of holes = {0}".format(nholes))
        print("    # of closures = {0}".format(nclosures))

    # Verify if requested
    if verify_only:
        _verify_counts(input_model, nholes, nclosures)
        return

    # Copy mesh to new file and replace holes/closures
    _write_holes(out_file, input_model, out['boundary_holes'] + out['internal_holes'], out['closures'])

    # Plot final cycles
    if plot_final:
        vertex_full = input_model['r']
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_trisurf(vertex_full[:,0], vertex_full[:,1], vertex_full[:,2], triangles=input_model['lc'], color=[0.0, 0.0, 0.0, 0.1])
        for cycle in out['boundary_holes']:
            cycle = np.append(cycle, cycle[0])
            ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:blue')
        if show_omitted:
            for cycle in out['omitted']:
                cycle = np.append(cycle, cycle[0])
                ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='k')
        for cycle in out['internal_holes']:
            cycle = np.append(cycle, cycle[0])
            ax.plot(vertex_full[cycle,0], vertex_full[cycle,1], vertex_full[cycle,2], color='tab:orange')
        for closure_cell in out['closures']:
            closure = input_model['lc'][closure_cell,0]
            ax.plot(vertex_full[closure,0], vertex_full[closure,1], vertex_full[closure,2], 'o', color='k')
        ax.set_aspect('equal','box')
        plt.show()


def script_entry():
    '''! Command line interface for computing holes and closures for ThinCurr meshes

    See @ref thincurr_holes_tool
    '''
    parser = argparse.ArgumentParser()
    parser.description = "Compute holes and closures for ThinCurr meshes"
    parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    parser.add_argument("--out_file", type=str, default=None, help='Ouput mesh file (default: Input file name with "-homology" appended)')
    parser.add_argument("--keep_nodeset_start", type=int, default=None, help="REMOVED")
    parser.add_argument("--plot_final", action="store_true", default=False, help="Show final homology basis?")
    parser.add_argument("--plot_steps", action="store_true", default=False, help="Show intermediate bases for each distinct surface (legacy method only)?")
    parser.add_argument("--show_omitted", action="store_true", default=False, help="Show boundary cycles that are omitted?")
    parser.add_argument("--show_covering", action="store_true", default=False, help="Show covering triangles for boundary cycles (legacy method only, only used if `--plot_steps`)?")
    parser.add_argument("--debug", action="store_true", default=False, help="Print additional debug information?")
    parser.add_argument("--ref_point", default=None, type=float, nargs='+', help='Reference location for base point [x,y,z] (default: [0,0,0])')
    parser.add_argument("--optimize_holes", action="store_true", default=False, help="Sample additional points to attempt to optimize holes (little effect with the default method)?")
    parser.add_argument("--use_legacy", action="store_true", default=False, help="Use the legacy Greedy Homotopy method?")
    parser.add_argument("--verbose", action="store_true", default=False, help="Print progress for each search base point (tree-cotree method only)?")
    parser.add_argument("--split_regions", action="store_true", default=False, help="Also split at all region interfaces (only needed to insulate regions)?")
    options = parser.parse_args()

    if options.keep_nodeset_start is not None:
        parser.exit(-1, '"--keep_nodeset_start" has been removed, please use `ThinCurr_mesh_tool.py` to convert non-hole nodsets before running this script.')

    compute_homology(
        in_file=options.in_file,
        out_file=options.out_file,
        plot_final=options.plot_final,
        plot_steps=options.plot_steps,
        show_omitted=options.show_omitted,
        show_covering=options.show_covering,
        debug=options.debug,
        ref_point=options.ref_point,
        optimize_holes=options.optimize_holes,
        use_legacy=options.use_legacy,
        verbose=options.verbose,
        split_regions=options.split_regions
    )


## @page thincurr_holes_tool `OFT_ThinCurr_holes`: ThinCurr hole/closure script
#
# @section OFT_ThinCurr_holes_desc Description and options
# This script identifies holes and closures in a ThinCurr surface mesh using
# homology computations. By default holes are selected using a tree-cotree approach
# (see \ref OpenFUSIONToolkit.ThinCurr._tree_homology); the previous Greedy Homotopy approach
# is available using `--use_legacy`.
#
# If the mesh has junctions (edges shared by more than two cells, or sheets of different regions touching only at a
# point), or `--split_regions` is given, it is split at the junctions instead (not supported for periodic meshes or
# with `--use_legacy`) and the output file holds the split mesh with the junction data (`thincurr/junctions`): see
# \ref doc_tw_main_junctions. Each region must be a separate sheet. See \ref OFT_ThinCurr_holes_split for when
# to use `--split_regions`.
#
#```shell
# usage: OFT_ThinCurr_holes [-h] --in_file IN_FILE [--out_file OUT_FILE] [--keep_nodeset_start KEEP_NODESET_START] [--plot_final] [--plot_steps] [--show_omitted]
#                                 [--show_covering] [--debug] [--ref_point REF_POINT [REF_POINT ...]] [--optimize_holes] [--use_legacy] [--verbose] [--split_regions]
#
# options:
#     -h, --help            show this help message and exit
#     --in_file IN_FILE     Input mesh file
#     --out_file OUT_FILE   Ouput mesh file (default: Input file name with "-homology" appended)
#     --keep_nodeset_start KEEP_NODESET_START
#                           REMOVED
#     --plot_final          Show final homology basis?
#     --plot_steps          Show intermediate bases for each distinct surface (legacy method only)?
#     --show_omitted        Show boundary cycles that are omitted?
#     --show_covering       Show covering triangles for boundary cycles (legacy method only, only used if `--plot_steps`)?
#     --debug               Print additional debug information?
#     --ref_point REF_POINT [REF_POINT ...]
#                           Reference location for base point [x,y,z] (default: [0,0,0])
#     --optimize_holes      Sample additional points to attempt to optimize holes (little effect with the default method)?
#     --use_legacy          Use the legacy Greedy Homotopy method?
#     --verbose             Print progress for each search base point (tree-cotree method only)?
#     --split_regions       Also split at all region interfaces (only needed to insulate regions)?
#```
#
# @section OFT_ThinCurr_holes_split Splitting at region interfaces (`--split_regions`)
# By default a mesh is only split at its junction edges. With `--split_regions` it is also split at every edge where
# two regions meet, with or without junctions, and those edges become two-cell junction edges.
#
# On its own this does not change the model: with no insulation, current conservation across a two-cell junction edge
# is ordinary continuity, so the sheets are joined back together and the model carries exactly the same currents as
# without the split (same number of unknowns, same resistance and inductance up to quadrature-level differences in the
# inductance). The split only matters for insulation, which acts at junction edges only: it allows sheets that meet
# two at a time (e.g. the sectors of a torus) to be insulated from each other when the model is set up, with
# `setup_model(insulated_junctions=...)`, without generating a new mesh file. Without the split, regions or pairs of
# regions that do not meet at a junction edge cannot be insulated (an error is raised if this is requested).
#
# Only use `--split_regions` if this insulation is needed, as the split mesh has more vertices, junction constraints
# and hole elements, which increases the cost of setting up the model, and jumpers cannot run along the split region
# interfaces. The choice is stored in the output file and used when the mesh is verified or insulated.

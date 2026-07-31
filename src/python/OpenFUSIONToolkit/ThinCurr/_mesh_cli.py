#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''!Utility for manipulating and combining ThinCurr surface meshes.

This script operates on 3-node (linear) triangular surface meshes stored in the
Open FUSION Toolkit "native" HDF5 mesh format. Mesh I/O is delegated to
@ref OpenFUSIONToolkit.ThinCurr.meshing.read_ThinCurr_mesh "read_ThinCurr_mesh"
and @ref OpenFUSIONToolkit.ThinCurr.meshing.write_ThinCurr_mesh "write_ThinCurr_mesh",
so ThinCurr properties (holes, closures, conductivity, thickness, periodicity) are
read and written through the standard library routines. It supports two top-level
workflows via sub-commands:

Combine two (or more) existing mesh files into a single mesh:

    python ThinCurr_mesh_tool.py combine --in_files a.h5 b.h5 --out_file combined.h5

Modify a single mesh file by removing regions, applying a rigid transform, and/or
generating multiple shifted/rotated copies:

    # Remove regions 2 and 3
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --remove_regions 2 3

    # Shift the whole mesh by (0,0,0.5)
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --shift 0 0 0.5

    # Stretch the mesh by 2x along x (scale about the origin)
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --scale 2 1 1

    # Build a 6-fold toroidal array (60 deg spacing about z):
    # the original plus 5 rotated copies at 60,120,180,240,300 deg
    python ThinCurr_mesh_tool.py modify --in_file segment.h5 --copies 5 --rotate z 60

    # Set a uniform surface resistivity on all regions
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --eta_surf 1.257E-5

    # Set per-region volumetric resistivity with thickness (2 regions)
    python ThinCurr_mesh_tool.py modify --in_file a.h5 \
        --eta_vol 1.257E-5 2.5E-5 --thickness 1.0E-3 2.0E-3

    # Read resistivity from an <oft><thincurr> XML input file
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --eta_from_xml oft_in.xml

    # Treat the last two nodesets as jumpers; all others become holes
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --jumpers -1 -2

    # Treat nodesets 2 and 3 as jumpers (half-open range [2:4))
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --jumper_range 2 4

    # Treat nodesets 3 onward as jumpers (STOP omitted -> to the end)
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --jumper_range 3

    # Mark every nodeset as a jumper
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --jumper_range 0

By default the modify workflow treats every native NODESET loop as a hole and
writes them as ThinCurr holes; --jumpers/--jumper_range reclassify selected
nodesets (0-based, negative indexing allowed) as jumpers instead. The combine
workflow aborts if any input still contains native NODESET entries.

Only one of --shift, --rotate, or --scale may be given at a time. When --copies
is given, --shift or --rotate defines the per-copy increment (--scale is not
allowed); otherwise the transform is applied once to the whole mesh.

@note Only linear 3-node triangular meshes are supported. High-order meshes or
non-triangular cells will be rejected. ThinCurr properties (conductivity,
thickness, periodicity) are preserved through geometry-only transforms but are
dropped by structural operations (region removal, combination, replication).

@authors Open FUSION Toolkit contributors
@date July 2026
'''
import os
import re
import xml.etree.ElementTree as ET
import numpy as np
from OpenFUSIONToolkit.ThinCurr.meshing import read_ThinCurr_mesh, write_ThinCurr_mesh
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_XML


class ThinCurrMesh:
    '''!Container for a 3-node triangular ThinCurr surface mesh.

    Connectivity is stored 0-based, matching the read_ThinCurr_mesh /
    write_ThinCurr_mesh library routines, which handle the 1-based on-disk
    conversion internally.

    Nodeset loops (holes and jumpers) are held in a single ordered `nodesets`
    list with a parallel boolean `is_jumper` list recording each one's role;
    on save they are split into ThinCurr holes and jumpers accordingly.

    Attributes:
      r          Vertex coordinate list, shape [np,3] (float64); a 2D [np,2]
                 list supplied at construction is upgraded to 3D with Z=0
      lc         Cell (triangle) vertex list, shape [nc,3], 0-based (int32)
      reg        Per-cell region index, shape [nc], values in 1..nregions (int32)
      nodesets   List of node-index arrays (0-based) for hole/jumper loops
      is_jumper  Parallel bool list; True marks a nodeset as a jumper (else hole)
      sidesets   List of cell-index arrays (0-based); ThinCurr "closures"
      props      Dict of ThinCurr properties carried through unchanged
                 (eta_surf, eta_vol, thickness, pmap, nfp)
    '''

    #: ThinCurr properties preserved verbatim on a geometry-only round-trip
    _PROP_KEYS = ('eta_surf', 'eta_vol', 'thickness') # Ignore 'pmap' and 'nfp' for now (only limited uses)
    _MESH_PROP_KEYS = ('reg_names', 'reg_attrs')

    def __init__(self, r, lc, reg, nodesets=None, sidesets=None, props=None,
                 mesh_props=None, is_jumper=None, had_native_nodesets=False, coil_sets=None):
        '''!Construct a mesh from its component arrays

        @param r Vertex list [np,3]
        @param lc Triangle vertex list [nc,3], 0-based
        @param reg Per-cell region index [nc]
        @param nodesets List of 0-based node-index arrays (optional)
        @param sidesets List of 0-based cell-index arrays (optional)
        @param props Dict of ThinCurr properties to carry through (optional)
        @param mesh_props Dict of base mesh properties to carry through (optional)
        @param is_jumper Parallel bool list marking jumper nodesets (optional)
        @param had_native_nodesets Whether the source file had native NODESETs
        @param coil_sets List of ThinCurr coil sets (optional)
        '''
        self.r = np.ascontiguousarray(r, dtype=np.float64)
        self.lc = np.ascontiguousarray(lc, dtype=np.int32)
        self.reg = np.ascontiguousarray(reg, dtype=np.int32).reshape(-1)
        if self.lc.ndim != 2 or self.lc.shape[1] != 3:
            raise ValueError("Cell list must have shape [nc,3] (3-node triangles)")
        if self.r.ndim != 2 or self.r.shape[1] not in (2, 3):
            raise ValueError("Vertex list must have shape [np,2] or [np,3]")
        if self.r.shape[1] == 2:  # Upgrade 2D point list to 3D (Z=0)
            self.r = np.hstack((self.r, np.zeros((self.r.shape[0], 1))))
        if self.reg.shape[0] != self.lc.shape[0]:
            raise ValueError("Region list length must match number of cells")
        self.nodesets = [np.asarray(ns, dtype=np.int32).reshape(-1) for ns in (nodesets or [])]
        self.sidesets = [np.asarray(ss, dtype=np.int32).reshape(-1) for ss in (sidesets or [])]
        if is_jumper is None:
            self.is_jumper = [False] * len(self.nodesets)
        else:
            self.is_jumper = [bool(x) for x in is_jumper]
            if len(self.is_jumper) != len(self.nodesets):
                raise ValueError("is_jumper length must match number of nodesets")
        self.props = dict(props) if props else {}
        self.mesh_props = dict(mesh_props) if mesh_props else {}
        self._had_native_nodesets = had_native_nodesets
        self.coil_sets = coil_sets

    @property
    def np(self):
        '''!Number of vertices'''
        return self.r.shape[0]

    @property
    def nc(self):
        '''!Number of cells (triangles)'''
        return self.lc.shape[0]

    @property
    def nregions(self):
        '''!Number of regions (maximum region index)'''
        if self.nc == 0:
            return 0
        return int(self.reg.max())

    @property
    def nholes(self):
        '''!Number of hole nodesets'''
        return sum(1 for j in self.is_jumper if not j)

    @property
    def njumpers(self):
        '''!Number of jumper nodesets'''
        return sum(1 for j in self.is_jumper if j)

    def copy(self):
        '''!Return a deep copy of the mesh'''
        new = ThinCurrMesh(
            self.r.copy(), self.lc.copy(), self.reg.copy(),
            [ns.copy() for ns in self.nodesets],
            [ss.copy() for ss in self.sidesets],
            {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.props.items()},
            mesh_props={k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.mesh_props.items()},
            is_jumper=list(self.is_jumper), had_native_nodesets=self._had_native_nodesets, coil_sets=self.coil_sets)
        return new

    # ------------------------------------------------------------------ I/O
    @classmethod
    def load(cls, filename):
        '''!Load a mesh using @ref OpenFUSIONToolkit.ThinCurr.meshing.read_ThinCurr_mesh

        Native ``NODESETXXXX`` loops are read (defaulting to holes) along with any
        ThinCurr holes/jumpers already stored in the file; ThinCurr closures or
        native side sets are read as closures. All indices are already 0-based.

        @param filename Path to input mesh file
        @result New @ref ThinCurrMesh instance
        '''
        print()
        print("Reading mesh: {0}".format(filename))
        mesh_info = read_ThinCurr_mesh(filename)
        if 'ho_info' in mesh_info:
            raise ValueError("High-order meshes are not supported (found high-order node info)")
        mesh_type = mesh_info.get('type', 'tri')
        if mesh_type != 'tri':
            raise ValueError("Only 3-node triangular ('tri') meshes are supported "
                             "(found '{0}')".format(mesh_type))
        r = np.asarray(mesh_info.pop('r'))
        lc = np.asarray(mesh_info.pop('lc'))  # already 0-based
        reg = np.asarray(mesh_info.pop('reg'))
        if lc.ndim != 2 or lc.shape[1] != 3:
            raise ValueError("Only 3-node triangular meshes are supported "
                             "(found {0} nodes/cell)".format(lc.shape[1]))
        if r.ndim == 2 and r.shape[1] == 2:
            print("  Note: input mesh is 2D; padding point list to 3D (Z=0)")
        tc = mesh_info.get('thincurr', {})

        def _nonempty(seq):
            return [np.asarray(a) for a in seq if np.asarray(a).shape[0] > 0]

        # Ordered loop list: native NODESETs (holes by default), then any stored
        # ThinCurr holes, then any stored ThinCurr jumpers.
        native = _nonempty(mesh_info.get('nodesets', []))
        stored_holes = _nonempty(tc.get('holes', []))
        stored_jumpers = _nonempty(tc.get('jumpers', []))
        nodesets = native + stored_holes + stored_jumpers
        is_jumper = ([False] * (len(native) + len(stored_holes))) + ([True] * len(stored_jumpers))
        # Closures: prefer ThinCurr closures, else native side sets
        if 'closures' in tc and np.asarray(tc['closures']).shape[0] > 0:
            sidesets = [np.asarray(tc['closures'])]
        else:
            sidesets = _nonempty(mesh_info.get('sidesets', []))
        props = {k: tc[k] for k in cls._PROP_KEYS if k in tc}
        mesh_props = {k: mesh_info[k] for k in cls._MESH_PROP_KEYS if k in mesh_info}
        mesh = cls(r, lc, reg, nodesets, sidesets, props,
                   coil_sets=tc.get('coil_sets'),
                   mesh_props=mesh_props,
                   is_jumper=is_jumper, had_native_nodesets=len(native) > 0)
        mesh.print_info()
        return mesh

    def save(self, filename):
        '''!Save the mesh using @ref OpenFUSIONToolkit.ThinCurr.meshing.write_ThinCurr_mesh

        Nodesets are split into ThinCurr holes and jumpers per `is_jumper`, and
        side sets are merged into a single ThinCurr closure set. All indices are
        passed 0-based (the library converts to 1-based for storage).

        @param filename Path to output mesh file
        '''
        holes = [ns for ns, j in zip(self.nodesets, self.is_jumper) if not j]
        jumpers = [ns for ns, j in zip(self.nodesets, self.is_jumper) if j]
        if len(self.sidesets) == 0:
            closures = None
        elif len(self.sidesets) == 1:
            closures = self.sidesets[0]
        else:
            print("  Note: merging {0} side sets into a single closure set".format(len(self.sidesets)))
            closures = np.concatenate(self.sidesets)
        write_ThinCurr_mesh(filename, self.r, self.lc, self.reg,
                            reg_attrs=self.mesh_props.get('reg_attrs'),
                            reg_names=self.mesh_props.get('reg_names'),
                            holes=holes, jumpers=jumpers, closures=closures,
                            eta_surf=self.props.get('eta_surf'),
                            eta_vol=self.props.get('eta_vol'),
                            thickness=self.props.get('thickness'),
                            coil_sets=self.coil_sets,
                            pmap=self.props.get('pmap'),
                            nfp=self.props.get('nfp'))
        self.print_info()

    def print_info(self):
        '''!Print a short summary of the mesh contents'''
        print("  # of points   = {0}".format(self.np))
        print("  # of cells    = {0}".format(self.nc))
        print("  # of regions  = {0}".format(self.nregions))
        print("  # of holes    = {0}".format(self.nholes))
        print("  # of jumpers  = {0}".format(self.njumpers))
        print("  # of closures = {0}".format(len(self.sidesets)))
        print("  # of coils    = {0}".format(len(self.coil_sets) if self.coil_sets is not None else 0))
        if self.props:
            print("  ThinCurr props: {0}".format(", ".join(sorted(self.props))))
        if self.mesh_props:
            print("  Mesh props: {0}".format(", ".join(sorted(self.mesh_props))))

    def _drop_props(self, operation):
        '''!Clear ThinCurr properties that a structural operation would invalidate'''
        if self.props:
            print("  Warning: dropping ThinCurr properties ({0}) that cannot be "
                  "remapped through {1}".format(", ".join(sorted(self.props)), operation))
            self.props = {}

    # ------------------------------------------------------------- editing
    def transform(self, shift=None, rotate=None, scale=None, center=None):
        '''!Apply a transform to the vertex coordinates (in place)

        Operations, if supplied, are applied in the order rotate, scale, shift:
        ``r' = S R (r - center) + center + shift`` (scaling is about the origin).

        @param shift Translation vector [3] (optional)
        @param rotate Tuple ``(axis, angle_deg)`` where axis is 'x', 'y', or 'z' (optional)
        @param scale Per-axis scale factors [3] applied about the origin (optional)
        @param center Center of rotation [3] (default: origin)
        '''
        if rotate is not None:
            axis, angle_deg = rotate
            rmat = rotation_matrix(axis, angle_deg)
            c = np.zeros(3) if center is None else np.asarray(center, dtype=np.float64)
            self.r = (self.r - c) @ rmat.T + c
        if scale is not None:
            self.r = self.r * np.asarray(scale, dtype=np.float64)
        if shift is not None:
            self.r = self.r + np.asarray(shift, dtype=np.float64)
        return self

    def set_resistivity(self, eta_surf=None, eta_vol=None, thickness=None):
        '''!Set per-region resistivity (and optional thickness) on the mesh (in place)

        Resistivity may be stored as surface resistivity (`eta_surf`, with
        `thickness` optional) or volumetric resistivity (`eta_vol`, which requires
        `thickness`); the two forms are mutually exclusive. Each argument is a
        per-region list; a single value is broadcast to every region. Values are
        stored in the ThinCurr property table and written to the output mesh.

        @param eta_surf Surface resistivity per region (scalar or [nregions])
        @param eta_vol Volumetric resistivity per region (scalar or [nregions])
        @param thickness Thickness per region (scalar or [nregions])
        '''
        if eta_surf is None and eta_vol is None:
            raise ValueError("Resistivity requires either eta_surf or eta_vol")
        if eta_surf is not None and eta_vol is not None:
            raise ValueError("Specify only one of eta_surf or eta_vol, not both")
        if eta_vol is not None and thickness is None:
            raise ValueError("eta_vol requires thickness to also be specified")
        nreg = self.nregions

        def _resolve(name, vals):
            if vals is None:
                return None
            arr = np.atleast_1d(np.asarray(vals, dtype=np.float64)).reshape(-1)
            if arr.shape[0] == 1:  # Broadcast a single value to all regions
                arr = np.full((nreg,), arr[0], dtype=np.float64)
            if arr.shape[0] != nreg:
                raise ValueError("{0} has {1} value(s) but mesh has {2} region(s)".format(
                    name, arr.shape[0], nreg))
            if np.any(arr <= 0.0):
                raise ValueError("All {0} values must be > 0".format(name))
            return arr

        resolved = {'eta_surf': _resolve('eta_surf', eta_surf),
                    'eta_vol': _resolve('eta_vol', eta_vol),
                    'thickness': _resolve('thickness', thickness)}
        # Replace any previously carried resistivity with the requested values
        for key in ('eta_surf', 'eta_vol', 'thickness'):
            self.props.pop(key, None)
            if resolved[key] is not None:
                self.props[key] = resolved[key]
        return self

    def set_jumpers(self, indices):
        '''!Mark the given nodesets as jumpers (all others become holes; in place)

        Indices are 0-based positions into the ordered `nodesets` list and may be
        negative (Python-style backward indexing).

        @param indices Iterable of nodeset indices to mark as jumpers
        '''
        n = len(self.nodesets)
        flags = [False] * n
        for raw in indices:
            i = int(raw)
            if i < -n or i >= n:
                raise ValueError("jumper index {0} is out of range for {1} nodeset(s)".format(raw, n))
            flags[i % n if n > 0 else 0] = True
        self.is_jumper = flags
        return self

    def remove_regions(self, regions):
        '''!Remove all cells belonging to the specified regions (in place)

        Surviving regions are renumbered contiguously (1..N), unreferenced
        vertices are dropped, and nodesets/sidesets are remapped, dropping any
        entries that reference removed nodes/cells.

        @param regions Iterable of region indices to remove
        '''
        remove = set(int(x) for x in regions)
        missing = remove - set(int(x) for x in np.unique(self.reg))
        if missing:
            print("  Warning: region(s) {0} not present in mesh".format(sorted(missing)))
        keep_cell = np.array([r not in remove for r in self.reg], dtype=bool)
        n_removed = int(np.sum(~keep_cell))
        print("  Removing {0} cell(s) in region(s) {1}".format(n_removed, sorted(remove)))
        self._drop_props("region removal")
        # Renumber surviving regions to a contiguous 1..N range
        surviving = np.unique(self.reg[keep_cell]) if np.any(keep_cell) else np.array([], dtype=np.int32)
        remap = {int(old): new + 1 for new, old in enumerate(surviving)}
        self.lc = self.lc[keep_cell, :]
        self.reg = np.array([remap[int(r)] for r in self.reg[keep_cell]], dtype=np.int32)
        # Remap sidesets (cell indices)
        cell_new = -np.ones((keep_cell.shape[0],), dtype=np.int64)
        cell_new[keep_cell] = np.arange(int(np.sum(keep_cell)))
        self.sidesets = _remap_index_sets(self.sidesets, cell_new, "sideset")
        self._reindex_vertices()
        return self

    def _reindex_vertices(self):
        '''!Drop unreferenced vertices and remap connectivity/nodesets (in place)'''
        used = np.zeros((self.np,), dtype=bool)
        if self.nc > 0:
            used[self.lc.reshape(-1)] = True
        new_of_old = -np.ones((self.np,), dtype=np.int64)
        new_of_old[used] = np.arange(int(np.sum(used)))
        self.r = self.r[used, :]
        if self.nc > 0:
            self.lc = new_of_old[self.lc].astype(np.int32)
        # Remap nodesets, keeping the parallel is_jumper flags in sync
        self.nodesets, self.is_jumper = _remap_index_sets(
            self.nodesets, new_of_old, "nodeset", parallel=self.is_jumper)

    def append(self, other, distinct_regions=True):
        '''!Append another mesh onto this one (in place)

        @param other The @ref ThinCurrMesh to append
        @param distinct_regions If True, offset the appended mesh's region indices
            so its regions remain distinct from this mesh's regions
        @result self
        '''
        np_offset = self.np
        nc_offset = self.nc
        reg_offset = self.nregions if distinct_regions else 0
        self.r = np.vstack((self.r, other.r)) if self.np > 0 else other.r.copy()
        self.lc = np.vstack((self.lc, other.lc + np_offset)) if nc_offset > 0 else (other.lc + np_offset)
        self.reg = np.concatenate((self.reg, other.reg + reg_offset))
        for ns, j in zip(other.nodesets, other.is_jumper):
            self.nodesets.append(ns + np_offset)
            self.is_jumper.append(j)
        for ss in other.sidesets:
            self.sidesets.append(ss + nc_offset)
        self._had_native_nodesets = self._had_native_nodesets or other._had_native_nodesets
        if self.props or other.props:
            self._drop_props("mesh combination")
        return self


def rotation_matrix(axis, angle_deg):
    '''!Build a rotation matrix about a principal Cartesian axis

    @param axis Rotation axis: 'x', 'y', or 'z' (case-insensitive)
    @param angle_deg Rotation angle in degrees (right-handed about `axis`)
    @result 3x3 rotation matrix
    '''
    a = str(axis).lower()
    t = np.deg2rad(float(angle_deg))
    c, s = np.cos(t), np.sin(t)
    if a == 'x':
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=np.float64)
    if a == 'y':
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=np.float64)
    if a == 'z':
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=np.float64)
    raise ValueError("Rotation axis must be one of 'x', 'y', or 'z' (got '{0}')".format(axis))


def read_resistivity_xml(filename):
    '''!Read per-region resistivity values from an ``<oft><thincurr>`` XML file

    Values are read from child elements of ``<thincurr>``: surface resistivity
    from ``<eta_surf>`` (or ``<eta>`` as an alias), volumetric resistivity from
    ``<eta_vol>``, and thickness from ``<thickness>``. Each element holds a
    whitespace- or comma-separated list of values (one per region).

    @param filename Path to the XML input file
    @result Tuple ``(eta_surf, eta_vol, thickness)``; each is a list of floats or None
    '''
    xml_doc = ThinCurr_XML.load(filename)

    return xml_doc.eta, xml_doc.eta_vol, xml_doc.thickness


def read_coils_xml(filename):
    '''!Read per-region resistivity values from an ``<oft><thincurr>`` XML file

    Values are read from child elements of ``<thincurr>``: surface resistivity
    from ``<eta_surf>`` (or ``<eta>`` as an alias), volumetric resistivity from
    ``<eta_vol>``, and thickness from ``<thickness>``. Each element holds a
    whitespace- or comma-separated list of values (one per region).

    @param filename Path to the XML input file
    @result Tuple ``(eta_surf, eta_vol, thickness)``; each is a list of floats or None
    '''
    xml_doc = ThinCurr_XML.load(filename)

    return xml_doc.icoils+xml_doc.vcoils


def _remap_index_sets(sets, old_to_new, label, parallel=None):
    '''!Remap and filter a list of index arrays using an old->new lookup

    Entries mapped to a negative value (removed) are dropped. Empty sets that
    result are discarded with a warning. An optional `parallel` list (one entry
    per set) is filtered identically to keep it aligned with the surviving sets.

    @param sets List of index arrays (0-based)
    @param old_to_new Lookup array mapping old index -> new index (<0 = removed)
    @param label Name used in warning messages
    @param parallel Optional list of per-set metadata to filter alongside `sets`
    @result Filtered/remapped list of index arrays; or ``(sets, parallel)`` if
        `parallel` was supplied
    '''
    out = []
    out_parallel = []
    for i, s in enumerate(sets):
        mapped = old_to_new[s]
        keep = mapped >= 0
        if not np.all(keep):
            print("  Warning: dropping {0} entrie(s) from {1} {2}".format(
                int(np.sum(~keep)), label, i + 1))
        mapped = mapped[keep].astype(np.int32)
        if mapped.shape[0] > 0:
            out.append(mapped)
            if parallel is not None:
                out_parallel.append(parallel[i])
        else:
            print("  Warning: {0} {1} became empty and was dropped".format(label, i + 1))
    if parallel is not None:
        return out, out_parallel
    return out


def resolve_jumper_indices(nnodesets, explicit=None, index_range=None):
    '''!Resolve jumper nodeset indices from CLI options

    @param nnodesets Number of nodesets available to categorize
    @param explicit Iterable of explicit 0-based indices (negative allowed)
    @param index_range ``(start, stop)`` half-open range with Python slice
        semantics (negative allowed)
    @result Sorted list of resolved 0-based jumper indices
    '''
    if nnodesets == 0:
        raise ValueError("no nodesets are available to mark as jumpers")
    if explicit is not None:
        out = set()
        for raw in explicit:
            i = int(raw)
            if i < -nnodesets or i >= nnodesets:
                raise ValueError("jumper index {0} is out of range for {1} "
                                 "nodeset(s)".format(raw, nnodesets))
            out.add(i % nnodesets)
        return sorted(out)
    if index_range is not None:
        start = int(index_range[0])
        stop = int(index_range[1]) if len(index_range) > 1 else None  # None -> to end
        return list(range(nnodesets))[start:stop]  # Python slice semantics
    return []


def combine_meshes(filenames, distinct_regions=True):
    '''!Combine multiple mesh files into a single mesh

    Aborts if any input file contains native ``NODESETXXXX`` entries, since such
    uncategorized loops cannot be meaningfully combined; run the `modify` workflow
    first to convert them into ThinCurr holes/jumpers.

    @param filenames List of input file paths
    @param distinct_regions Keep each input's regions distinct by offsetting IDs
    @result Combined @ref ThinCurrMesh
    '''
    meshes = []
    for fn in filenames:
        mesh = ThinCurrMesh.load(fn)
        if mesh._had_native_nodesets:
            raise ValueError(
                "'{0}' contains native NODESET entries; combine does not support "
                "uncategorized nodesets. Run 'modify' on it first to convert them "
                "into ThinCurr holes/jumpers.".format(fn))
        meshes.append(mesh)
    combined = meshes[0].copy()
    for mesh in meshes[1:]:
        combined.append(mesh, distinct_regions=distinct_regions)
    return combined


def replicate_mesh(mesh, ncopies, shift=None, rotate=None, center=None,
                   distinct_regions=True):
    '''!Generate an output mesh consisting of the original plus transformed copies

    The untransformed original is always kept. Copy ``k`` (k = 1 .. ncopies) is
    the input mesh with the incremental transform applied ``k`` times: shifted by
    ``k*shift`` and/or rotated by ``k*angle`` about `center`. The output therefore
    contains ``ncopies + 1`` instances.

    @param mesh Input @ref ThinCurrMesh
    @param ncopies Number of transformed copies to add (excludes the original)
    @param shift Per-copy incremental translation [3] (optional)
    @param rotate Per-copy incremental rotation ``(axis, angle_deg)`` (optional)
    @param center Center of rotation [3] (default: origin)
    @param distinct_regions If True, offset region IDs so every copy is distinct
    @result Replicated @ref ThinCurrMesh
    '''
    if ncopies < 1:
        raise ValueError("Number of copies must be >= 1")
    if shift is None and rotate is None:
        raise ValueError("Replication requires a --shift or --rotate increment")
    out = None
    for i in range(ncopies + 1):  # i=0 is the untransformed original
        piece = mesh.copy()
        step_shift = None if shift is None else np.asarray(shift, dtype=np.float64) * i
        step_rotate = None if rotate is None else (rotate[0], float(rotate[1]) * i)
        piece.transform(shift=step_shift, rotate=step_rotate, center=center)
        if out is None:
            out = piece
        else:
            out.append(piece, distinct_regions=distinct_regions)
    return out


def _add_common_out(parser):
    parser.add_argument("--out_file", type=str, default=None, help="Output mesh file")


def build_parser():
    '''!Construct the command-line argument parser'''
    import argparse
    parser = argparse.ArgumentParser(
        description="Manipulate and combine native-format 3-node triangular ThinCurr meshes")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- combine ---
    p_comb = sub.add_parser("combine", help="Combine two or more mesh files")
    p_comb.add_argument("--in_files", type=str, nargs='+', required=True,
                        help="Input mesh files (two or more)")
    p_comb.add_argument("--merge_regions", action="store_true", default=False,
                        help="Merge identical region IDs across inputs instead of "
                             "keeping them distinct (default: keep distinct)")
    _add_common_out(p_comb)

    # --- modify ---
    p_mod = sub.add_parser("modify", help="Remove regions, transform, and/or replicate a mesh")
    p_mod.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    p_mod.add_argument("--remove_regions", type=int, nargs='+', default=None,
                       help="Region indices to remove")
    xform = p_mod.add_mutually_exclusive_group()
    xform.add_argument("--shift", type=float, nargs=3, default=None,
                       metavar=('X', 'Y', 'Z'),
                       help="Translation [X Y Z]. Applied once to the whole mesh, or as "
                            "the per-copy increment when --copies is given")
    xform.add_argument("--rotate", nargs=2, default=None, metavar=('AXIS', 'ANGLE'),
                       help="Rotation of ANGLE degrees about AXIS (x|y|z). Applied once to "
                            "the whole mesh, or as the per-copy increment when --copies is given")
    xform.add_argument("--scale", type=float, nargs=3, default=None,
                       metavar=('SX', 'SY', 'SZ'),
                       help="Scale factors along X, Y, Z about the origin. Cannot be "
                            "combined with --copies")
    p_mod.add_argument("--rotate_center", type=float, nargs=3, default=None,
                       metavar=('X', 'Y', 'Z'), help="Center of rotation (default: origin)")
    p_mod.add_argument("--copies", type=int, default=None,
                       help="Add this many transformed copies, not counting the original "
                            "(which is always kept), applying --shift or --rotate "
                            "incrementally to each (requires --shift or --rotate; not "
                            "valid with --scale)")
    p_mod.add_argument("--persist_regions", action="store_true", default=False,
                       help="Keep the original region IDs on every copy instead of "
                            "offsetting them so each copy is distinct (default: offset "
                            "so each copy gets its own distinct regions)")
    p_mod.add_argument("--eta_surf", type=float, nargs='+', default=None,
                       help="Surface resistivity per region (one value, or one per "
                            "region); thickness is optional")
    p_mod.add_argument("--eta_vol", type=float, nargs='+', default=None,
                       help="Volumetric resistivity per region (one value, or one per "
                            "region); requires --thickness")
    p_mod.add_argument("--thickness", type=float, nargs='+', default=None,
                       help="Region thickness (one value, or one per region)")
    p_mod.add_argument("--eta_from_xml", type=str, default=None,
                       help="Read eta/eta_surf, eta_vol, and thickness from an "
                            "<oft><thincurr> XML file instead of the flags above")
    p_mod.add_argument("--coils_from_xml", type=str, default=None,
                        help="Read coil sets from an "
                            "<oft><thincurr> XML file and add to the output mesh")
    jumps = p_mod.add_mutually_exclusive_group()
    jumps.add_argument("--jumpers", type=int, nargs='+', default=None, metavar='IDX',
                       help="Nodeset indices to treat as jumpers instead of holes "
                            "(0-based; negative indices count from the end). All other "
                            "nodesets are holes")
    jumps.add_argument("--jumper_range", type=int, nargs='+', default=None,
                       metavar='N',
                       help="Range of nodeset indices to treat as jumpers, given as "
                            "START [STOP] with Python slice semantics [START:STOP) "
                            "(0-based; negative allowed; omit STOP to run to the end). "
                            "Use '0' to mark every nodeset as a jumper")
    _add_common_out(p_mod)
    return parser


def _parse_rotate(val):
    '''!Validate and normalize a (axis, angle) argument pair'''
    if val is None:
        return None
    axis, angle = val
    rotation_matrix(axis, angle)  # validates axis and angle
    return (str(axis).lower(), float(angle))


def script_entry(argv=None):
    '''!Command-line entry point'''
    parser = build_parser()
    options = parser.parse_args(argv)

    if options.command == "combine":
        if len(options.in_files) < 2:
            parser.error("combine requires at least two --in_files")
        out_file = options.out_file
        if out_file is None:
            out_file = os.path.splitext(options.in_files[0])[0] + "-combined.h5"
        try:
            result = combine_meshes(options.in_files,
                                    distinct_regions=not options.merge_regions)
        except ValueError as err:
            parser.error(str(err))
        result.save(out_file)

    elif options.command == "modify":
        try:
            rotate = _parse_rotate(options.rotate)
        except ValueError as err:
            parser.error(str(err))
        # --shift/--rotate/--scale are mutually exclusive via argparse; only the
        # additional constraints tied to --copies need checking here.
        if options.copies is not None:
            if options.scale is not None:
                parser.error("--scale cannot be combined with --copies")
            if options.shift is None and rotate is None:
                parser.error("--copies requires --shift or --rotate")
        if options.jumper_range is not None and not 1 <= len(options.jumper_range) <= 2:
            parser.error("--jumper_range takes START or START STOP (1 or 2 values)")
        # Resolve resistivity specification (direct flags or from an XML file)
        direct_eta = (options.eta_surf, options.eta_vol, options.thickness)
        res_spec = None
        if options.eta_from_xml is not None:
            if any(v is not None for v in direct_eta):
                parser.error("--eta_from_xml cannot be combined with "
                             "--eta_surf/--eta_vol/--thickness")
            try:
                res_spec = read_resistivity_xml(options.eta_from_xml)
            except (ValueError, OSError, ET.ParseError) as err:
                parser.error("Failed to read resistivity from XML: {0}".format(err))
            if res_spec[0] is None and res_spec[1] is None:
                parser.error("No <eta>/<eta_surf> or <eta_vol> found in "
                             "'{0}'".format(options.eta_from_xml))
        elif any(v is not None for v in direct_eta):
            res_spec = direct_eta
        #
        coil_sets = None
        if options.coils_from_xml is not None:
            try:
                coil_sets = read_coils_xml(options.coils_from_xml)
            except (ValueError, OSError, ET.ParseError) as err:
                parser.error("Failed to read coils from XML: {0}".format(err))
        out_file = options.out_file
        if out_file is None:
            out_file = os.path.splitext(options.in_file)[0] + "-modified.h5"

        mesh = ThinCurrMesh.load(options.in_file)
        if options.remove_regions is not None:
            print()
            print("Removing regions...")
            mesh.remove_regions(options.remove_regions)
        if options.copies is not None:
            print()
            print("Replicating mesh...")
            mesh = replicate_mesh(mesh, options.copies,
                                  shift=options.shift, rotate=rotate,
                                  center=options.rotate_center,
                                  distinct_regions=not options.persist_regions)
        else:
            if options.shift is not None or rotate is not None or options.scale is not None:
                print()
                print("Applying transform...")
                mesh.transform(shift=options.shift, rotate=rotate,
                               scale=options.scale, center=options.rotate_center)
        # Categorize jumpers against the final nodeset list (default: all holes)
        if options.jumpers is not None or options.jumper_range is not None:
            print()
            print("Categorizing jumpers...")
            try:
                jidx = resolve_jumper_indices(len(mesh.nodesets),
                                              explicit=options.jumpers,
                                              index_range=options.jumper_range)
                mesh.set_jumpers(jidx)
            except ValueError as err:
                parser.error(str(err))
            print("  {0} hole(s), {1} jumper(s)".format(mesh.nholes, mesh.njumpers))
        # Set resistivity last, matched against the final region count
        if res_spec is not None:
            print()
            print("Setting resistivity...")
            try:
                mesh.set_resistivity(eta_surf=res_spec[0], eta_vol=res_spec[1],
                                     thickness=res_spec[2])
            except ValueError as err:
                parser.error(str(err))
        # Set coils
        if coil_sets is not None:
            print()
            print("Adding coil sets...")
            mesh.coil_sets = mesh.coil_sets + coil_sets if mesh.coil_sets is not None else coil_sets
        mesh.save(out_file)

#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! ThinCurr mesh manipulation and combination utility script and supporting functions

See @ref thincurr_mesh_tool

@authors Open FUSION Toolkit contributors
@date July 2026
'''
import os
import xml.etree.ElementTree as ET
import numpy as np
from .meshing import read_ThinCurr_mesh, write_ThinCurr_mesh
from .coils import ThinCurr_XML


class ThinCurrMesh:
    '''! Container for a 3-node triangular ThinCurr surface mesh'''

    #: ThinCurr properties preserved verbatim on a geometry-only round-trip
    _PROP_KEYS = ('eta_surf', 'eta_vol', 'thickness') # Ignore 'pmap' and 'nfp' for now (only limited uses)
    _MESH_PROP_KEYS = ('reg_names', 'reg_attrs')

    def __init__(self, r, lc, reg, holes=None, jumpers=None, closures=None,
                 nodesets=None, tc_props=None, mesh_props=None,
                 coil_sets=None):
        '''! Construct a mesh from its component arrays

        @param r Vertex list [np,3]
        @param lc Triangle vertex list [nc,3], 0-based
        @param reg Per-cell region index [nc]
        @param holes List of 0-based node-index arrays for holes (optional)
        @param jumpers List of 0-based node-index arrays for jumpers (optional)
        @param closures List of 0-based cell-index arrays for closures (optional)
        @param nodesets List of 0-based node-index arrays for unclassified nodesets (optional)
        @param tc_props Dict of ThinCurr properties to carry through (optional)
        @param mesh_props Dict of base mesh properties to carry through (optional)
        @param coil_sets List of ThinCurr coil sets (optional)
        '''
        self.r = np.ascontiguousarray(r, dtype=np.float64)
        if self.r.ndim != 2 or self.r.shape[1] not in (2, 3):
            raise ValueError("Vertex list must have shape [np,2] or [np,3]")
        if self.r.shape[1] == 2:  # Upgrade 2D point list to 3D (Z=0)
            self.r = np.hstack((self.r, np.zeros((self.r.shape[0], 1))))
        #
        self.lc = np.ascontiguousarray(lc, dtype=np.int32)
        if self.lc.ndim != 2 or self.lc.shape[1] != 3:
            raise ValueError("Cell list must have shape [nc,3] (3-node triangles)")
        #
        self.reg = np.ascontiguousarray(reg, dtype=np.int32).reshape(-1)
        if self.reg.shape[0] != self.lc.shape[0]:
            raise ValueError("Region list length must match number of cells")
        #
        if nodesets is not None and (holes is not None or jumpers is not None):
            raise ValueError("Cannot supply both nodesets and holes/jumpers")
        #
        self.holes = [np.ascontiguousarray(v, dtype=np.int32).reshape(-1) for v in holes] if holes is not None else None
        self.jumpers = [np.ascontiguousarray(v, dtype=np.int32).reshape(-1) for v in jumpers] if jumpers is not None else None
        self.closures = np.ascontiguousarray(closures, dtype=np.int32) if closures is not None else None
        self.nodesets = [np.ascontiguousarray(v, dtype=np.int32).reshape(-1) for v in nodesets] if nodesets is not None else None
        #
        self.tc_props = dict(tc_props) if tc_props else {}
        self.mesh_props = dict(mesh_props) if mesh_props else {}
        self.coil_sets = coil_sets

    @property
    def np(self):
        '''! Number of vertices'''
        return self.r.shape[0]

    @property
    def nc(self):
        '''! Number of cells (triangles)'''
        return self.lc.shape[0]

    @property
    def nregions(self):
        '''! Number of regions (maximum region index)'''
        return self.reg.max()

    @property
    def nholes(self):
        '''! Number of hole nodesets'''
        if self.nodesets is not None:
            return len(self.nodesets)
        else:
            return len(self.holes) if self.holes is not None else 0

    @property
    def nclosures(self):
        '''! Number of closure nodesets'''
        return self.closures.shape[0] if self.closures is not None else 0

    @property
    def njumpers(self):
        '''! Number of jumper nodesets'''
        return len(self.jumpers) if self.jumpers is not None else 0

    def copy(self):
        '''! Return a deep copy of the mesh'''
        new = ThinCurrMesh(
            self.r.copy(),
            self.lc.copy(),
            self.reg.copy(),
            holes=[h.copy() for h in self.holes] if self.holes is not None else None,
            jumpers=[j.copy() for j in self.jumpers] if self.jumpers is not None else None,
            closures=self.closures.copy() if self.closures is not None else None,
            nodesets=[ns.copy() for ns in self.nodesets],
            tc_props={k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.tc_props.items()},
            mesh_props={k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.mesh_props.items()},
            coil_sets=self.coil_sets)
        return new

    # ------------------------------------------------------------------ I/O
    @classmethod
    def load(cls, filename):
        '''! Load a mesh using @ref OpenFUSIONToolkit.ThinCurr.meshing.read_ThinCurr_mesh "read_ThinCurr_mesh" and
        create instance of @ref ThinCurrMesh.

        @param filename Path to input mesh file
        @result New @ref ThinCurrMesh instance
        '''
        print()
        print("Reading mesh: {0}".format(filename))
        mesh_info = read_ThinCurr_mesh(filename)
        if 'ho_info' in mesh_info:
            raise ValueError("High-order meshes are not supported (found high-order node info)")
        # mesh_type = mesh_info.get('type', 'tri')
        r = mesh_info.pop('r')
        lc = mesh_info.pop('lc')
        reg = mesh_info.pop('reg')
        if r.ndim == 2 and r.shape[1] == 2:
            print("  Note: input mesh is 2D; padding point list to 3D (Z=0)")
        tc = mesh_info.get('thincurr', {})

        # Read holes and jumpers, falling back to nodesets for old meshes
        nodesets = None
        holes = None
        jumpers = None
        if 'nodesets' in mesh_info:
            if ('holes' in tc) or ('jumpers' in tc):
                raise ValueError("Input mesh contains both native NODESETs and ThinCurr holes/jumpers")
            nodesets = mesh_info.get('nodesets')
        else:
            holes = tc.get('holes')
            jumpers = tc.get('jumpers')

        # Read closures, falling back to sideset 1 for old meshes
        closures = None
        if 'sidesets' in mesh_info:
            if 'closures' in tc:
                raise ValueError("Input mesh contains both native side sets and ThinCurr closures")
            closures = mesh_info['sidesets'][0]
        else:
            closures = tc.get('closures')

        # Copy other fields to keep
        tc_props = {k: tc[k] for k in cls._PROP_KEYS if k in tc}
        mesh_props = {k: mesh_info[k] for k in cls._MESH_PROP_KEYS if k in mesh_info}

        # Create mesh object
        mesh = cls(r, lc, reg,
                   holes=holes, jumpers=jumpers, closures=closures,
                   nodesets=nodesets,
                   tc_props=tc_props,
                   mesh_props=mesh_props,
                   coil_sets=tc.get('coil_sets'))

        # Display info and return
        mesh.print_info()
        return mesh

    def save(self, filename):
        '''! Save the mesh contained in this object using
        @ref OpenFUSIONToolkit.ThinCurr.meshing.write_ThinCurr_mesh "write_ThinCurr_mesh"

        @param filename Path to output mesh file
        '''
        # Convert from legacy to modern ThinCurr representation if needed
        if self.holes is None:
            self.holes = self.nodesets
            self.nodesets = None
        write_ThinCurr_mesh(filename, self.r, self.lc, self.reg,
                            reg_attrs=self.mesh_props.get('reg_attrs'),
                            reg_names=self.mesh_props.get('reg_names'),
                            holes=self.holes, jumpers=self.jumpers, closures=self.closures,
                            eta_surf=self.tc_props.get('eta_surf'),
                            eta_vol=self.tc_props.get('eta_vol'),
                            thickness=self.tc_props.get('thickness'),
                            coil_sets=self.coil_sets,
                            pmap=self.tc_props.get('pmap'),
                            nfp=self.tc_props.get('nfp'))
        self.print_info()

    def print_info(self):
        '''! Print a short summary of the mesh contents'''
        print("  # of points   = {0}".format(self.np))
        print("  # of cells    = {0}".format(self.nc))
        print("  # of regions  = {0}".format(self.nregions))
        print("  # of holes    = {0}".format(self.nholes))
        print("  # of jumpers  = {0}".format(self.njumpers))
        print("  # of closures = {0}".format(self.nclosures))
        print("  # of coils    = {0}".format(len(self.coil_sets) if self.coil_sets is not None else 0))
        if self.tc_props:
            print("  ThinCurr props: {0}".format(", ".join(sorted(self.tc_props))))
        if self.mesh_props:
            print("  Mesh props: {0}".format(", ".join(sorted(self.mesh_props))))

    def _drop_props(self, operation):
        '''! Clear ThinCurr properties that a structural operation would invalidate'''
        if self.tc_props:
            print("  Warning: dropping ThinCurr properties ({0}) that cannot be "
                  "remapped through {1}".format(", ".join(sorted(self.tc_props)), operation))
            self.tc_props = {}

    # ------------------------------------------------------------- editing
    def transform(self, shift=None, rotate=None, scale=None, center=None):
        '''! Apply a transform to the mesh coordinates (in place)

        Multiple operations, if supplied, are applied in the order `rotate`, `scale`, `shift`
        (scaling is always about the origin)

        @param shift Translation vector [3] (optional)
        @param rotate Tuple `(axis, angle_deg)` where axis is `'x'`, `'y'`, or `'z'` (optional)
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
        '''! Set per-region resistivity (and optional thickness) on the mesh (in place)

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
            self.tc_props.pop(key, None)
            if resolved[key] is not None:
                self.tc_props[key] = resolved[key]

    def set_jumpers(self, indices):
        '''! Mark the given nodesets as jumpers (all others become holes; in place)

        Indices are 0-based positions into the ordered `nodesets` list and may be
        negative (Python-style backward indexing).

        @param indices Iterable of nodeset indices to mark as jumpers
        '''
        if self.nodesets is None:
            raise ValueError("Cannot set jumpers on a mesh with no nodesets")
        n = len(self.nodesets)
        flags = [False] * n
        for raw in indices:
            i = int(raw)
            if i < -n or i >= n:
                raise ValueError("jumper index {0} is out of range for {1} nodeset(s)".format(raw, n))
            flags[i % n if n > 0 else 0] = True
        # Perform conversion
        self.holes = [ns for i, ns in enumerate(self.nodesets) if not flags[i]]
        self.jumpers = [ns for i, ns in enumerate(self.nodesets) if flags[i]]
        self.nodesets = None

    def remove_regions(self, regions):
        '''! Remove all cells belonging to the specified regions (in place)

        Surviving regions are renumbered contiguously (1..N), unreferenced
        vertices are dropped, and nodesets/sidesets are remapped, dropping any
        entries that reference removed nodes/cells.

        @param regions Iterable of region indices to remove
        '''
        def _remap_index_sets(sets, old_to_new, label):
            out = []
            for i, s in enumerate(sets):
                mapped = old_to_new[s]
                keep = mapped >= 0
                if not np.all(keep):
                    print("  Warning: dropping {0} entrie(s) from {1} {2}".format(
                        int(np.sum(~keep)), label, i + 1))
                mapped = mapped[keep].astype(np.int32)
                if mapped.shape[0] > 0:
                    out.append(mapped)
                else:
                    print("  Warning: {0} {1} became empty and was dropped".format(label, i + 1))
            return out
        #
        if self.nodesets is not None:
            raise ValueError(
                "Mesh contains native NODESET entries; region removal does not support "
                "uncategorized nodesets. Run 'modify' on it first to convert them "
                "into ThinCurr holes/jumpers.")
        remove = set(int(x) for x in regions)
        missing = remove - set(int(x) for x in np.unique(self.reg))
        if missing:
            print("  Warning: region(s) {0} not present in mesh".format(sorted(missing)))
        keep_cell = np.array([r not in remove for r in self.reg], dtype=bool)
        n_removed = int(np.sum(~keep_cell))
        print("  Removing {0} cell(s) in region(s) {1}".format(n_removed, sorted(remove)))

        # Renumber surviving regions to a contiguous 1..N range
        surviving = np.unique(self.reg[keep_cell]) if np.any(keep_cell) else np.array([], dtype=np.int32)
        remap = {int(old): new + 1 for new, old in enumerate(surviving)}
        self.lc = self.lc[keep_cell, :]
        self.reg = np.array([remap[int(r)] for r in self.reg[keep_cell]], dtype=np.int32)
        # Remap eta_surf, eta_vol, thickness properties (if any)
        for key in ('eta_surf', 'eta_vol', 'thickness'):
            if key in self.tc_props:
                arr = self.tc_props[key]
                self.tc_props[key] = np.array([arr[old - 1] for old in surviving], dtype=np.float64)
        # Remap sidesets (cell indices)
        cell_new = -np.ones((keep_cell.shape[0],), dtype=np.int64)
        cell_new[keep_cell] = np.arange(int(np.sum(keep_cell)))
        self.closures = _remap_index_sets(self.closures, cell_new, "closures")

        # Remove unreferenced vertices and remap connectivity/nodesets (in place)
        used = np.zeros((self.np,), dtype=bool)
        if self.nc > 0:
            used[self.lc.reshape(-1)] = True
        new_of_old = -np.ones((self.np,), dtype=np.int64)
        new_of_old[used] = np.arange(int(np.sum(used)))
        self.r = self.r[used, :]
        if self.nc > 0:
            self.lc = new_of_old[self.lc].astype(np.int32)
        # Remap holes/jumpers
        self.holes = _remap_index_sets(self.holes, new_of_old, "holes")
        self.jumpers = _remap_index_sets(self.jumpers, new_of_old, "jumpers")

    def append(self, other, distinct_regions=True):
        '''! Append another mesh onto this one (in place)

        @param other The @ref ThinCurrMesh to append
        @param distinct_regions If True, offset the appended mesh's region indices
            so its regions remain distinct from this mesh's regions
        '''
        np_offset = self.np
        nc_offset = self.nc
        reg_offset = self.nregions if distinct_regions else 0
        self.r = np.vstack((self.r, other.r)) if self.np > 0 else other.r.copy()
        self.lc = np.vstack((self.lc, other.lc + np_offset)) if nc_offset > 0 else (other.lc + np_offset)
        self.reg = np.concatenate((self.reg, other.reg + reg_offset))
        # Append holes
        if self.holes is None:
            self.holes = [ns.copy() for ns in other.holes] if other.holes is not None else None
        else:
            self.holes = [ns+np_offset for ns in other.holes] if other.holes is not None else self.holes
        # Append jumpers
        if self.jumpers is None:
            self.jumpers = [ns.copy() for ns in other.jumpers] if other.jumpers is not None else None
        else:
            self.jumpers = [ns+np_offset for ns in other.jumpers] if other.jumpers is not None else self.jumpers
        # Append closures
        if self.closures is None:
            self.closures = other.closures.copy() if other.closures is not None else None
        else:
            self.closures = np.concatenate((self.closures, other.closures + nc_offset)) if other.closures is not None else self.closures
        # Append resistivity properties (if any) and warn if they exist in either mesh
        if distinct_regions:
            if 'eta_surf' in self.tc_props and 'eta_surf' in other.tc_props:
                self.tc_props['eta_surf'] = np.concatenate((self.tc_props['eta_surf'], other.tc_props['eta_surf']))
            if 'eta_vol' in self.tc_props and 'eta_vol' in other.tc_props:
                self.tc_props['eta_vol'] = np.concatenate((self.tc_props['eta_vol'], other.tc_props['eta_vol']))
            if 'thickness' in self.tc_props and 'thickness' in other.tc_props:
                self.tc_props['thickness'] = np.concatenate((self.tc_props['thickness'], other.tc_props['thickness']))
        else:
            self._drop_props("append with non-distinct regions")


def rotation_matrix(axis, angle_deg):
    '''! Build a rotation matrix about a principal Cartesian axis

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


def resolve_jumper_indices(nnodesets, explicit=None, index_range=None):
    '''! Resolve jumper nodeset indices from CLI options

    @param nnodesets Number of nodesets available to categorize
    @param explicit Iterable of explicit 0-based indices (negative allowed)
    @param index_range `(start, stop)` half-open range with Python slice
        semantics (negative allowed)
    @result Sorted list of resolved 0-based jumper indices
    '''
    if nnodesets == 0:
        raise ValueError("no nodesets are available to mark as jumpers")
    if explicit is not None:
        out = set()
        for raw in explicit:
            i = int(raw)
            if (i < -nnodesets) or (i >= nnodesets):
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
    '''! Combine multiple mesh files into a single mesh

    Aborts if any input file contains native `NODESET` entries, since such
    uncategorized loops cannot be meaningfully combined; run the `modify` workflow
    first to convert them into ThinCurr holes/jumpers.

    @param filenames List of input file paths
    @param distinct_regions Keep each input's regions distinct by offsetting IDs
    @result Combined @ref ThinCurrMesh
    '''
    meshes = []
    for fn in filenames:
        mesh = ThinCurrMesh.load(fn)
        if mesh.nodesets is not None:
            raise ValueError(
                "'{0}' contains native NODESET entries; combine does not support "
                "uncategorized nodesets. Run 'modify' on it first to convert them "
                "into ThinCurr holes/jumpers.".format(fn))
        meshes.append(mesh)
    combined = meshes[0].copy()
    for mesh in meshes[1:]:
        combined.append(mesh, distinct_regions=distinct_regions)
        # Warn about uncombined coils
        if combined.coil_sets is not None:
            print("  Warning: Multiple input meshes contain ThinCurr coil sets; "
                  "these are not combined and only the first mesh's coil sets are preserved.")
    return combined


def replicate_mesh(mesh, ncopies, shift=None, rotate=None, center=None,
                   distinct_regions=True):
    '''! Generate an output mesh consisting of the original plus transformed copies

    The untransformed original is always kept. Copy `k` (`k = 1..ncopies`) is
    the input mesh with the incremental transform applied `k` times: shifted by
    `k*shift` and/or rotated by `k*angle` about `center`. The output therefore
    contains `ncopies + 1` instances.

    @param mesh Input @ref ThinCurrMesh
    @param ncopies Number of transformed copies to add (excludes the original)
    @param shift Per-copy incremental translation [3] (optional)
    @param rotate Per-copy incremental rotation `(axis, angle_deg)` (optional)
    @param center Center of rotation [3] (default: origin)
    @param distinct_regions If True, offset region IDs so every copy is distinct
    @result Replicated @ref ThinCurrMesh
    '''
    if mesh.nodesets is not None:
        raise ValueError(
            "Mesh contains native NODESET entries; replication does not support "
            "uncategorized nodesets. Run 'modify' on it first to convert them "
            "into ThinCurr holes/jumpers.")
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
    p_comb.add_argument("--out_file", type=str, default=None, help="Output mesh file")
    p_comb.add_argument("--merge_regions", action="store_true", default=False,
                        help="Merge identical region IDs across inputs instead of "
                             "keeping them distinct (default: keep distinct)")

    # --- modify ---
    p_mod = sub.add_parser("modify", help="Remove regions, transform, and/or replicate a mesh")
    p_mod.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    p_mod.add_argument("--out_file", type=str, default=None, help="Output mesh file")
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
    return parser


def _parse_rotate(val):
    '''!Validate and normalize a (axis, angle) argument pair'''
    if val is None:
        return None
    axis, angle = val
    rotation_matrix(axis, angle)  # validates axis and angle
    return (str(axis).lower(), float(angle))


def script_entry(argv=None):
    '''!Command-line entry point

    See @ref thincurr_mesh_tool
    '''
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
        res_spec = (options.eta_surf, options.eta_vol, options.thickness)
        if options.eta_from_xml is not None:
            if any(v is not None for v in res_spec):
                parser.error("--eta_from_xml cannot be combined with "
                             "--eta_surf/--eta_vol/--thickness")
            try:
                xml_doc = ThinCurr_XML.load(options.eta_from_xml)
            except (ValueError, OSError, ET.ParseError) as err:
                parser.error("Failed to read resistivity from XML: {0}".format(err))
            if (xml_doc.eta is None) and (xml_doc.eta_vol is None):
                parser.error("No <eta>/<eta_surf> or <eta_vol> found in "
                             "'{0}'".format(options.eta_from_xml))
            res_spec = (xml_doc.eta, xml_doc.eta_vol, xml_doc.thickness)
        #
        coil_sets = None
        if options.coils_from_xml is not None:
            try:
                xml_doc = ThinCurr_XML.load(options.coils_from_xml)
                coil_sets = xml_doc.icoils+xml_doc.vcoils
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


## @page thincurr_mesh_tool `OFT_ThinCurr_mesh_tool`: ThinCurr mesh manipulation utility
#
# @tableofcontents
#
# @section thincurr_mesh_tool_desc Description
# This script supports a variety of operations to add/modify information about ThinCurr models
# as well as combining multiple models. It supports two top-level
# workflows via sub-commands:
#
# Modify a single mesh file by removing regions, applying a rigid transform, and/or
# generating multiple shifted/rotated copies:
#
#```shell
# # Convert every nodeset in the file to a jumper
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --jumper_range 0
#
# # Set a uniform surface resistivity on all regions
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --eta_surf 1.257E-5
#
# # Copy resistivity from a ThinCurr XML input file
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --eta_from_xml oft_in.xml
#
# # Set per-region volumetric resistivity with thickness (2 regions)
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --eta_vol 1.257E-5 2.5E-5 --thickness 1.0E-3 2.0E-3
#
# # Remove regions 2 and 3
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --remove_regions 2 3
#
# # Shift the whole mesh by (0,0,0.5)
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --shift 0 0 0.5
#
# # Stretch the mesh by 2x along x (scale about the origin)
# OFT_ThinCurr_mesh_tool modify --in_file input_model.h5 --scale 2 1 1
#
# # Build a 6-fold toroidal array (60 deg spacing about z):
# # the original plus 5 rotated copies at 60,120,180,240,300 deg
# OFT_ThinCurr_mesh_tool modify --in_file segment_model.h5 --copies 5 --rotate z 60
#```
#
# Combine two (or more) existing mesh files into a single mesh:
#
#```shell
#OFT_ThinCurr_mesh_tool combine --in_files model1.h5 model2.h5 --out_file combined_model.h5
#```
#
# In general, the \ref thincurr_compute_holes "OFT_thincurr_holes" utility should be used on the final
# model to compute holes and colusures from the combined geometry. For compatibility with older models,
# this script treats all NODESETs as holes by default, which will be overriden by the `OFT_thincurr_holes`
# if run. For future meshes only NODESETs corresponding to jumpers should be specified, which can
# be converted to jumpers instead of holes using the `--jumpers` or `--jumper_range` options.
#
# Only one of `--shift`, `--rotate`, or `--scale` may be given at a time. When `--copies`
# is given, `--shift` or `--rotate` defines the per-copy increment (`--scale` is not
# allowed); otherwise the transform is applied once to the whole mesh.
#
# @section thincurr_mesh_tool_opts Script Options
# Script options for `modify` workflow:
#
#```shell
#usage: OFT_ThinCurr_mesh_tool modify [-h] --in_file IN_FILE [--out_file OUT_FILE] [--remove_regions REMOVE_REGIONS [REMOVE_REGIONS ...]] [--shift X Y Z |
#                                        --rotate AXIS ANGLE | --scale SX SY SZ] [--rotate_center X Y Z] [--copies COPIES] [--persist_regions]
#                                        [--eta_surf ETA_SURF [ETA_SURF ...]] [--eta_vol ETA_VOL [ETA_VOL ...]] [--thickness THICKNESS [THICKNESS ...]]
#                                        [--eta_from_xml ETA_FROM_XML] [--coils_from_xml COILS_FROM_XML] [--jumpers IDX [IDX ...] | --jumper_range N [N ...]]
#
#options:
#  -h, --help            show this help message and exit
#  --in_file IN_FILE     Input mesh file
#  --out_file OUT_FILE   Output mesh file
#  --remove_regions REMOVE_REGIONS [REMOVE_REGIONS ...]
#                        Region indices to remove
#  --shift X Y Z         Translation [X Y Z]. Applied once to the whole mesh, or as the per-copy increment when --copies is given
#  --rotate AXIS ANGLE   Rotation of ANGLE degrees about AXIS (x|y|z). Applied once to the whole mesh, or as the per-copy increment when --copies is given
#  --scale SX SY SZ      Scale factors along X, Y, Z about the origin. Cannot be combined with --copies
#  --rotate_center X Y Z
#                        Center of rotation (default: origin)
#  --copies COPIES       Add this many transformed copies, not counting the original (which is always kept), applying --shift or --rotate incrementally to each
#                        (requires --shift or --rotate; not valid with --scale)
#  --persist_regions     Keep the original region IDs on every copy instead of offsetting them so each copy is distinct (default: offset so each copy gets its own
#                        distinct regions)
#  --eta_surf ETA_SURF [ETA_SURF ...]
#                        Surface resistivity per region (one value, or one per region); thickness is optional
#  --eta_vol ETA_VOL [ETA_VOL ...]
#                        Volumetric resistivity per region (one value, or one per region); requires --thickness
#  --thickness THICKNESS [THICKNESS ...]
#                        Region thickness (one value, or one per region)
#  --eta_from_xml ETA_FROM_XML
#                        Read eta/eta_surf, eta_vol, and thickness from an <oft><thincurr> XML file instead of the flags above
#  --coils_from_xml COILS_FROM_XML
#                        Read coil sets from an <oft><thincurr> XML file and add to the output mesh
#  --jumpers IDX [IDX ...]
#                        Nodeset indices to treat as jumpers instead of holes (0-based; negative indices count from the end). All other nodesets are holes
#  --jumper_range N [N ...]
#                        Range of nodeset indices to treat as jumpers, given as START [STOP] with Python slice semantics [START:STOP) (0-based; negative allowed;
#                        omit STOP to run to the end). Use '0' to mark every nodeset as a jumper
#```
#
# Script options for `combine` workflow:
#
#```shell
#usage: OFT_ThinCurr_mesh_tool combine [-h] --in_files IN_FILES [IN_FILES ...] [--out_file OUT_FILE] [--merge_regions]
#
#options:
#  -h, --help            show this help message and exit
#  --in_files IN_FILES [IN_FILES ...]
#                        Input mesh files (two or more)
#  --out_file OUT_FILE   Output mesh file
#  --merge_regions       Merge identical region IDs across inputs instead of keeping them distinct (default: keep distinct)
#```
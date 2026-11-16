# ThinCurr model reference

Detail for the `thincurr-model` skill: file formats, the Python and command-line
tools for each step, and troubleshooting. `SKILL.md` has the quick path.

## Native mesh file

HDF5, written by `write_ThinCurr_mesh()`, the converters, and `OFT_ThinCurr_holes`.
Indices are 1-based on disk; the Python readers/writers take and return 0-based.

| Dataset | Meaning |
|---|---|
| `mesh/R` [np,3] | Vertex positions [m] (a [np,2] list is treated as Z=0) |
| `mesh/LC` [nc,3] | Triangle vertex indices |
| `mesh/REG` [nc] | Region ID per cell (1..N) |
| `mesh/NODESETnnnn` | Legacy Cubit nodesets, used as holes (or jumpers) when `thincurr/holes` is absent |
| `thincurr/holes/HOLEnnnn` | Hole loops (ordered vertex lists) |
| `thincurr/CLOSURES` | Closure cells (one per closed surface) |
| `thincurr/HOLES_VERSION` | Written by `OFT_ThinCurr_holes` |
| `thincurr/junctions/...` | Present only in meshes split at junctions |
| `thincurr/ETA_SURF` or `ETA_VOL` + `THICKNESS` [nreg] | Resistivity [Ohm] or [Ohm-m] + [m] |
| `thincurr/coils/coilsetnnnn/...` | Coil sets |
| `thincurr/periodicity/NFP`, `PMAP` | Periodic (single field period) models |

Read everything back with `read_ThinCurr_mesh(filename)` (dictionary with `r`,
`lc`, `reg`, and a `thincurr` sub-dictionary).

## Building meshes

### gmsh (new geometry)

Use the Python API with the OCC kernel; it handles booleans (ports, cut-outs) and
curved surfaces. Rules that matter for ThinCurr:

- Mesh **surfaces only** (`gmsh.model.mesh.generate(2)`), and if you build a solid
  to get its faces, remove the volume and any faces you do not want
  (`occ.remove`). `getType(dim, tag)` distinguishes `Plane`/`Cylinder`/... faces.
- Put every surface you want in a **physical group**; the tag becomes the region
  ID. One region per distinct sheet or material, IDs 1..N.
- Surfaces that should be electrically connected must share their boundary curves
  (`occ.fragment` on touching parts makes the mesh conformal). Sheets that only
  overlap or nearly touch are **not** connected.
- Write format 2.2: `gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)`.
- Size control: `Mesh.MeshSizeMin`/`MeshSizeMax`, or `gmsh.model.mesh.setSize()` on
  points near small features.

```python
from OpenFUSIONToolkit._gmsh_cli import convert_gmsh_to_native   # CLI: OFT_convert_gmsh
convert_gmsh_to_native('part.msh', physical_regions=True)      # -> part.h5
```

`src/examples/ThinCurr/junctions/junction_meshes.py` builds sheets that meet on a
shared edge (a T-plate and a flanged tube) with the `geo` kernel.

### Structured grids (parametric surfaces)

```python
from OpenFUSIONToolkit.ThinCurr.meshing import build_triangles_from_grid, write_ThinCurr_mesh
r, lc, nodeset_n, nodeset_m = build_triangles_from_grid(grid)  # grid [n, m, 3]
write_ThinCurr_mesh('shell.h5', r, lc, reg, eta_vol=..., thickness=...)
```

`wrap_n`/`wrap_m` (default True) join the last row/column to the first; do not
repeat the end point in a wrapped direction. Leave a direction unwrapped for an
open edge (eg. a toroidal sector). Prefer running `OFT_ThinCurr_holes` over passing
the returned nodesets as holes.

For non-axisymmetric (stellarator) surfaces, `build_regcoil_grid(file, 'coil', ntheta, nphi)`
reads a REGCOIL output and `build_torus_bnorm_grid(file, ntheta, nphi)` reads a
B-normal mode file; `ThinCurr_periodic_toroid(r_grid, nfp, ntheta, nphi)` builds the
full or periodic mesh with holes and `write_to_file()` (see
`src/examples/ThinCurr/torus/torus_bnorm_mode.ipynb` and `REGCOIL/`).

`build_ThinCurr_dummy(center, size)` makes a tiny square sheet, useful as a
placeholder model that carries only coils.

### Combining and editing meshes

`OFT_ThinCurr_mesh_tool` (module `OpenFUSIONToolkit.ThinCurr._mesh_cli`):

```bash
OFT_ThinCurr_mesh_tool combine --in_files vv.h5 ports.h5 --out_file all.h5      # regions kept distinct
OFT_ThinCurr_mesh_tool modify --in_file sector.h5 --copies 7 --rotate z 45 --out_file full.h5
OFT_ThinCurr_mesh_tool modify --in_file m.h5 --remove_regions 3 --out_file m2.h5
OFT_ThinCurr_mesh_tool modify --in_file m.h5 --eta_vol 7e-7 --thickness 1e-3 --out_file m2.h5
OFT_ThinCurr_mesh_tool modify --in_file m.h5 --coils_from_xml coils.xml --circular_coils --out_file m2.h5
```

The tool remaps existing holes, closures, and per-region resistivity through
these edits (it warns when it has to drop resistivity), but the result is not
guaranteed to be topologically complete, eg. after removing a region a hole loop
passed through. Run `inspect_model.py` afterwards and rerun `OFT_ThinCurr_holes`
if the topology check fails. Combined or replicated parts are separate sheets:
coincident vertices are not merged.

## Holes, closures, jumpers, junctions

```bash
OFT_ThinCurr_holes --in_file mesh.h5 [--out_file X] [--plot_final] [--split_regions]
```
Python: `compute_homology('mesh.h5')` from `OpenFUSIONToolkit.ThinCurr._homology_cli`.
If the console script is not on `PATH`, `OFT_ThinCurr_holes.py` is in the install's
`python/` directory.

- **Holes** let current circulate around a boundary loop or handle (net current
  around a cylinder, poloidal/toroidal current on a torus). Count per connected
  piece: (boundary loops - 1) + 2 x genus.
- **Closures** remove the gauge freedom of a closed surface: one per closed piece.
- **Jumpers** (legacy Cubit nodesets) are cuts used to measure or drive current
  across a line, not holes. `setup_model(jumper_start=-K)` marks the last K nodesets
  as jumpers. Holes made by `OFT_ThinCurr_holes` never need this.
- **Junctions** are edges shared by three or more cells (a rib on a plate, a flange
  on a tube). `OFT_ThinCurr_holes` splits the mesh there and writes a
  `thincurr/junctions` group; ThinCurr enforces current conservation across the
  junction. Unsplit junction meshes are rejected. `--split_regions` also splits
  where two regions meet on an ordinary edge, so they can be insulated later.
  `setup_model(insulated_junctions=[3, (1, 2)])` insulates region 3 at every
  junction and regions 1 and 2 from each other, without a new mesh file. Junctions
  are not supported in periodic models.

## Resistivity and coils (XML)

```python
from OpenFUSIONToolkit.io import write_oft_xml
from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_Vcoil, ThinCurr_XML
xml = ThinCurr_XML()
xml.set_eta_vol([7.0E-7, 1.7E-8])          # [Ohm-m], one per region
xml.set_thickness([1.0E-3, 5.0E-3])        # [m]
# or: xml.set_eta([7.0E-4, 3.4E-6])        # surface resistivity [Ohm]
pf = ThinCurr_Icoil(name='PF1')            # fixed-current coil set
pf.add_subcoil(RZ=numpy.r_[1.2, 0.5])      # circular coil at (R, Z)
pf.add_subcoil(RZ=numpy.r_[1.2, -0.5], scale=-1.0)
xml.add_Icoil(pf)
ring = ThinCurr_Vcoil(name='ring', resistivity_per_length=2.2E-4, radius=0.005)  # passive conductor
ring.add_subcoil(RZ=numpy.r_[1.1, 0.25])
xml.add_Vcoil(ring)
write_oft_xml([xml], 'oft_in.xml')
```

Equivalent XML:

```xml
<oft>
  <thincurr>
    <eta_vol>7.0E-7, 1.7E-8</eta_vol>
    <thickness>1.0E-3, 5.0E-3</thickness>
    <icoils>
      <coil_set name="PF1">
        <coil>1.2, 0.5</coil>
        <coil scale="-1.0">1.2, -0.5</coil>
      </coil_set>
    </icoils>
  </thincurr>
</oft>
```

- `<eta>` (alias `<eta_surf>`) is surface resistivity [Ohm]; `<eta_vol>` [Ohm-m]
  needs `<thickness>` [m]. The list length must equal the largest region ID.
- A `coil` with two values is a circular loop at (R, Z) about the Z axis; general
  3D coils give `npts="N"` and N lines of `X, Y, Z`. All coils in a `coil_set` carry
  one current (times `scale`).
- `icoils` have prescribed current (drivers); `vcoils` are passive or
  voltage-driven conductors with `res_per_len` [Ohm/m] and `radius` [m]. vcoils
  are unknowns of the model and appear in its eigenmodes; icoils do not.
- `<sens_mask>` (one integer per region) removes regions from sensor signals.

Storing in the mesh file instead: pass `eta_surf=`, or `eta_vol=` and
`thickness=`, and `coil_sets=[...]` to `write_ThinCurr_mesh()`, or use
`OFT_ThinCurr_mesh_tool modify --eta_from_xml / --coils_from_xml`. The mesh file
takes precedence over the XML for each of resistivity and coils.

## Resolution

- `nelems` is roughly the number of interior vertices plus holes and vcoils
  (vertices on a boundary loop share one value). A dense L matrix needs 8*nelems^2 bytes
  (10k -> 0.8 GB, 20k -> 3.2 GB); above ~10k use HODLR in `thincurr-run`.
- Resolve the smallest feature that matters with several cells across it.
  To check convergence, rebuild at roughly half `h` and compare the leading decay
  times.
- Avoid very large element size ratios between adjacent regions; use gmsh size
  fields to grade instead.

## Troubleshooting

**Process exits with `ERROR: ... SUBROUTINE: tw_setup...`.** Fortran-side errors
abort the whole Python process (exit code 255 or similar), so they cannot be caught.
Run `inspect_model.py` first: it checks the topology before calling ThinCurr.

**`Error building hole mesh, path is not periodic` / "More nodesets than holes".**
A legacy mesh whose extra nodesets are jumpers: use `jumper_start` (`-1` for the
last nodeset). Or replace the nodesets by running `OFT_ThinCurr_holes`.

**`Hole points are not stored as an ordered loop` warning.** Legacy nodeset holes;
usually harmless. Rerun `OFT_ThinCurr_holes` to remove it.

**`Eta size mismatch`.** The `<eta>`/`<eta_vol>`/`<thickness>` list length is not
the largest region ID. `inspect_model.py` lists the region IDs.

**`eta_surf` negative in the report / "Cannot gather or infer surface
resisitivity".** No resistivity was given (or `eta_vol` without `thickness`).

**Unsplit junction mesh rejected.** Run `OFT_ThinCurr_holes`; it writes the split
mesh. Check that the gmsh surfaces meeting at the junction actually share the
edge (`occ.fragment`).

**Two sheets that should be connected carry independent currents.** They do not
share vertices along the seam. Rebuild with conformal geometry
(`occ.fragment`/`removeAllDuplicates`) rather than combining separate meshes.

**Too many holes found.** Small gaps or slivers in the geometry create extra
boundary loops. Check `Boundary loops` in the report and the plot.

## Source examples in this repository

| Topic | Location |
|---|---|
| Plate (Cubit), eigenmodes | `src/examples/ThinCurr/plate/` |
| Cylinder with jumper (Cubit), time domain | `src/examples/ThinCurr/cyl/` |
| Torus with ports, large model (HODLR) | `src/examples/ThinCurr/ports/` |
| Torus, mode-based surface (`ThinCurr_periodic_toroid`) | `src/examples/ThinCurr/torus/torus_bnorm_mode.ipynb` |
| REGCOIL surface | `src/examples/ThinCurr/REGCOIL/` |
| Junctions (gmsh, physical regions, insulation) | `src/examples/ThinCurr/junctions/` |
| Coil/XML creation | `src/examples/ThinCurr/xml/ex_oft_creation.ipynb` |

Python API: `src/python/OpenFUSIONToolkit/ThinCurr/meshing.py` (mesh I/O, grid
builders), `coils.py` (coils/XML), `_homology_cli.py`, `_mesh_cli.py`, and
`src/docs/ThinCurr/doc_thincurr_main.md` (numerics, holes, junctions).

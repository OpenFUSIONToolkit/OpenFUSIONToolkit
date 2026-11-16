---
name: thincurr-model
description: Build, convert, and verify ThinCurr (Open FUSION Toolkit) thin-wall electromagnetic models - triangular surface meshes of conducting structures plus per-region resistivity and coils. Use when asked to create a ThinCurr mesh or model of a vacuum vessel, wall, plate, shell, port, or coil set; to convert gmsh or Cubit meshes to ThinCurr's native format; to define holes/closures (OFT_ThinCurr_holes) or junctions between sheets; to set resistivity or thickness; or when ThinCurr setup fails on hole verification, "Eta size mismatch", or missing resistivity. Hands off to the thincurr-run skill for calculations.
---

# ThinCurr models

ThinCurr represents conducting structures as 2D sheets in 3D, discretized with a
triangular surface mesh. Building a correct model is most of the work; running it
is covered by the `thincurr-run` skill.

## What a model is (the handoff)

A model is a set of files that `ThinCurr.setup_model()` accepts and that
`scripts/inspect_model.py` reports as **Model OK**:

| Piece | Where it lives | Required |
|---|---|---|
| Triangular surface mesh, one region per distinct sheet/material | native `.h5` (`mesh/R`, `mesh/LC`, `mesh/REG`) | yes |
| Holes and closures (topology) | same `.h5`, written by `OFT_ThinCurr_holes` | yes, if the surface is not simply connected |
| Surface resistivity per region | `.h5` (`thincurr/ETA_SURF`, or `ETA_VOL`+`THICKNESS`) **or** XML `<eta>` | yes |
| Coils (`icoils` fixed current, `vcoils` passive/fixed voltage) | `.h5` (`thincurr/coils`) **or** XML | no |
| Jumper nodesets (legacy Cubit meshes) | `.h5` nodesets, selected with `jumper_start` | only if present |

When resistivity or coils are in **both** the mesh file and the XML, the mesh file
wins silently. Keep each in one place.

Report the model to the user (and to `thincurr-run`) as: mesh file, XML file (if
any), `jumper_start` (if not 0), `insulated_junctions` (if any), and the
`inspect_model.py` summary. The first four map one-to-one onto the `model` block
of a `thincurr-run` case file (`mesh_file`, `xml_file`, `jumper_start`,
`insulated_junctions`).

## Before anything else

```bash
python -c "from OpenFUSIONToolkit.ThinCurr import ThinCurr; print('ok')"
```

If that fails, set `OFT_ROOTPATH` to the install prefix (scripts here append
`$OFT_ROOTPATH/python` to `sys.path`) or add the install's `python/` directory to
`PYTHONPATH`. Inside the OFT repository the `oft-build` skill covers finding an
install. Needs `numpy`, `scipy`, `h5py`; `matplotlib` for plots; the `gmsh` Python
module for new geometry. Only one `OFT_env` may exist per Python process.

## Choosing a path

| Starting point | Path |
|---|---|
| New geometry described in words/dimensions | **gmsh Python API** (OCC kernel) -> `.msh` -> `convert_gmsh_to_native(..., physical_regions=True)` -> holes. See `examples/tube_with_port.py`. |
| Smooth parametric surface (torus, shaped shell, plasma boundary) | **Structured grid** -> `build_triangles_from_grid()` -> `write_ThinCurr_mesh()` -> holes. See `examples/torus_shell.py`. |
| Existing gmsh `.msh` | `OFT_convert_gmsh --in_file X.msh --physical_regions` -> holes |
| Existing Cubit/Exodus mesh | `OFT_convert_cubit --in_file X.g` -> holes (Cubit itself is not scriptable here) |
| Existing native `.h5` | Inspect it first; it may already be complete |
| Several parts, or N identical copies (eg. port extensions, coils sectors) | `OFT_ThinCurr_mesh_tool combine` / `modify --copies --rotate z ANGLE`, then holes |
| Surface from a REGCOIL or B-normal (mode) file | `build_regcoil_grid()` / `build_torus_bnorm_grid()` + `ThinCurr_periodic_toroid` (see `reference.md`) |

## Standard pipeline

1. **Geometry and mesh.** Give each distinct sheet or material its own region
   (gmsh physical surface). Choose the edge length `h` from the smallest feature
   that matters (port width, gap), not the overall size.
2. **Native format.** `.h5` via the converters or `write_ThinCurr_mesh()`.
3. **Topology.** Always run `OFT_ThinCurr_holes --in_file mesh.h5` (Python:
   `compute_homology('mesh.h5')`), which writes `mesh-homology.h5`. It finds holes
   for every boundary loop but one per piece, two per handle (torus), a closure per
   closed surface, and splits the mesh at junctions (edges shared by 3+ cells).
   Use the `-homology.h5` file from here on.
4. **Resistivity and coils.** Either an XML (`ThinCurr_XML`, see `reference.md`)
   or arrays passed to `write_ThinCurr_mesh()` before step 3 (they are carried
   through). Give exactly one value per region.
5. **Verify.**
   ```bash
   python <skill-dir>/scripts/inspect_model.py mesh-homology.h5 \
       --xml oft_in.xml --plot model.png
   ```
   Look at the plot. Fix every `ERROR`; read every `WARNING`. (`<skill-dir>` is this
   skill's directory, the one holding this `SKILL.md`.)

## Units

`eta` / `eta_surf` is **surface** resistivity in Ohm (Ohm per square) = volumetric
resistivity / thickness. Prefer giving `eta_vol` [Ohm-m] and `thickness` [m] so the
numbers are recognizable (stainless ~7e-7 Ohm-m, copper ~1.7e-8 Ohm-m). Positions
are in meters. The L/R time of a sheet scales as `mu0 * size / eta_surf`.

## Checking the model

`inspect_model.py` reports, and you should relay:

- points, cells, unknowns (`nelems`), regions with cell count and area;
- edge length min/mean/max against the features you meant to resolve;
- `Topology needs N hole(s), M closure(s)` (from the Euler characteristic) against
  what the file defines; it stops before ThinCurr setup when they cannot match,
  because a Fortran-side setup error aborts the whole Python process;
- where resistivity and coils come from, and the resolved `eta_surf` per region;
- `nelems` above ~10k: tell `thincurr-run` to use HODLR (dense L is 8*N^2 bytes).

Exit code 0 and `Model OK` means the model is ready for `thincurr-run`.

## When it fails

`reference.md` has the full list. The common ones: unsplit junctions (run
`OFT_ThinCurr_holes`); legacy Cubit meshes with more nodesets than holes (the extras
are jumpers, pass `--jumper-start -K`); `Eta size mismatch` (give one value for
every region ID up to the largest, including unused IDs, so keep IDs 1..N);
missing resistivity (negative `eta_surf` in the report).

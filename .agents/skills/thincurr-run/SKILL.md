---
name: thincurr-run
description: Run ThinCurr (Open FUSION Toolkit) thin-wall eddy-current calculations on an existing model - currently L/R eigenmode analysis (wall time constants, current decay times, and eddy-current mode structures), with dense or HODLR-compressed inductance matrices. Use when asked for the L/R time, wall time constant, eddy-current decay modes, or eigenvalues of a vacuum vessel, shell, plate, or other conducting structure, or for help with ThinCurr eigensolver runs, HODLR settings, or matrix cache files. Model building is the thincurr-model skill.
---

# Running ThinCurr

ThinCurr solves `L dI/dt + R I = V` for currents in thin conducting sheets.
This skill runs calculations on a model that already exists. If there is no model,
or `inspect_model.py` has not reported **Model OK** for it, use the
`thincurr-model` skill first.

Supported analyses:

| `analysis` | Computes |
|---|---|
| `eig` | Eigenmodes of `L I = tau R I`: decay times `tau` and current patterns of the slowest-decaying eddy currents |

## Before anything else

```bash
python -c "from OpenFUSIONToolkit.ThinCurr import ThinCurr; print('ok')"
```

If that fails, set `OFT_ROOTPATH` to the install prefix (the driver appends
`$OFT_ROOTPATH/python` to `sys.path`) or add the install's `python/` directory to
`PYTHONPATH`; inside the OFT repository the `oft-build` skill covers finding an
install. Needs `numpy`, `scipy`, `h5py`, and `matplotlib` for plots. Only one
`OFT_env`/`ThinCurr` per Python process: run each case in its own process.

## Running a case

Write a JSON case file and run the driver:

```bash
python <skill-dir>/scripts/run_thincurr.py case.json --output-dir out/
```

`<skill-dir>` is this skill's directory (the one holding this `SKILL.md`), wherever it
is installed.

A minimal eigenvalue case:

```json
{
  "model": {"mesh_file": "vessel-homology.h5", "xml_file": "oft_in.xml"},
  "analysis": "eig",
  "eig": {"neigs": 10},
  "output": {"prefix": "vessel"}
}
```

Paths in `model` are relative to the case file. Outputs go to `--output-dir`, else
`output.dir` (relative to the case file), else the working directory:

- `<prefix>_eig.json` - decay times, solver, timings, model size, resistivity;
- `<prefix>_modes.png` - current magnitude and direction of the leading modes;
- `<prefix>_plot/` - XDMF/HDF5 fields `J_01`, `J_02`, ... for ParaView/VisIt
  (or `OpenFUSIONToolkit.io.build_XDMF` + pyvista);
- `<prefix>_L.save` / `<prefix>_L_HODLR.save` - inductance cache, reused on
  reruns (validated against the mesh, coils, and HODLR settings, so it is safe).

Examples (they read meshes from `src/examples/ThinCurr/`, so they run as-is only
from the root of an OFT source checkout, with `--output-dir`):

- `examples/plate_eig.json` - square plate, dense, direct solver.
- `examples/cyl_eig.json` - cylinder with a legacy jumper nodeset (`jumper_start`).
- `examples/ports_eig_hodlr.json` - torus with ports, ~22k unknowns, HODLR + ARPACK
  (computing L takes about a minute on 2 threads; reruns use the cache).

## Choices the driver makes (override in the case file)

- **HODLR**: `"hodlr": "auto"` uses HODLR above 10,000 unknowns (dense L is
  8*N^2 bytes: 20k unknowns is 3.2 GB). Force with `true`/`false`.
- **Eigensolver**: `"eig": {"direct": "auto"}` uses LAPACK (all modes, robust) for
  dense models up to 4,000 unknowns, ARPACK otherwise. HODLR requires ARPACK.
- **Mode order**: slowest-decaying (largest `tau`) first.
- **Model options**: `jumper_start`, `insulated_junctions`, and an `eta` override
  (`eta_surf`, or `eta_vol` + `thickness`, one value per region) for scans without
  rebuilding the model.

## Checking and reporting the answer

Report:

- the model (files, unknowns, regions, resistivity per region) and whether HODLR
  was used;
- the leading decay times in ms, and what the slowest modes look like (from the
  PNG, or `region_dissipation_fraction` in the JSON for multi-region models);
- anything that limits confidence: mesh resolution not checked, HODLR tolerance,
  resistivity assumed rather than given.

Sanity checks before reporting:

- `tau` should scale like `mu0 * L_char / eta_surf` (`L_char` a characteristic
  size, eg. the minor radius). Off by orders of magnitude usually means a
  resistivity unit mix-up (`eta` is **surface** resistivity in Ohm, not Ohm-m).
- Symmetric structures give degenerate pairs (equal `tau`) for non-axisymmetric
  modes; a broken pair can indicate an asymmetric mesh.
- Topologically required modes (net current around a cylinder or torus, around a
  port) are usually the slowest; if they are missing, check the holes.
- For a quantitative answer, rerun with the model rebuilt at about half the edge
  length (thincurr-model) and compare the leading `tau`. With HODLR, also compare
  against a tighter `L_svd_tol` if the answer matters to better than a few percent.

## When it fails

See `reference.md`. A Fortran-side error (lines starting `ERROR:` followed by
`SUBROUTINE:`) aborts the Python process; most come from the model, so rerun
`inspect_model.py` from the thincurr-model skill.

## Going beyond the driver

Time-domain (`run_td`), frequency-response (`compute_freq_response`), sensors,
B-field reconstruction, and reduced models are not in the driver yet: write
the script directly against the API. `reference.md` has the call sequence and the
repository's notebook examples for each.

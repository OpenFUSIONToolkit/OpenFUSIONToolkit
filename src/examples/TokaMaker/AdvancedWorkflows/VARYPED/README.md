# VARYPED Equilibrium Scan

_Scanning pedestal pressure and edge current density with TokaMaker_

## Overview

This example uses `OpenFUSIONToolkit.TokaMaker.varyped` to generate free-boundary
equilibria with scaled pressure and current density profiles. It starts from a
DIII-D-like equilibrium and kinetic profiles, and demonstrates single-point and
small two-dimensional scans, output files, and profile comparisons.

The module also provides options for scans with an OMFIT VARYPED-compatible
grid, boundary constraints, tolerances, and output layout.

## File overview

| File | Role |
|---|---|
| `D3Dlike_VARYPED_ex.ipynb` | Main notebook; run from top to bottom |
| `DIIID_mesh.h5` | Link to the DIII-D TokaMaker mesh |
| `TkMkr_D3Dlike_Hmode_BSamp=1.0.geqdsk` | Link to the input equilibrium g-file |
| `TkMkr_D3Dlike_Hmode_BSamp=1.0.peqdsk` | Link to the input kinetic profiles p-file |

The data files are symbolic links to other TokaMaker examples in this repository.
Keep the repository directory structure intact so that the links resolve.

## Getting Started

1. Install or build OFT with the TokaMaker Python modules, and use a Python
   environment with Jupyter, NumPy, SciPy, Matplotlib, and h5py.
2. Open `D3Dlike_VARYPED_ex.ipynb` in JupyterLab or VS Code, with this folder as
   the notebook's working directory.
3. In the import cell, replace `/PATH_TO_OFT` with the OFT install root. The
   notebook adds its `python` subdirectory to `sys.path`. Alternatively, set
   `OFT_ROOTPATH` to the install root and uncomment the corresponding line in
   that cell.
4. Run the cells top to bottom. The notebook sets up the mesh, boundary and coil
   constraints, runs the scans, and plots and summarizes the results.

## Scan and outputs

- The notebook runs single scaling values of `1` and `1.5`, then a small scan
  using `[0.5, 1.5]` for both pressure and current scaling. Omitting
  `scaling_values` uses the default 11-by-11 grid from 0.5 to 1.5, which takes
  longer to run.
- The `1.5` single-point example enables file and results output. It creates a
  `varyped<shot>_<time>_<runid>` directory in the working directory containing
  g-files, p-files, a `results_*` file, and a `scan_status_*` file. The other
  example scans disable file output and use the returned results in memory.
- The plotting helpers accept the returned results; `plot_pfiles` and
  `plot_gfiles` instead read a scan directory containing a `results_*` file.

Input g-file directions are interpreted using the COCOS specified when reading
the file (the notebook uses the default COCOS 1). TokaMaker solves with positive
plasma current and `F0`, so coil constraints must match that solver convention.
Output g-files restore the input current and field directions and default to
COCOS 2; COCOS 7 can also be selected. Specify the output convention when
reading or plotting the files, since it is not stored in them.

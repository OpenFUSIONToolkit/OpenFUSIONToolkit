# NSTX-U Isoflux Control Notebooks

Main author: Matthew Parsons, PPPL

This example contains NSTX-U geometry files, Python libraries, and Jupyter notebooks for mesh generation, plasma shape creation, and shape control simulation workflows.

## File Overview

### Data & Geometry Files
- **NSTXU_mesh.h5**  
  HDF5 file containing the NSTX-U computational mesh.

- **NSTXU_geom.json**  
  JSON file describing machine geometry, boundaries, and parameters used by the mesh and shape generation tools.

---

### Python Modules
- **library_NSTXU_TM.py**  
  Tools and helper functions for NSTX-U tokamak modeling (coil definitions, equilibria utilities, etc.).

- **library_outputs.py**  
  Functions for exporting, formatting, and analyzing outputs from simulations.

- **library_PCS.py**  
  Plasma Control System utilities, including control matrices, actuator/control handling, and simulation helpers.

---

### Jupyter Notebooks
- **NSTXU_mesh_generator.ipynb**  
  Creates and visualizes NSTX-U meshes using the geometry and libraries.

- **NSTXU_shape_generator.ipynb**  
  Generates plasma shapes and explores shape parameterizations.

- **NSTXU_shape_control_simulator.ipynb**  
  Runs plasma shape-control simulations using PCS routines and NSTX-U geometry.

---

## Purpose

These files support workflows for:

- NSTX-U mesh generation  
- Plasma shape design and exploration  
- Shape-control simulation  
- Geometry and model processing  

The Python libraries provide shared utilities used across the notebooks.

---

## Getting Started

1. Ensure you have **Python 3.9+** and common scientific packages (`numpy`, `scipy`, `matplotlib`, `h5py`, etc.).
    - To run **NSTXU_shape_generator.ipynb** you also need additional packages (see notebook header).
2. Open the notebooks in **JupyterLab** or **VS Code**.
3. Run `NSTXU_shape_control_simulator.ipynb`!

---
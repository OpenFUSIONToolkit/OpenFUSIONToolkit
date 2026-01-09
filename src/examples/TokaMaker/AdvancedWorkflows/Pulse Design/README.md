# Pulse Design Notebooks

This folder contains examples of simple pulse design notebooks created in TokaMaker to illustrate the requisite workflows.

## File Overview

### Jupyter Notebooks
- **CUTE_pulse_ex.ipynb**  
  Shows a simple pulse design framework for the CUTE machine.

### Python Scripts
- **d3d_pulse_ex.py**  
  Shows a simple pulse design for General Atomics' DIII-D reactor.
---

## Purpose

These files support workflows for:

- Pulse design and optimization

---

## Getting Started
### Jupyter Notebooks
1. Ensure you have **Python 3.9+** and common scientific packages (`numpy`, `scipy`, `matplotlib`, `h5py`, etc.).
2. Open the notebooks in **JupyterLab** or **VS Code**.
3. Place all files in the same working directory so notebooks can import the libraries and load the geometry files.

### DIII-D Python Script
1. Ensure you have **Python 3.9+**, OFT, and common scientific packages (`numpy`, `scipy`, `matplotlib`, `h5py`, etc.).
2. Create a directory `163303` with subdirectories `eqs` and `profs`. Put geqdsk files for 163303 with naming convention `g163303-[time].eqdsk`, where `time` is in ms, in `eqs`. Put peqdsk files with naming convention `p163303.[time]` in the `profs` folder. Also add JSON files `ip.json` and `betap.json` in the `163303` directory and a mesh file `DIIID_mesh.h5`. To access this data, please contact the Columbia University Fusion Research Center.
3. Run the script using Python. If desired, set the variable `graph_shape=True` to view the plasma shape at each timestep.
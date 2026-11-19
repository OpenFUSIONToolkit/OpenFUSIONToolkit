# Breakdown and Start-up Design

This folder contains an example of tokamak breakdown and start-up design with TokaMaker. It covers the start-up of an ARC V3A discharge, from premagnetised coils in an empty vessel (from t = -2 s), through breakdown and the early plasma current rise, to a ramp of free-boundary equilibria reaching 0.57 MA at 1.8 s.

## File Overview

### Jupyter Notebooks
- **ARC_breakdown_startup_ex.ipynb**
  Walks through the workflow step by step, explaining the code and the physics as it goes:
  1. Machine model: the vacuum vessel and blanket are cut into segments, which TokaMaker treats as coils, so it computes every inductance and field map.
  2. Circuit model, integrated exactly between time steps.
  3. The start-up scenario: breakdown criteria (toroidal electric field and the Lloyd field-line criterion), the early plasma channel with its force balance and a limit on the stray field during burn-through, and the hardware limits (coil currents, supply voltages and slew rate).
  4. Stage 1, solved as one convex optimisation problem: premagnetisation, breakdown and the current rise to 100 kA, maximising the flux at the plasma boundary at the hand-off.
  5. A check of the breakdown against the avalanche physics (connection length and Townsend criterion).
  6. Stage 2: a ramp of TokaMaker free-boundary equilibria, coupled to the circuit by Picard iteration and optimised jointly with stage 1, with the plasma kept limited on the inboard wall.
  7. Figures of the whole start-up, and an animation written to `ARC_breakdown_startup.gif`.

### Python Modules
- **startup_model.py**
  Segmentation of the passive structure, the machine model (inductances, unit flux responses and coil limits from TokaMaker), the circuit model with exact time stepping, the early plasma channel, and field-line tracing with the Townsend criterion.
- **startup_optimizer.py**
  A small conic-program builder on the Clarabel solver, the stage-1 problem (`StartupProblem`), and the joint GS ramp (`GSRamp`).
- **startup_plots.py**
  Figures: machine, traces, breakdown maps, snapshots and the whole pulse.

### Other Files
- **requirements_Breakdown_Startup.txt**
  Python packages needed in addition to OFT (`clarabel` is the only one specific to this example).

The ARC V3A mesh and gEQDSK are read from the `ARC` example folder (`../../ARC/`). The gEQDSK is downloaded there from [Zenodo](https://zenodo.org/records/19498373) on first run if it is not already present.

---

## Purpose

These files support workflows for:

- Breakdown design with eddy currents in thick passive structures
- Coil and supply-voltage optimisation for premagnetisation, breakdown and the early current ramp
- Coupling a convex circuit optimisation to TokaMaker free-boundary equilibria

---

## Getting Started
1. Ensure you have **Python 3.9+**, OFT newer than the v26.9 release (for `TokaMaker.setup_model()`) with `OFT_ROOTPATH` set, and common scientific packages (`numpy`, `scipy`, `matplotlib`, `h5py`, etc.).
2. Install the conic solver Clarabel: `pip install clarabel` (or `pip install -r requirements_Breakdown_Startup.txt`).
3. Open `ARC_breakdown_startup_ex.ipynb` in **JupyterLab** or **VS Code** and run it from this folder, so it can import the modules and find `../../ARC/`. It takes about 12 minutes on a laptop.

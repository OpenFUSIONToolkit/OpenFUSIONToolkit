# Breakdown and Start-up

This folder contains a breakdown and start-up model built on TokaMaker. The model covers the passive structure as eddy-current segments, an exact circuit model, the breakdown criteria (Lloyd criterion, field-line tracing and the Townsend avalanche) and an early plasma channel. The folder uses the model in two ways, each with its own notebook:

| category | notebook | what it does |
|---|---|---|
| **optimization** | `optimization/ARC_breakdown_startup_ex.ipynb` | designs the start-up of an ARC V3A discharge: supply voltages from premagnetisation through breakdown to a 0.57 MA ramp of free-boundary equilibria, as convex optimization problems |
| **validation** | `validation/DIIID_breakdown_validation_ex.ipynb` | runs the same model feed-forward on DIII-D shot 200480 (the DIII-D reference discharge of Kim et al., Nucl. Fusion 66 (2026) 036043) from approximate trajectories of the discharge, and compares its loop voltage, plasma current and breakdown conditions with the measurements published there |

Both notebooks use the same machine model, circuit model, breakdown checks and figure style. The shared Python modules sit in this folder, and each notebook imports them with `sys.path.insert(0, '..')`.

## File Overview

### Shared Python Modules
- **startup_model.py**
  - Segmentation of the passive structure into segments that TokaMaker treats as coils.
  - The machine model: inductances, unit flux responses and coil limits from TokaMaker. Coil sets made of several elements, such as the DIII-D E coil, are supported.
  - The circuit model with exact time stepping, and `feed_forward`, which drives the structure with given coil currents (used for validation).
  - The early plasma channel.
  - Field-line tracing and the Townsend criterion.
- **startup_optimizer.py**
  - A small conic-program builder on the Clarabel solver.
  - The stage-1 problem (`StartupProblem`) and the joint GS ramp (`GSRamp`).
  - Used by the optimization notebook only.
- **startup_plots.py**
  - Figures: machine, traces, breakdown maps, snapshots, the whole pulse and the start-up animation.
  - `STYLE` and `COLORS`, the shared figure style.

### optimization/
- **ARC_breakdown_startup_ex.ipynb**
  Walks through the design step by step, explaining the code and the physics as it goes:
  1. Machine model: the vacuum vessel and blanket are cut into segments, which TokaMaker treats as coils, so it computes every inductance and field map.
  2. Circuit model, integrated exactly between time steps.
  3. The start-up scenario: the breakdown criteria (toroidal electric field and the Lloyd field-line criterion); the early plasma channel, with its force balance and a limit on the stray field during burn-through; and the hardware limits (coil currents, supply voltages and slew rate).
  4. Stage 1, solved as one convex optimization problem: premagnetisation, breakdown and the current rise to 100 kA, maximising the flux at the plasma boundary at the hand-off.
  5. A check of the breakdown against the avalanche physics (connection length and Townsend criterion).
  6. Stage 2: a ramp of TokaMaker free-boundary equilibria, coupled to the circuit by Picard iteration and optimized jointly with stage 1, with the plasma kept limited on the inboard wall.
  7. Figures of the whole start-up, and an animation written to `ARC_breakdown_startup.gif`.

  The ARC V3A mesh and gEQDSK are read from the `ARC` example folder (`../../../ARC/`). If the gEQDSK is not already there, it is downloaded from [Zenodo](https://zenodo.org/records/19498373) on first run.

### validation/
- **DIIID_breakdown_validation_ex.ipynb**
  A minimum working example of a feed-forward test. It runs with approximate trajectories of DIII-D 200480 (coil currents, plasma current and loop voltage as piecewise-linear keyframes in the notebook) and cites Kim et al. for every trace taken from the published figures. Nothing is optimized or fitted. It compares:
  1. the loop voltage at the inboard midplane, from the vessel eddy-current model, with the published measurement;
  2. the breakdown criteria (Lloyd, connection length and Townsend) with the breakdown and with the pressure scan of Kim et al.;
  3. the early-channel force balance and stray field, using the measured plasma current;
  4. the plasma current predicted with a prescribed electron temperature, with the published measurement.

  The model has **no plasma (burn-through) model**: the plasma current is replayed from the measurement where the loop voltage needs it, and the predicted plasma current depends on the assumed temperature.
- **g200480.00100**
  The EFIT01 gEQDSK of shot 200480 at 100 ms. The notebook uses its limiter and toroidal field.
  The coils and the vacuum vessel are read from `../../../DIIID/DIIID_geom.json`, whose E coil has one entry per turn (the graded solenoid winding sets the E coil's stray field at breakdown).

### Other Files
- **requirements_Breakdown_Startup.txt**
  Python packages needed in addition to OFT. `clarabel` is needed by the optimization notebook only.

---

## Purpose

These files support workflows for:

- Breakdown design with eddy currents in thick passive structures
- Coil and supply-voltage optimization for premagnetisation, breakdown and the early current ramp
- Coupling a convex circuit optimization to TokaMaker free-boundary equilibria
- Validating the electromagnetic and breakdown model against a measured discharge

---

## Getting Started
1. Ensure you have:
   - **Python 3.9+**;
   - OFT newer than the v26.9 release (for `TokaMaker.setup_model()`), with `OFT_ROOTPATH` set;
   - common scientific packages (`numpy`, `scipy`, `matplotlib`, `h5py`, etc.).
2. For the optimization notebook, install the conic solver Clarabel: `pip install clarabel` (or `pip install -r requirements_Breakdown_Startup.txt`).
3. Open a notebook in **JupyterLab** or **VS Code** and run it from its own folder, so that it can import the shared modules from `..`.
   - `optimization/ARC_breakdown_startup_ex.ipynb` takes about 12 minutes on a laptop.
   - `validation/DIIID_breakdown_validation_ex.ipynb` takes about 30 seconds.

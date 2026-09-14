# ARC V3A

_Mesh generation, free-boundary equilibrium, and vertical stability control for the ARC V3A tokamak_

## Overview

These examples are built on the ARC V3A digital assets described in Appendix A of Hillesheim et al, "Overview of the physics basis for the ARC fusion power plant", Journal of Plasma Physics (2026), [doi:10.1017/S0022377826101706](https://doi.org/10.1017/S0022377826101706). The underlying files are published on [Zenodo](https://zenodo.org/records/19498373) and each notebook downloads what it needs automatically if it is not already present locally.

## File overview

| File | Role |
|---|---|
| `ARCV3A_make_mesh.ipynb` | Builds the computational mesh from the device description (`doc_tMaker_ARC_ex1`) |
| `ARCV3A_baseline_ex.ipynb` | Reproduces the reference free-boundary equilibrium (`doc_tMaker_ARC_ex2`) |
| `ARCV3A_control_ex.ipynb` | Vertical stability and feedback control (`doc_tMaker_ARC_ex3`) |
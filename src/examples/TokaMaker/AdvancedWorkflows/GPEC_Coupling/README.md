# GPEC_Coupling

_Ideal kink and tearing stability of a DIII-D plasma during a pressure scan with TokaMaker + GPEC_

## Overview

This workflow reconstructs a DIII-D equilibrium (shot 147131 at 2300 ms) in TokaMaker, scales up
its pressure at fixed plasma current, shape and current density profile, and passes each
equilibrium to [GPEC](https://github.com/OpenFUSIONToolkit/GPEC). For each equilibrium, GPEC
computes:

- The $n=1$ ideal energy $\delta W$, which gives the ideal (no-wall) $\beta_N$ limit.
- The tearing stability index $\Delta'$ at the 2/1 surface.
- The Glasser correction $\Delta_{crit}$ at the same surface, using GPEC's SLAYER inner-layer model.

## File overview

| File | Role |
|---|---|
| `kink_tearing_stability.ipynb` | Main notebook. Run it from top to bottom |
| `g147131.02300_DIIID_KEFIT` | gEQDSK file from a kinetic EFIT reconstruction |
| `g147131.02300_DIIID_KEFIT_corrected.kin` | Kinetic profiles ($n_{e,i}$, $T_{e,i}$, $\omega_{E \times B}$) vs. $\hat{\psi}$ |
| `requirements_GPEC_Coupling.txt` | Python packages needed in addition to OFT |

The notebook also uses the DIII-D mesh from `examples/TokaMaker/DIIID`. At runtime it creates a
`gpec_runs/` directory with the gEQDSK, `gpec.toml` and `gpec.h5` files for each scan point.

## Setup

1. Install the Python packages and make OFT available, either by installing it or by setting
   `OFT_ROOTPATH` to your OFT install folder:
   ```bash
   pip install -r requirements_GPEC_Coupling.txt
   export OFT_ROOTPATH=/path/to/OpenFUSIONToolkit
   ```
2. Install [Julia](https://julialang.org/) and make sure `julia` is on your `PATH`.
3. Clone GPEC, install its dependencies, and install `PythonCall` (used by `juliacall`) in your
   default Julia environment:
   ```bash
   git clone https://github.com/OpenFUSIONToolkit/GPEC.git
   julia --project=/path/to/GPEC -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
   julia -e 'using Pkg; Pkg.add("PythonCall")'
   ```
4. Point the notebook to your copy of GPEC:
   ```bash
   export GPEC_PROJECT=/path/to/GPEC
   ```

## Running

```bash
jupyter lab kink_tearing_stability.ipynb
```

The first GPEC call includes Julia's compilation time, which can take several minutes. Each scan
point after that takes less than a minute.

This example is not part of the automatically tested example suite (`test_examples.sh`), since
it depends on external tools.

## Limitations

- The pressure is scaled by scaling the density, with fixed temperature profiles.
- $\Delta_{crit}$ uses an assumed $Z_{eff} = 1.5$ and perpendicular heat diffusivity
  $\chi_\perp = 1$ m²/s.
- The tearing metrics are not reliable close to the ideal limit, so they are only shown for
  $\beta < \beta_{lim}$.

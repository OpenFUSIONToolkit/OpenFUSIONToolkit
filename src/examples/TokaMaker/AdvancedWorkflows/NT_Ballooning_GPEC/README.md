# NT_Ballooning_GPEC

_Estimating the ballooning-limited pedestal height of a negative triangularity plasma with TokaMaker + GPEC_

## Overview

This workflow predicts how large an edge pedestal a **negative triangularity** (NT) plasma can
support before it becomes unstable to ideal **infinite-n ballooning modes**.

TokaMaker generates a family of bootstrap-consistent free-boundary equilibria that scan the
pedestal height at fixed shape and plasma current, and
[GPEC](https://github.com/OpenFUSIONToolkit/GPEC) evaluates the infinite-n ballooning stability
boundary for each one. The pedestal at which the edge first crosses that boundary is the
prediction. NT is the interesting case here because an NT edge has no access to the second
stable region, so the first ballooning boundary is the operative pedestal limit.

The example uses the same DIII-D NT discharge (shot 192185 at 2440 ms) and mesh as the other
DIII-D examples in `examples/TokaMaker/DIIID`, which it reads directly from that folder.

## File overview

| File | Role |
|---|---|
| `DIIID_NT_ballooning_GPEC_ex.ipynb` | Main notebook — run top to bottom |
| `compute_ballooning_boundary.jl` | Julia driver that calls GPEC for the ballooning boundary of each gEQDSK |
| `requirements_NT_Ballooning_GPEC.txt` | Python packages needed in addition to OFT |

The notebook creates a `ballooning_scan/` directory at runtime holding the gEQDSK written for
each equilibrium and the `ballooning_boundary_*.json` file GPEC writes back for it. These
outputs are not shipped with the example.

## Prerequisites

### 1. Python and OpenFUSIONToolkit

- The packages in `requirements_NT_Ballooning_GPEC.txt`:
  ```bash
  pip install -r requirements_NT_Ballooning_GPEC.txt
  ```
- A build of OFT providing `OpenFUSIONToolkit.TokaMaker`, including the `bootstrap` and `util`
  submodules. As in the other examples, set `OFT_ROOTPATH` to the OFT install root before
  launching Jupyter, or have OFT on your `PYTHONPATH`.

### 2. Julia and GPEC

GPEC is developed alongside OFT but distributed as a separate Julia package, so it is treated
here the way the TokaMaker+TORAX pulse design workflow treats TORAX: as an external module that
you install yourself and point the notebook at.

1. Install [Julia](https://julialang.org/) and make sure `julia` is on your `PATH`.
2. Obtain GPEC. Either add it to this repository as a submodule:
   ```bash
   # from the root of the OpenFUSIONToolkit repository
   git submodule add https://github.com/OpenFUSIONToolkit/GPEC.git external/GPEC
   git submodule update --init --recursive
   ```
   or clone it anywhere you like:
   ```bash
   git clone https://github.com/OpenFUSIONToolkit/GPEC.git
   ```
3. Instantiate the Julia project once, so its dependencies are installed:
   ```bash
   julia --project=/path/to/GPEC -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'
   ```
   `Pkg.resolve()` is what regenerates GPEC's `Manifest.toml` when it has fallen out of date with
   respect to its `Project.toml`. If you skip it and the notebook fails with

   ```
   ERROR: LoadError: ArgumentError: Package DelaunayTriangulation [...] is required but does
   not seem to be installed: - Run `Pkg.instantiate()` to install all recorded dependencies.
   ```

   (or the same message for any other package), run the command above — `Pkg.instantiate()` on
   its own cannot install a direct dependency that is missing from the manifest. GPEC's
   `Manifest.toml` is git-ignored, so regenerating it is a local operation only.

The notebook locates GPEC by looking, in order, at the `GPEC_PROJECT` environment variable,
`external/GPEC` in any parent directory (the submodule location above), and a `GPEC` clone
sitting next to the OpenFUSIONToolkit checkout. Setting `GPEC_PROJECT` explicitly is the
unambiguous option:

```bash
export GPEC_PROJECT=/path/to/GPEC
```

## Running

```bash
export OFT_ROOTPATH=/path/to/OpenFUSIONToolkit   # or have OFT on PYTHONPATH
export GPEC_PROJECT=/path/to/GPEC
jupyter lab DIIID_NT_ballooning_GPEC_ex.ipynb
```

Run the cells top to bottom. The two cells that shell out to Julia are the slow steps: the first
GPEC call pays Julia's JIT compilation cost (a few minutes), after which each equilibrium is
much faster. The TokaMaker equilibrium solves take roughly ten seconds each.

This example is not part of the automatically tested example suite (`test_examples.sh`), since
it depends on external tools.

## Method

1. **Baseline equilibrium.** The DIII-D NT gEQDSK is reproduced with an inverse solve, using the
   LCFS as an isoflux constraint and the EFIT coil currents as regularization targets, exactly as
   in the `DIIID_baseline_ex` example.
2. **Local NT proxy.** `eval_NT_H_proxy()` gives a cheap geometric read on edge ballooning
   stability — whether the maximum local shear on a near-edge surface sits in good or bad
   curvature ([A.O. Nelson et al., Nucl. Fusion 62, 096020 (2022)](https://doi.org/10.1088/1741-4326/ac8064)).
   It says whether the surface is favourably arranged, not how much pressure gradient it holds,
   which is what the GPEC calculation adds.
3. **Pedestal family.** `Hmode_profiles()` builds kinetic profiles at fixed density and fixed
   pedestal width, scanning the pedestal temperature with a stiff core ($T_{core}/T_{ped}$ fixed).
4. **Bootstrap-consistent equilibria.** `solve_with_bootstrap()` evaluates the Redl/Sauter
   bootstrap current from the kinetic profiles and reconstructs $FF'$ self-consistently, so the
   pedestal bootstrap current — and the edge magnetic shear it produces — feeds into the
   stability calculation.
5. **Ballooning boundary.** For each gEQDSK, GPEC's `LocalStability.ballooning_alpha_boundary`
   returns the experimental pressure-gradient drive $\alpha(\psi_N)$ and the first critical
   $\alpha_{crit}(\psi_N)$ at which the ballooning $\Delta'$ changes sign.
6. **Marginal pedestal.** The scan is reduced to $\max(\alpha/\alpha_{crit})$ over
   $\psi_N \geq 0.8$ per case, and interpolated to one to give the ballooning-limited pedestal
   height. A final equilibrium is solved at that pedestal and checked against GPEC.

## Limitations

- Only the **infinite-n** ballooning limit is imposed. A complete pedestal prediction (EPED-like)
  also needs the finite-n peeling-ballooning boundary and a pedestal width scaling; here the
  width is fixed by hand.
- The core profile shape, density profile, $Z_{eff}$, and pedestal width are prescribed rather
  than predicted, so the absolute pedestal height depends on those choices. The workflow is most
  useful for comparing configurations (shape, collisionality, current) under a fixed set of
  assumptions.

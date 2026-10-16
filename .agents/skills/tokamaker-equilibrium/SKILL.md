---
name: tokamaker-equilibrium
description: Compute axisymmetric Grad-Shafranov MHD equilibria with TokaMaker (Open FUSION Toolkit) from a target plasma boundary. Use when asked to solve, generate, or match a tokamak equilibrium - fixed-boundary from an LCFS shape, or free-boundary on an existing machine mesh where coil currents are solved for. Covers shaping moments (R0, a, kappa, delta, squareness), X-points and isoflux/saddle constraints, FF' and P' flux-function profiles, Ip/pax/beta targets, coil regularization and bounds, and gEQDSK output. Also use for TokaMaker convergence failures.
---

# TokaMaker equilibria

TokaMaker solves the Grad-Shafranov equation on an unstructured 2D triangular mesh.
Two modes:

- **Fixed boundary** — the LCFS is the mesh boundary. No coils, no machine geometry. Fast, and the right choice when only the plasma interior matters (profile studies, stability input, transport coupling).
- **Free boundary** — the plasma sits inside a machine mesh with coils and conductors. Given a target LCFS, TokaMaker solves for the coil currents that produce it (an "inverse" solve). Requires an existing mesh file.

## Before anything else

Check that the Python package imports. It needs a built/installed OFT:

```bash
python -c "from OpenFUSIONToolkit.TokaMaker import TokaMaker; print('ok')"
```

If that fails, set `OFT_ROOTPATH` to the install prefix (the driver appends
`$OFT_ROOTPATH/python` to `sys.path`), or add the install's `python/` directory to
`PYTHONPATH`. In this repository the install usually lives under `builds/*/install_release`.
`numpy` is required; `h5py` is required to read mesh files and `matplotlib` to plot.

Only **one** `OFT_env`/`TokaMaker` instance may exist per Python process. To run
several cases, use a separate process for each.

## Running a case

Write a JSON case file and run the driver:

```bash
python .agents/skills/tokamaker-equilibrium/scripts/run_equilibrium.py case.json
```

It writes `<prefix>_stats.json`, a gEQDSK file `g<prefix>`, and a flux-surface plot
`<prefix>_psi.png`, and prints the usual `print_info()` summary plus coil currents.

Start from the example closest to the request and edit it:

- `examples/fixed_boundary.json` — analytic LCFS shape, mesh generated on the fly.
- `examples/free_boundary.json` — ITER L-mode on the repository's `ITER_mesh.h5`, solving for coil currents against an isoflux + X-point target.

Both run as-is from the repository root.

A minimal fixed-boundary case:

```json
{
  "mode": "fixed",
  "F0": 0.10752,
  "lcfs": {"type": "shape", "R0": 0.42, "a": 0.15, "kappa": 1.4, "delta": 0.4},
  "targets": {"Ip": 120000.0, "Ip_ratio": 1.0},
  "profiles": {"ffp": {"type": "power", "alpha": 1.5, "gamma": 2.0},
               "pp":  {"type": "power", "alpha": 4.0, "gamma": 1.0}}
}
```

Free boundary adds `"mesh_file"`, a `"coils"` block, and usually X-points. See
`reference.md` for the full schema.

## Specifying the boundary

Three sources, set by `lcfs.type`:

| `type` | Give it | Use when |
|---|---|---|
| `shape` | `R0`, `Z0`, `a`, `kappa`, `delta` (or separate `*_upper`/`*_lower`) | The target is described by shaping parameters |
| `points` | `points: [[R, Z], ...]` | You have an explicit contour |
| `eqdsk` | `filename` | Matching the boundary of an existing equilibrium |

For a diverted plasma, add X-points. For a **single null** with a smooth opposite
side, give them explicitly and choose the moments so the Miller contour passes
through the X-point (`R_x = R0 - a*delta_lower`, `Z_x = Z0 - a*kappa_lower` for a
lower null):

```json
"lcfs": {"type": "shape", "R0": 6.2, "Z0": 0.5, "a": 2.0,
         "kappa_upper": 1.7, "delta_upper": 0.33,
         "kappa_lower": 1.95, "delta_lower": 0.5375,
         "xpoints": [[5.125, -3.4]]}
```

For a **double null**, use `"xpoints": "auto"`, which builds an X-pointed contour
with corners at both nulls and constrains both.

## Checking the answer

Always compare the achieved shape against what was asked. `print_info()` and
`<prefix>_stats.json` report `kappa`/`kappaU`/`kappaL`, `delta`/`deltaU`/`deltaL`,
`R_geo`, `a_geo`, `q_0`, `q_95`, `beta_n`, `l_i`, and `Topology`. Shape constraints
are soft, so a converged solve can still miss the target — say so rather than
reporting success on convergence alone. Report:

- whether it converged, and in how many iterations;
- achieved vs. requested shape and global parameters;
- for free boundary, the coil currents and whether any hit its bounds;
- `Topology` (Limited vs. Diverted) when a diverted shape was requested.

## When it does not converge

`reference.md` has the full list. The first four things to try: raise
`settings.maxits` to 60–100; soften the `init` seed shape (smaller `a`, lower
`kappa`, zero `delta`); lower `settings.urf` to 0.1; and for free boundary make sure
a VSC is defined, since elongated equilibria are vertically unstable and will not
converge without one.

## Going beyond the driver

The case file covers a single static solve. For anything else — forward solves from
known coil currents, H-mode or bootstrap profiles, reconstruction from diagnostics,
time-dependent VDEs, linear stability, pulse design — write the script directly
against the API. `reference.md` has the call sequence and a table of the
repository's worked notebook examples.

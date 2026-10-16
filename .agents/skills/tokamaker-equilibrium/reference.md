# TokaMaker equilibrium reference

Detail for the `tokamaker-equilibrium` skill: the case-file schema, the equivalent
direct-API calls, and troubleshooting. `SKILL.md` has the quick path.

## Case-file schema

All distances are in meters, currents in amperes, pressures in pascals, fields in tesla.
Paths inside a case file are resolved relative to the case file's own directory.

### Top level

| Key | Required | Meaning |
|---|---|---|
| `mode` | no | `"free"` (default) or `"fixed"` |
| `mesh_file` | free only | HDF5 mesh from `gs_Domain.save_mesh()` / the `*_mesh_ex.ipynb` examples |
| `mesh_dx` | no | Target edge length when generating a fixed-boundary mesh (default: minor diameter / 40) |
| `F0` | one of | Vacuum `F = B0 * R0` |
| `B0`, `R0_tor` | one of | Alternative to `F0`; `F0 = B0 * R0_tor` |
| `order` | no | FE order, 2–4 (default: 2) |
| `nthreads` | no | OpenMP threads (default: 2) |
| `lcfs` | yes | Boundary definition, see below |
| `targets` | yes | Global targets passed straight to `set_targets()` |
| `profiles` | no | `ffp` and `pp` flux functions |
| `coils` | free only | VSC, bounds, regularization |
| `init` | no | `init_psi()` overrides |
| `settings` | no | Any attribute of `mygs.settings` |
| `output` | no | Output selection |

### `lcfs`

Three ways to say what shape you want.

**`"type": "shape"`** — Miller-style shaping moments:

```json
{"type": "shape", "R0": 6.2, "Z0": 0.5, "a": 2.0,
 "kappa": 1.8, "delta": 0.45, "npts": 80, "xpoints": "none"}
```

Use `kappa`/`delta` for an up-down symmetric shape, or `kappa_upper`/`kappa_lower`
and `delta_upper`/`delta_lower` for an asymmetric one. `npts` is the number of
isoflux constraint points (default 80).

`xpoints` selects the boundary topology:

- `"none"` (default) — smooth Miller contour from `create_isoflux()`, no saddle constraints. Use for limited or fixed-boundary plasmas.
- `"auto"` / `"upper"` / `"lower"` — X-pointed contour from `create_isoflux_xpts()`, with saddle constraints at both / the upper / the lower X-point. The contour has a corner at **both** X-points regardless, so this suits double-null or near-double-null shapes.
- an explicit list `[[R, Z], ...]` — smooth Miller contour plus saddle constraints at exactly those points. This is the right choice for a single-null shape with a smooth opposite side: pick the moments so the Miller contour passes through the X-point, i.e. `R_x = R0 - a*delta_lower`, `Z_x = Z0 - a*kappa_lower` for a lower null.

With `"auto"`/`"upper"`/`"lower"` you can also set `r_inner_mid` and a `squareness`
block (`outer_upper`, `outer_lower`, `inner_upper`, `inner_lower`, all INTOR
squareness, default 0).

**`"type": "points"`** — an explicit contour, plus optional explicit `xpoints`:

```json
{"type": "points", "points": [[8.2, 0.41], [8.06, 1.46]], "xpoints": [[5.125, -3.4]]}
```

**`"type": "eqdsk"`** — take the boundary (`rzout`) from a gEQDSK file:

```json
{"type": "eqdsk", "filename": "g147131.02300", "xpoints": [[5.125, -3.4]]}
```

### `targets`

Passed verbatim as keyword arguments to `set_targets()`. Available keys: `Ip`,
`Ip_ratio`, `pax`, `estore`, `Dflux`, `R0`, `Z0`.

Pick **one** current-source constraint besides `Ip`:

- `Ip_ratio` — ratio of the FF' to P' contribution to `Ip`; roughly `1/beta_p - 1`.
- `pax` — on-axis pressure. Most natural when you know the target pressure.
- `estore` — stored energy.

`R0`/`Z0` target the magnetic axis position. In a free-boundary solve `R0` is used to
control plasma pressure and `Z0` drives the VSC, so use those instead of `pax` when
reproducing a known equilibrium from known coil currents (a "forward" solve).

### `profiles`

```json
{"ffp": {"type": "power", "alpha": 1.5, "gamma": 2.0},
 "pp":  {"type": "power", "alpha": 4.0, "gamma": 1.0}}
```

- `"power"` — `((1-psi_norm)^alpha)^gamma` via `create_power_flux_fun()`. L-mode-like defaults are peaked FF' (`alpha=1.5, gamma=2.0`) and broad P' (`alpha=4.0, gamma=1.0`).
- `"linterp"` — piecewise linear, give `x` and `y` arrays on normalized flux [0,1]. Use this to copy profiles out of a gEQDSK (`ffprim`, `pprime`; note the sign flip — TokaMaker's FF' and P' are the negative of the gEQDSK values).
- `"spline"` — cubic spline through `x`/`y` via `create_spline_flux_fun()`.
- `"default"` or omitted — leaves TokaMaker's built-in profile in place. Rarely what you want.

Profile amplitudes do not matter; they are rescaled to satisfy the global targets.
Only the shape matters. Forcing `y[-1] = 0` at the edge is usual.

### `coils` (free boundary)

```json
{"vsc": {"F9A": 1.0, "F9B": -1.0},
 "bounds": {"default": [-50e6, 50e6], "PF1": [-20e6, 20e6]},
 "regularization": {"default_weight": 1e-2,
                    "weights": {"CS1U": 2e-2},
                    "targets": {"PF6": 1.2e7},
                    "vsc_weight": 1e-2}}
```

- `vsc` — gains defining a vertical stability coil set. Nearly every elongated equilibrium needs one or it will not converge. Pair real up/down coils with opposite gains (`{"F9A": 1.0, "F9B": -1.0}`), or name a dedicated VS coil set (`{"VS": 1.0}`).
- `bounds` — hard current limits. `default` applies to every coil set; named entries override.
- `regularization` — one soft term per coil set, `weight * (I_coil - target)^2`. The driver adds every coil set automatically, so you only name the ones that differ. Higher weight ⇒ the coil stays nearer its target and the shape error grows.
- `vsc_weight` — weight on the virtual `#VSC` coil. Use a **small** value (~1e-2) to let it adjust vertical position, or a **large** value (~1e2) to disable it when you have defined a real VS coil set that already does that job.

Coil set names come from the mesh file. Print them with:

```python
print(list(mygs.coil_sets))
```

### `init`

`init_psi()` seeds the nonlinear iteration. Defaults are derived from the LCFS
moments: same `R0`/`Z0`, `0.8*a`, elongation relaxed 40% toward circular, zero
triangularity — deliberately rounder and smaller than the target, which converges
more reliably. Override any of `R0`, `Z0`, `a`, `kappa`, `delta`, or set
`{"uniform": true}` for a uniform current over the whole plasma region (the usual
choice for fixed boundary).

### `settings`

Any attribute of `mygs.settings`. The ones that matter here:

| Setting | Default | Use |
|---|---|---|
| `maxits` | 40 | Raise to 60–100 for hard cases |
| `urf` | 0.2 | Under-relaxation; raise toward 0.3–0.5 to speed up, lower to 0.1 to stabilize |
| `nl_tol` | 1e-6 | Convergence tolerance |
| `limited_only` | false | Skip X-point search — forces a limited solution |
| `lim_zmax` | 1e99 | Exclude divertor structures from the limiter point set |
| `rmin` | 0.0 | Minimum axis major radius, catches "lost" equilibria |
| `pm` | true | Per-iteration solver output |

`free_boundary` is set from `mode`; do not set it here.

### `output`

```json
{"prefix": "my_case", "dir": ".", "eqdsk": true, "plot": true, "stats": true,
 "run_info": "my scenario"}
```

Writes `<prefix>_stats.json`, `g<prefix>` (gEQDSK), and `<prefix>_psi.png`.

## Equivalent direct API calls

When a case file is too rigid, write the script directly. The order of operations
matters:

```python
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh
from OpenFUSIONToolkit.TokaMaker.util import create_isoflux, create_power_flux_fun

myOFT = OFT_env(nthreads=2)          # once per process, never twice
mygs = TokaMaker(myOFT)

mesh_pts, mesh_lc, mesh_reg, coil_dict, cond_dict = load_gs_mesh('machine_mesh.h5')
mygs.setup_mesh(mesh_pts, mesh_lc, mesh_reg)
mygs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)

mygs.settings.free_boundary = True   # must precede setup()
mygs.settings.maxits = 60
mygs.update_settings()
mygs.setup(order=2, F0=5.3*6.2)      # F0 = B0*R0

mygs.set_coil_vsc({'VS': 1.0})
mygs.set_coil_bounds({name: [-50.E6, 50.E6] for name in mygs.coil_sets})
mygs.set_coil_reg(reg_terms=[mygs.coil_reg_term({name: 1.0}, target=0.0, weight=1.E-2)
                             for name in mygs.coil_sets])

mygs.set_targets(Ip=15.6E6, pax=6.2E5)

lcfs = create_isoflux(60, 6.2, 0.5, 2.0, 1.7, 0.33, 1.95, 0.5375)
xpts = np.array([[5.125, -3.4]])
mygs.set_isoflux_constraints(np.vstack((lcfs, xpts)))
mygs.set_saddle_constraints(xpts)

mygs.set_profiles(ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
                  pp_prof=create_power_flux_fun(40, 4.0, 1.0))

mygs.init_psi(6.3, 0.5, 2.0, 1.4, 0.0)
_, nl_its = mygs.solve(return_its=True)   # raises ValueError if it fails
mygs.print_info()
```

Fixed boundary differs in three places: mesh the LCFS instead of loading a machine
mesh, set `free_boundary = False`, and skip all coil and shape constraints.

```python
from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
gs_mesh = gs_Domain()
gs_mesh.define_region('plasma', 0.015, 'plasma')
gs_mesh.add_polygon(lcfs, 'plasma')
mesh_pts, mesh_lc, _ = gs_mesh.build_mesh()
mygs.setup_mesh(mesh_pts, mesh_lc)
mygs.settings.free_boundary = False
mygs.setup(order=2, F0=0.10752)
mygs.set_targets(Ip=120.E3, Ip_ratio=1.0)
mygs.init_psi()
mygs.solve()
```

To run a second case in the same process, call `mygs.reset()` before loading the
new mesh. Creating a second `OFT_env` or `TokaMaker` in one process is not
supported — use a `multiprocessing.Process` per case, as the repository's tests do.

### Useful post-solve methods

| Call | Returns |
|---|---|
| `get_stats()` | Dictionary: `Ip`, `kappa`, `kappaU/L`, `delta`, `deltaU/L`, `R_geo`, `a_geo`, `vol`, `q_0`, `q_95`, `P_ax`, `W_MHD`, `beta_pol`, `beta_tor`, `beta_n`, `l_i`, `dflux`, `tflux` |
| `print_info()` | Same, formatted |
| `get_coil_currents()` | `(dict of currents [A], array)` |
| `get_profiles()` | `psi, F, F', P, P'` |
| `get_q()` | `psi, q, <R>/<1/R> averages, dl, rbounds, zbounds` |
| `get_xpoints()` | X-point locations; last entry is the active one |
| `trace_surf(psi)` | Contour of one normalized flux surface |
| `save_eqdsk(fn)` | gEQDSK output (COCOS 7 by default) |
| `save_ifile(fn)` | iFile output |
| `plot_machine`, `plot_psi`, `plot_constraints` | Matplotlib plots onto `(fig, ax)` |

## Troubleshooting

**`ValueError: Error in solve: Exceeded "maxits"`.** The nonlinear iteration did not
converge. In rough order of what to try:

1. Raise `settings.maxits` to 60–100. Some profile choices just converge slowly.
2. Check the targets are physical together. `Ip` with a large `pax` at low `F0` may have no solution. Try `Ip_ratio` instead of `pax` first to find a working ballpark, then read off `P_ax` from `print_info()`.
3. Soften `init`: smaller `a`, lower `kappa`, zero `delta`. A seed that is too close to a strongly shaped target often diverges.
4. Lower `settings.urf` to 0.1.
5. Free boundary only: make sure a VSC is defined. Elongated equilibria are vertically unstable and will not converge without one.

**Converged, but the shape is wrong.** The isoflux and saddle constraints are *soft* —
they compete with coil regularization. Lower `regularization.default_weight` (or the
weights on the coils nearest the region that is off) to let the coils work harder,
and check the resulting currents stay within `bounds`.

**Coil currents are implausibly large.** The opposite problem: raise
`default_weight`. There is no universally right value; it is a shape-error vs.
coil-current trade-off, and 1e-2 is only a starting point. If one coil set runs away,
raise its weight specifically.

**`Topology = Limited` when a diverted plasma was wanted.** The X-point fell outside
the plasma or no saddle constraint was set. Check the X-point location is inside the
vacuum vessel, that it is included in both the isoflux list and `set_saddle_constraints`
(the driver does both automatically), and that `settings.limited_only` is false.

**Axis drifts vertically / equilibrium is "lost".** Set `settings.rmin` above zero, add
or strengthen the VSC, and consider targeting `Z0` explicitly.

**`KeyError: Unknown coil "..."`.** Coil set names come from the mesh file's
`coil_dict`, and `coil_set` grouping may merge several mesh regions into one set.
Print `list(mygs.coil_sets)` after `setup_regions()`.

**Mesh resolution.** For a new geometry, re-run at roughly half the `mesh_dx` and
confirm `q_95`, `beta_n`, and `l_i` move by less than your tolerance.

## Source examples in this repository

| Topic | Notebook |
|---|---|
| Fixed boundary, analytic shape and from gEQDSK | `src/examples/TokaMaker/fixed_boundary/fixed_boundary_ex1.ipynb` |
| Free boundary, inverse and forward (ITER) | `src/examples/TokaMaker/ITER/ITER_baseline_ex.ipynb` |
| Shape helpers, X-points, squareness (DIII-D) | `src/examples/TokaMaker/DIIID/DIIID_shapehelper_ex.ipynb` |
| Mesh generation for a new machine | `src/examples/TokaMaker/*/*_mesh_ex.ipynb` |
| H-mode profiles, bootstrap | `src/examples/TokaMaker/ITER/ITER_Hmode_ex.ipynb` |
| Reconstruction from diagnostics | `src/examples/TokaMaker/ITER/ITER_reconstruction_ex.ipynb` |
| Time-dependent / VDE, stability | `src/examples/TokaMaker/CUTE/CUTE_VDE_ex.ipynb` |

Python API source: `src/python/OpenFUSIONToolkit/TokaMaker/_core.py` (solver class)
and `util.py` (shape and profile helpers).

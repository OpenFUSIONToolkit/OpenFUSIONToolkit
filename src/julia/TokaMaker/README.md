# TokaMaker.jl

Native Julia bindings for [TokaMaker](https://openfusiontoolkit.github.io/OpenFUSIONToolkit/),
the 2D axisymmetric Grad-Shafranov solver in the Open FUSION Toolkit (OFT).
This package binds **directly** to the same `liboftpy` shared library that
the Python `OpenFUSIONToolkit.TokaMaker` package uses, via Julia `ccall`.
There is no Python interpreter in the loop — Julia owns its FFI boundary.

## Status

This is a near-complete port mirroring the Python module at
`src/python/OpenFUSIONToolkit/TokaMaker/`. The core solver workflow and most
physics are functional; the remaining gaps are session I/O, equilibrium
copy/replace, a few utility/machine-file helpers, and test-suite parity. See
[`PORT_PLAN.md`](PORT_PLAN.md) for the milestone-by-milestone completion plan.

### Phase coverage

- [x] Phase 0: Scaffolding, library resolution, OFT env init
- [x] Phase 1: Settings struct + raw `ccall` wrappers (61/65 C symbols)
- [x] Phase 2: `TokaMaker` core type, `setup_mesh!`/`setup!`/`solve!`/`vac_solve!`
- [x] Phase 3: Profiles, constraints, coil management, util helpers
- [x] Phase 4: High-level meshing (`GsDomain` via Triangulate.jl) — no Cubit backend
- [x] Phase 5: Field eval, IO (eqdsk/ifile/mug), globals, q-profile, stats
- [x] Phase 6: Time-dependent + eigenvalue (`setup_td!`/`step_td!`/`set_psi_dt!`/`eig_td`/`eig_wall`)
- [x] Phase 7: Reconstruction (all 7 constraint types + `run_reconstruction!`)
- [x] Phase 8: Bootstrap current (Redl 2021 model + edge parameterization)
- [x] Phase 9: Plotting (Makie ext) — 4/5 (`plot_mesh` not yet ported)
- [~] Phase 10: Test suite mirroring `test_TokaMaker.py` (~25–35% parity)
- [x] Phase 11: ABI sync infra (`dev/check_abi_sync.jl`)

### Known gaps vs the Python module

Tracked against `src/python/OpenFUSIONToolkit/TokaMaker/`:

- **Session I/O**: `save_tokamaker`/`load_tokamaker` (HDF5) C wrappers unbound;
  equilibrium `save_TokaMaker` not ported.
- **Equilibrium lifecycle**: `copy_eq` / `replace_eq` / equilibrium `copy()` missing.
- **Bootstrap**: self-consistent `solve_with_bootstrap` loop not ported.
- **Utilities**: `read_mhdin`, `read_kfile`, `compute_forces_components`,
  `get_jphi_from_GS` missing.
- **Smaller helpers**: `set_resistivity`, `set_coil_current_dist`,
  `abspsi_to_normalized`/`psinorm_to_absolute`, `reconstruction.read_fit_in`,
  `plot_mesh`, `GsDomain.save_json`.
- **Tests**: spheromak, coil-field, LTX (eq/eig/stability), ITER bootstrap/Redl
  validation, EQDSK/COCOS round-trip, and concurrent multi-threaded solve not
  yet mirrored.

Run `julia --project=. dev/check_abi_sync.jl` to re-check FFI drift; currently
4 unbound C symbols (`oft_smesh_get`, `oft_vmesh_get`, `tokamaker_save_tokamaker`,
`tokamaker_load_tokamaker`).

## Library resolution

`TokaMaker.jl` looks for `liboftpy` in this order:

1. `ENV["OFT_LIBRARY_PATH"]/liboftpy.{so,dylib}`
2. `<package>/../../python/OpenFUSIONToolkit/liboftpy.{so,dylib}` (dev tree)
3. `<package>/../../../bin/liboftpy.{so,dylib}` (installed tree)
4. System search (`Libdl.find_library("oftpy")`)

Build OFT with `cmake -DOFT_BUILD_PYTHON=ON` (which also produces `liboftpy`)
before using this package end-to-end.

## Quick start

```julia
using TokaMaker

env = OFTEnv(nthreads=-1)
gs = TokaMaker(env)
setup_mesh!(gs; r=mesh_r, lc=mesh_lc, reg=mesh_reg)
setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
setup!(gs; order=2, F0=6.2)
set_profiles!(gs;
    ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
    pp_prof=create_power_flux_fun(40, 1.5, 2.0))
set_targets!(gs; Ip=15e6)
init_psi!(gs; r0=6.2, z0=0.0, a=2.0, kappa=1.7, delta=0.33)
set_isoflux_constraints!(gs, create_isoflux(32, 6.2, 0.0, 2.0, 1.7, 0.33))
eq = solve!(gs)
psi = get_psi(gs)
```

## License

LGPL-3.0-only, matching the parent project.

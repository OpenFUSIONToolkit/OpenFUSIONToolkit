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

- **Session I/O**: done (`save_tokamaker`/`load_tokamaker`, `replace_eq!`); see
  the HDF5 setup note below — liboftpy and HDF5.jl must share one libhdf5.
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

Run `julia --project=. dev/check_abi_sync.jl` to re-check FFI drift;
`tokamaker_save_tokamaker`/`tokamaker_load_tokamaker` are now bound, leaving
`oft_smesh_get`/`oft_vmesh_get` deferred (tessellated-mesh getters).

## Library resolution

`TokaMaker.jl` looks for `liboftpy` in this order:

1. `ENV["OFT_LIBRARY_PATH"]/liboftpy.{so,dylib}`
2. `<package>/../../python/OpenFUSIONToolkit/liboftpy.{so,dylib}` (dev tree)
3. `<package>/../../../bin/liboftpy.{so,dylib}` (installed tree)
4. System search (`Libdl.find_library("oftpy")`)

Build OFT with `cmake -DOFT_BUILD_PYTHON=ON` (which also produces `liboftpy`)
before using this package end-to-end.

## HDF5 setup (required for session I/O)

`TokaMaker.jl` uses HDF5.jl (for mesh I/O), and liboftpy also calls HDF5
internally (`save_tokamaker`/`load_tokamaker`, `save_mug`). By default HDF5.jl
loads its own libhdf5 (from HDF5_jll) while liboftpy links a *different* one
(from the OFT build). Two libhdf5 copies in one process split HDF5's global
state, so liboftpy's HDF5 writes silently produce an empty file and reads fail
("Failed to read FE order") — even on valid files.

Fix: point HDF5.jl at liboftpy's libhdf5 by adding a `LocalPreferences.toml` in
the package directory (it is gitignored — the path is machine-specific):

```toml
[HDF5]
libhdf5 = "<repo>/builds/.../lib/libhdf5.dylib"
libhdf5_hl = "<repo>/builds/.../lib/libhdf5.dylib"
```

Find the exact path with `otool -L <liboftpy.dylib> | grep libhdf5` (macOS) or
`ldd` (Linux), then use the unversioned `libhdf5.{dylib,so}` in that directory
for both keys (the OFT build omits a separate high-level library; HDF5.jl
resolves HL symbols lazily). Restart Julia after creating the file. Without it,
the session-I/O tests skip with a reminder.

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

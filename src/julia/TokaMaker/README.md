# TokaMaker.jl

Native Julia bindings for [TokaMaker](https://openfusiontoolkit.github.io/OpenFUSIONToolkit/),
the 2D axisymmetric Grad-Shafranov solver in the Open FUSION Toolkit (OFT).
This package binds **directly** to the same `liboftpy` shared library that
the Python `OpenFUSIONToolkit.TokaMaker` package uses, via Julia `ccall`.
There is no Python interpreter in the loop — Julia owns its FFI boundary.

## Status

This is an in-progress port mirroring the Python module at
`src/python/OpenFUSIONToolkit/TokaMaker/`. See
`src/julia/TokaMaker/CLAUDE_PORT_NOTES.md` (if present) for the porting plan.

### Phase coverage

- [x] Phase 0: Scaffolding, library resolution, OFT env init
- [x] Phase 1: Settings struct + 49 raw `ccall` wrappers
- [x] Phase 2: `TokaMaker` core type, `setup_mesh!`/`setup!`/`solve!`
- [x] Phase 3: Profiles, constraints, coil management, util helpers
- [ ] Phase 4: High-level meshing (GsDomain via Triangulate.jl)
- [ ] Phase 5: Field eval, IO (eqdsk/ifile/mug), globals, q-profile, stats
- [ ] Phase 6: Time-dependent + eigenvalue
- [ ] Phase 7: Reconstruction
- [ ] Phase 8: Bootstrap current
- [ ] Phase 9: Plotting (Makie ext)
- [ ] Phase 10: Test suite mirroring `test_TokaMaker.py`
- [ ] Phase 11: ABI sync infra

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

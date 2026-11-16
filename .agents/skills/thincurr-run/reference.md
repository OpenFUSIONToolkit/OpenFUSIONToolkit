# ThinCurr run reference

Detail for the `thincurr-run` skill: the case-file schema, the equivalent
direct-API calls, HODLR settings, and troubleshooting. `SKILL.md` has the quick path.

## Case-file schema

Units: SI throughout; decay times in seconds.

### Top level

| Key | Default | Meaning |
|---|---|---|
| `model` | required | Model files and setup options, see below |
| `analysis` | `"eig"` | Calculation to run |
| `nthreads` | 2 | OpenMP threads (`--nthreads` overrides) |
| `hodlr` | `"auto"` | `true`/`false`, or `"auto"` (HODLR above 10,000 unknowns) |
| `hodlr_options` | none | Keyword arguments for `set_hodlr_options()`, see "HODLR" |
| `cache` | `true` | Cache L in `<prefix>_L[_HODLR].save` in the output directory |
| `eig` | | Eigenvalue settings (when `analysis` is `"eig"`) |
| `output` | | Output selection |

### `model`

```json
{"mesh_file": "vv-homology.h5", "xml_file": "oft_in.xml",
 "jumper_start": 0, "insulated_junctions": [3, [1, 2]],
 "eta": {"eta_vol": [7.0e-7, 1.7e-8], "thickness": [1.0e-3, 5.0e-3]}}
```

| Key | Meaning |
|---|---|
| `mesh_file` | Native ThinCurr mesh (`.h5`), required |
| `xml_file` | XML with `<oft><thincurr>` resistivity/coils, if the mesh file does not hold them |
| `jumper_start` | Index of the first legacy nodeset that is a jumper (`-1`: the last one) |
| `insulated_junctions` | For meshes split at junctions: a region (insulated at all its junctions) or a `[a, b]` pair (insulated from each other) |
| `eta` | Resistivity override applied after setup with `set_eta_values()`: `eta_surf` [Ohm], or `eta_vol` [Ohm-m] with `thickness` [m], one value per region |

The driver stops if any region ends up with non-positive resistivity.

### `eig`

| Key | Default | Meaning |
|---|---|---|
| `neigs` | 5 | Number of modes (must be < number of unknowns) |
| `direct` | `"auto"` | `true`: LAPACK on the dense matrices (all modes, then truncated); `false`: ARPACK (iterative, needed for HODLR). `"auto"`: direct for dense models up to 4,000 unknowns |
| `save_modes` | 4 | Number of leading modes written to the plot files and PNG |

### `output`

| Key | Default | Meaning |
|---|---|---|
| `prefix` | case file name | Output file prefix |
| `dir` | working directory | Output directory, relative to the case file (`--output-dir` overrides) |
| `plot` | `true` | Write `<prefix>_modes.png` |
| `xdmf` | `true` | Write mode currents to `<prefix>_plot/` |
| `vectors` | `false` | Save all eigenvectors to `<prefix>_eigvecs.npy` [neigs, nelems] |

### Result file (`<prefix>_eig.json`)

```json
{"analysis": "eig",
 "model": {"nelems": 412, "nholes": 0, "nregs": 1, "n_vcoils": 0, "n_icoils": 1, "eta_surf": [1.257e-05], ...},
 "operators": {"hodlr": false, "L_cache": ".../plate_L.save", "time_L": 0.2, "time_R": 0.01},
 "results": {"neigs": 5, "solver": "direct", "time_eig": 0.04,
             "modes": [{"index": 1, "tau_s": 0.00973, "tau_ms": 9.73,
                        "region_dissipation_fraction": [...]}, ...]}}
```

`region_dissipation_fraction` (multi-region models only) is the share of each
mode's resistive dissipation in each region, a quick way to tell which part of the
structure a mode lives in.

## Equivalent direct API calls

```python
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr

myOFT = OFT_env(nthreads=2)              # once per process
tw = ThinCurr(myOFT)
tw.setup_model(mesh_file='vv-homology.h5', xml_filename='oft_in.xml')  # jumper_start=, insulated_junctions=
tw.setup_io(basepath='vv_plot/')         # where plot files go

tw.compute_Lmat(cache_file='vv_L.save')                          # dense, tw.Lmat [nelems, nelems]
# tw.compute_Lmat(use_hodlr=True, cache_file='vv_L_HODLR.save')  # large models
tw.compute_Rmat()                                                # sparse, tw.Rmat

eig_vals, eig_vecs = tw.get_eigs(10, False)   # (neigs,), (neigs, nelems); second arg: direct?
tw.save_current(eig_vecs[0, :], 'J_01')
tw.build_XDMF()
```

Other useful calls on a set-up model:

| Call | Returns / does |
|---|---|
| `tw.nelems`, `tw.nholes`, `tw.nregs`, `tw.n_icoils`, `tw.n_vcoils`, `tw.coils` | Model size and coil names/types |
| `tw.get_eta_values(include_eta_vol=True)` | `(eta_surf, eta_vol or None)` per region |
| `tw.set_eta_values(eta_surf=..., eta_vol=..., thickness=...)` | Change resistivity (then recompute R) |
| `tw.reconstruct_current(vec, centering='cell')` | Surface current [A/m] per cell (or `'vertex'`) |
| `tw.compute_Mcoil()` | Coil-to-model mutual inductances (drivers) |
| `tw.compute_Msensor('floops.loc')` | Model/coil-to-sensor couplings; sensors from `ThinCurr.sensor` (`Mirnov`, `flux_loop`, `save_sensors`) |
| `tw.compute_Bmat()` | Operators to reconstruct B from currents |
| `tw.run_td(dt, nsteps, coil_currs=..., sensor_obj=...)` | Time-domain run; `tw.plot_td()` for plot files |
| `tw.compute_freq_response(fdriver=..., freq=...)` | Frequency-domain response to a driver |
| `tw.build_reduced_model(basis)` | Reduced model (eg. from eigenvectors) |

## HODLR

Hierarchical off-diagonal low-rank compression of L. Memory and time scale far
better than the dense `O(N^2)`, at a controlled accuracy cost. Settings
(`"hodlr_options"`, passed to `set_hodlr_options()`):

| Option | Default | Effect |
|---|---|---|
| `target_size` | 1200 | Elements per leaf block |
| `L_svd_tol` | -1e-5 | SVD truncation tolerance for L blocks (negative: relative) |
| `L_aca_rel_tol` | 0.05 | ACA tolerance relative to `L_svd_tol` |
| `aca_min_its` | 20 | Minimum ACA iterations |
| `aca_sep_thresh` | [1.0, 1.1] | Separation ratio deciding ACA vs SVD (lowest level, other levels) |
| `pre_single`, `pre_coarse`, `pre_neigs` | off | Preconditioner options for time-domain/frequency runs |

Tighten `L_svd_tol` (eg. `-1e-6`) to check that HODLR is not changing the answer.
Changing any option invalidates the cache, which is then recomputed.

## Troubleshooting

**Process exits with `ERROR: ...` / `SUBROUTINE: ...` and no Python traceback.**
A Fortran-side error aborted the process. During setup this is a model problem:
run `inspect_model.py` (thincurr-model skill).

**`Surface resistivity is missing or non-positive for regions [...]`.** The model
has no resistivity for those regions. Fix it in the model, or set `model.eta`.

**ARPACK not available / ARPACK does not converge.** For dense models below a few
thousand unknowns use `"direct": true`. Otherwise ask for fewer modes or check the
model (a missing closure on a closed surface makes L singular).

**Decay times far from expectations.** Check the resistivity units first (`eta` is
surface resistivity, Ohm). Then holes: without the hole on a cylinder the net
azimuthal current mode is missing, and the slowest mode will be much faster.

**Out of memory computing L.** Use `"hodlr": true`.

**Slow L computation.** Expected: dense L costs `O(N^2)` and is cached; rerunning
with the same model reuses `<prefix>_L.save`. Raise `nthreads` for large models.

## Source examples in this repository

| Topic | Notebook |
|---|---|
| Eigenmodes of a plate | `src/examples/ThinCurr/plate/plate_eigenvalues.ipynb` |
| HODLR eigenmodes and time domain (torus with ports) | `src/examples/ThinCurr/ports/ports_HODLR.ipynb` |
| Junctions, insulated junctions, eigenmodes | `src/examples/ThinCurr/junctions/thincurr_junctions.ipynb` |
| Time domain, coils and sensors (cylinder) | `src/examples/ThinCurr/cyl/cyl_time_domain.ipynb`, `cyl_time_domain-jumpers.ipynb` |
| Forces | `src/examples/ThinCurr/cyl/cyl_forces.ipynb` |
| Frequency response, time domain (torus) | `src/examples/ThinCurr/torus/torus_freq_response.ipynb`, `torus_time_domain.ipynb` |
| Reduced models | `src/examples/ThinCurr/ports/ports_reduction.ipynb` |
| Fourier analysis of sensor signals | `src/examples/ThinCurr/torus/torus_bnorm_fourier.ipynb` |

Python API source: `src/python/OpenFUSIONToolkit/ThinCurr/_core.py`; numerics and
HODLR background: `src/docs/ThinCurr/doc_thincurr_main.md`.

# Phase 4: xMHD Jacobian on the quadrature-point cache, CSR position map, REFERENCE default

Branch `claude/xmhd-qcache-phase4`, starting from `claude/xmhd-qcache-phase3` (`b6f2c8b`). Phase 4 makes the
REFERENCE cache the default for the NL residual and ports the Jacobian assembly (`xmhd_build_ops`) to the cache,
on every p-multigrid level, with precomputed CSR scatter positions.

**Summary**
- `xmhd_qcache=2` (REFERENCE) is the new default. The Jacobian uses the cache when `xmhd_qcache=2`, on all levels
  on the finest mesh (finest and coarse polynomial levels). Modes 0 and 1 and levels on coarse meshes use direct
  basis evaluation.
- The linearized physics exists once: `xmhd_jac_qp` returns a per-quadrature-point coefficient tensor `D`
  coupling the value/gradient (Lagrange), value/curl (H(Curl)) and gradient (Grad(H^1)) components of every test
  and trial field. Both paths form local matrices from `D` with the same batched block contractions: one
  `dgemm` per row space and test-component group per cell. They differ only in how basis functions and fields
  are evaluated.
- The direct path is **round-off identical**, not bit-identical, to the original code: ≤5.9e-15 relative to
  the row maximum over all 37 matrices (33 configurations, every level), 1 thread, against a build of the
  unmodified Phase 3 branch. The summation order changed from per-qp outer products to batched contractions.
- Cache vs direct: ≤6.5e-14 relative to the row maximum on every level of every configuration. The exception
  is curved cells: 2.4e-12 (p=3) and 7e-13 (p=2) in the electron-temperature/H(Curl)-curl couplings. Products
  with random vectors agree to ≤1.7e-15. The CSR position map vs searching: ≤5.5e-16 (assembly order only).
- Run mode: linear and nonlinear iteration counts are identical at every step for NK (tet p2, tet p3 p-MG)
  and MFNK (tet p3 p-MG), for the original, new direct and new cache Jacobians (see below).
- Speed (Jacobian build per call, original → new default, fine level):
  - 1 thread: 3.0x (tet p2), 4.5x (tet p3), 5.3x (tet p4), 7.9x (hex p3), 4.2x (curved p3); coarse p-levels 3.0x.
  - 4 threads: 2.7x, 3.5-4.1x, 4.5-4.8x, 7.3x and 4.1x; coarse levels 2.9-3.5x; 2 MPI ranks 4.2x.
  - The CSR position map alone gives 1.3-2.0x. The batched contraction (shared by both paths) gives most of the
    rest: the new direct path is as fast as the cache path on fine tet levels, and the cache is 1.2-1.3x faster
    on hexes and coarse levels.
  - NK/MFNK time steps get only 9-11% faster: Jacobian time per run drops 2.7-3.6x, but the linear solve
    dominates. Iteration counts are identical.
- Memory: the position map costs 7-75 KB/cell for tets p2-p4 and 275 KB/cell for hexes p3, ~8-13% of RSS.
  Coarse-level caches add 0.5-3 MB.
- Multi-threaded runs need an OpenMP-aware BLAS (OFT's `build_libs.py` OpenBLAS) or `OPENBLAS_NUM_THREADS=1`.
  With Ubuntu's pthreads OpenBLAS the per-cell `dgemm` calls are threaded and the build stops scaling.
- New options: `xmhd_jac_posmap` (0/1, default 1), `xmhd_jac_posmap_maxmb` (memory cap per level and rank) and
  the experimental, default-off `xmhd_coarse_quad`. The coarse quadrature makes coarse-level builds 1.5-2.3x
  faster with unchanged iteration counts in the one run case tested.

## Design

### Pointwise kernel `xmhd_jac_qp`

`xmhd_jac_qp(u0, shp, det, he, eta_cell, solid, co, D)` contains all the linearized physics of the old
`xmhd_build_ops` quadrature loop after field reconstruction: floors, transport coefficients (Braginskii,
`eta_temp`/`res_profile`, shaping through `shp(4)` as in Phase 3), and the solid-cell branch. It also covers
every Jacobian block: B (ideal, resistive, Hall, electron inertia, hyper-resistivity), V (JxB, advection, pressure,
`kin`/`iso` viscosity, `nu_par`, `nu_gyro`, upwinding), n (advection, diffusion, upwinding, hyper-diffusivity),
T and Te (advection, compression, Ohmic and viscous heating, anisotropic conduction, equilibration, upwinding),
and the `xmhd_linear` variants. The local matrix of block `(f,g)` is

`K^{fg}_{rc} = sum_q sum_{a,b} E^f_{ra}(x_q) D_{(f,a),(g,b)}(x_q) E^g_{cb}(x_q)`

with basis components `a,b` in the order value then derivative: Lagrange `[phi, grad phi]` (4), H(Curl)
`[phi, curl phi]` (6), Grad(H^1) `[grad phi]` (3). Component `a` of field `f` is row `co(f)+a` of `D`, so `D` is
`ncomp x ncomp`, at most 43 x 43 (two-temperature with hyper terms). `D` includes `det*w`. Each term of the old
code maps to a few entries, eg. `-curl(phi).(V x psi)` becomes `-[V x]` in the (B curl, B value) block. The
Hall-term temperature coupling is written directly into the Te column for two-temperature models. This replaces
the old pointer aliasing `mtmp(1,7)%m=>mtmp(1,te_ind)%m`. Like `xmhd_resid_qp`, the kernel reads only module
settings and has no allocation or type-bound calls.

### Block contraction (`xmhd_jac_plan_setup`, `xmhd_jac_wform`, `xmhd_jac_contract`)

For each row space `s` (Lagrange, H(Curl), Grad(H^1)) and cell kind (fluid, solid), the active blocks of
`xmhd_mat_mask` are columns of one local matrix `C_s(nce_s, ntot)`. Solid cells use B blocks only. Blocks are
ordered by the test components they use: value only, value and derivative, derivative only. The rule is
static. Lagrange gradients are used for V-V, the diagonal n/T/Te blocks and the n-n2 couplings. H(Curl) values
are used for B-B, B-Bg and J2-J2. Per qp, `xmhd_jac_wform` contracts each used row of `D` with the trial functions:
`W(jc,(a,q)) = sum_b D((f,a),(g,b)) E^g_{jc,b}(x_q)`, a vectorized loop over `jc` with zero rows skipped. Per
cell, `xmhd_jac_contract` does two GEMMs per row space:

`C_s(:, 1:colv) = Ev_s · Wv_s^T` (value components, `K = nval*nqp`) and
`C_s(:, cold0+1:ntot) += Ed_s · Wd_s^T` (derivative components, `K = nder*nqp`).

`mtmp(f,g)%m` are pointer views of the column ranges of `C_s`, so scaling, BCs (Dirichlet rows, `vbc`
rotations), revolved-periodic packing and assembly reuse the existing code. That code is now one routine,
`xmhd_jac_finish_cell`, shared by both paths. BLAS `dgemm` (OpenBLAS small-matrix kernels) beat the
cache-blocked hand-written loop by 1.6-4.1x (see the A/B section), so `dgemm` is the default (`xmhd_jac_blas`, a
public module variable, also settable from the bench driver).

### Direct path (`xmhd_jac_direct`)

The direct path evaluates the basis per qp with the standard evaluators, exactly as before. The fine-level arrays
are shared with the interpolator through the `xmhd_*_rop` pointers. It copies them into the per-cell `E` arrays
and evaluates the fields with `fullin%interp`. It is used for `xmhd_qcache=0,1`, for levels on coarse meshes
(cell-centered `cc_interp` fields), and when the finest level is not the current one.

### Cache path (`xmhd_jac_qcache`)

Per cell:
1. gather the fine-level coefficients from the interpolator's local arrays and reconstruct V, n, T, [Te], [n2],
   B, J, [J2] at all qp with `fem_qcache_interp_cell` on the **finest-level** cache. Pack/unpack exactly as
   Phase 3 (`vals(6+3(k-1):8+3(k-1)) = grad(V_k)`, then `xmhd_interp_unpack` applies all scale factors), so there is
   no separate convention to get wrong;
2. build the test/trial functions of the **current level** at all qp with `fem_qcache_eval_cell` and gathered
   reference rows (`fem_qcache_cell_rows`), exact for all spaces, including H(Curl) faces on tets, where the
   canonical frame does not apply. Measured: basis rebuild is 5% (tets) to 11% (hexes) of the build, so the
   canonical tables plus `fem_qcache_mat_from_ref` (per-block transforms of every local matrix) would gain little;
3. `xmhd_jac_qp` per qp with `det*w` from `geom%det*geom%wts`, shaping factors from `geom%pt`, then
   `xmhd_jac_wform`, `xmhd_jac_contract`, `xmhd_jac_finish_cell` (DOF maps are columns of the cache's
   `dofs_*` tables).

### Coarse p-levels

Coarse polynomial levels live on the finest mesh and, as before, use the finest level's quadrature rule.
`xmhd_jcache_setup` builds a cache per coarse level holding that level's Lagrange/Grad(H^1)/H(Curl) tables, DOF
tables and element sizes (`he` depends on the level's order) **on the finest level's geometry** (shared pointer,
`own_geom=F`). Fields are reconstructed with the finest-level cache. `xmhd_set_ops` builds or validates the
finest cache, then the coarse caches, before the level loop. Rebuilding the finest cache (eg. shaping switched on,
or a mode change in the check) first frees the coarse caches that share its geometry. The quadrature order of
coarse levels is unchanged by default.

**Experimental `xmhd_coarse_quad=T`** (off by default): coarse polynomial levels integrate with their own
H(Curl) space's rule. In the cache path that level then owns a geometry for its rule plus a reconstruction cache
of the finest spaces on it (`qc%rec`). In the direct path only the rule changes. This changes the
preconditioner's coarse operators, not the nonlinear system: the finest level is bit-identical, coarse levels
differ by O(10%) (`pmg`: 1.0e-1 row-relative). Effect on build time and MG iterations: see below.

### CSR position map (`xmhd_jac_posmap_setup`, `xmhd_jac_add_native`)

`mat_add_values` searches each row's column list (`search_array`) for every local entry. The existing
`mloc%ind` cache only reuses the search within one cell, across blocks with the same pair of FE spaces. The
position map stores those positions for every cell and every (row space, column space) pair used by the blocks,
row-major, as int32 relative to the block row's start (`map(f,g)%ext(1,row)`). The fast add is
`M(ext(1,row)-1+rel) += K(jr,jc)`, row by row. Phase 1's scheme is kept: no lock for cells that share no DOFs
with another thread partition, per-field locks otherwise.
- Built on the first Jacobian build of a level (one search pass, 0.1-4 s at 1 thread, reported in `XMHD_PROF jac`), and
  rebuilt if the matrix's column list, rows or non-zeros change (`xmhd_jac_posmap_valid`). Both the direct and the
  cache path use it.
- Native matrices only. Other matrix types (PETSc) keep `mat_add_local`.
- `xmhd_jac_posmap=0` disables it, and `xmhd_jac_posmap_maxmb` skips it for levels that would need more memory
  per rank. Memory: see below.

### Lifecycle and selection

- The residual cache and the per-level Jacobian caches live in `xmhd_ops%qc`, the position maps in
  `xmhd_ops%pm`. All are freed by `xmhd_qcache_free` (end of `xmhd_run`/`xmhd_lin_run`, profile and check modes).
- `xmhd_build_ops` uses the cache path only if the level's cache is valid (mode, FE spaces, mesh, shaping points)
  and shares the geometry of the reconstruction cache, whose spaces must match the interpolator's. Otherwise it
  falls back to the direct path, so any caller is safe.
- Profiling counters: `XMHD_PROF jac level,qc_builds,posmap_mb,posmap_setup_s` per level and
  `XMHD_PROF jcache coarse_mem_mb,blas,posmap` (parsed by `run_bench.py`).

## Equivalence

Two independent checks.

**1. Against the original Jacobian, across builds.** `xmhd_bench` mode `qcheck` with `jac_dump='<prefix>'` writes
the Jacobian of every level, built at the check state, as CSR (`xmhd_jac_dump`). `jac_compare.py` runs the 32
Phase 3 configurations of `test_MUG_qcache.py` plus a 3-level p-MG case (`pmg_p4`), 33 in total, with any
executable, and compares dumps entrywise. The dumps were made with commit `d813eb2`, before the W-formation
tuning. The in-process checks below were rerun on the final build. The reference is a build of the unmodified Phase 3 branch plus the dump hook only, so its
Jacobian is the original code. The runs used 1 thread because the projected initial state is thread
nondeterministic: CG reductions make separate 2-thread runs differ by up to 1e-6 in row-relative Jacobian entries,
so cross-run comparisons need 1 thread. In-process comparisons are not affected. Largest difference relative
to the row maximum (`results/phase4_jdump_*.txt`):

| comparison | max over 37 matrices | worst case |
|---|---|---|
| new direct (`xmhd_qcache=0`) vs original | 5.9e-15 | `hex_3` |
| cache (`xmhd_qcache=2`) vs new direct | 2.4e-12 | `curved_p3` (all others ≤5.4e-14) |
| cache vs original | 2.4e-12 | `curved_p3` (all others ≤5.4e-14) |

The curved difference is in the Te ← B(curl) entries (~5e-12 of the entry itself). For p=3 H(Curl) on curved
cells, the cache rebuilds physical curls from folded logical components with the per-qp `cgop`, and the
evaluator combines them in a different order with cancellation. Phase 3's curved residual check (block-max
relative) is not sensitive enough to show this.

**2. In-process, every level (`test_MUG_qcache.py`).** After the residual checks, `xmhd_jac_check` builds all
levels' Jacobians at the same state in three configurations: (1) `xmhd_qcache=0` with search assembly (the
reference), (2) `xmhd_qcache=2` with the position map (the default), (3) `xmhd_qcache=0` with the position map.
It reports per level the entrywise difference relative to the row maximum and to the matrix maximum, and the
difference of products with two pseudo-random vectors. The test requires `rel_row ≤ 1e-11` and
`rel_matvec ≤ 1e-12`, with the expected number of levels. New cases: `test_qcache_pmg_p4` (3 levels, Hall,
hyper terms, shaping, solid cells), `test_qcache_curved_p2` (curved Jacobian; the p=3 curved Jacobian check needs
~8 GB, so `test_qcache_curved_p3` keeps the residual check only, `jac_check=F`) and
`test_qcache_pmg_coarse_quad`.

Largest values over all levels of each test (final build, `results/phase4_qcheck.log`, 35 passed in 103 s with 2
threads). "cache": `xmhd_qcache=2` + position map vs the reference; "posmap": `xmhd_qcache=0` + position map vs the
reference (assembly order only).

| test | levels | cache: rel_row | cache: rel_matvec | posmap: rel_row | posmap: rel_matvec |
|---|---|---|---|---|---|
| `test_qcache_physics_tet_p2[base]` | 1 | 9.6e-15 | 4.0e-16 | 2.7e-16 | 2.0e-17 |
| `test_qcache_physics_tet_p2[single_temp_hall]` | 1 | 2.3e-15 | 4.0e-16 | 2.8e-16 | 3.5e-17 |
| `test_qcache_physics_tet_p2[two_temp_hall_equil]` | 1 | 7.3e-15 | 4.3e-16 | 2.8e-16 | 3.2e-17 |
| `test_qcache_physics_tet_p2[iso_gyro]` | 1 | 6.9e-15 | 4.2e-16 | 3.6e-16 | 2.1e-17 |
| `test_qcache_physics_tet_p2[kin_gyro_nopar]` | 1 | 7.8e-15 | 4.2e-16 | 3.4e-16 | 5.8e-18 |
| `test_qcache_physics_tet_p2[brag_two_temp]` | 1 | 2.6e-15 | 4.1e-16 | 2.8e-16 | 7.6e-18 |
| `test_qcache_physics_tet_p2[brag_single_temp]` | 1 | 2.6e-15 | 4.0e-16 | 3.6e-16 | 1.9e-17 |
| `test_qcache_physics_tet_p2[hyper]` | 1 | 9.0e-15 | 3.9e-16 | 2.7e-16 | 1.5e-17 |
| `test_qcache_physics_tet_p2[hyper_centered]` | 1 | 9.7e-15 | 4.1e-16 | 3.0e-16 | 2.8e-17 |
| `test_qcache_physics_tet_p2[shaping_eta_temp]` | 1 | 6.5e-15 | 4.3e-16 | 2.6e-16 | 2.0e-17 |
| `test_qcache_physics_tet_p2[eta_profile]` | 1 | 1.0e-14 | 4.0e-16 | 2.6e-16 | 3.1e-17 |
| `test_qcache_physics_tet_p2[no_upwind_iso]` | 1 | 9.2e-15 | 3.9e-16 | 5.5e-16 | 4.0e-17 |
| `test_qcache_physics_tet_p2[no_advec_den]` | 1 | 8.8e-15 | 3.9e-16 | 3.6e-16 | 1.6e-17 |
| `test_qcache_physics_tet_p2[no_temp]` | 1 | 5.1e-16 | 4.7e-17 | 2.8e-16 | 7.1e-18 |
| `test_qcache_physics_tet_p2[floors_solid]` | 1 | 5.2e-15 | 3.8e-16 | 2.8e-16 | 7.5e-17 |
| `test_qcache_physics_tet_p3[base]` | 1 | 2.2e-14 | 4.2e-16 | 2.6e-16 | 6.2e-18 |
| `test_qcache_physics_tet_p3[single_temp_hall]` | 1 | 5.9e-15 | 4.2e-16 | 4.0e-16 | 1.3e-17 |
| `test_qcache_physics_tet_p3[hyper]` | 1 | 2.3e-14 | 4.1e-16 | 2.6e-16 | 1.0e-17 |
| `test_qcache_physics_tet_p3[shaping_eta_temp]` | 1 | 2.7e-14 | 4.1e-16 | 3.3e-16 | 1.1e-17 |
| `test_qcache_physics_tet_p3[floors_solid]` | 1 | 2.4e-14 | 4.1e-16 | 3.1e-16 | 1.2e-17 |
| `test_qcache_tet[2]` | 1 | 9.5e-15 | 4.2e-16 | 2.8e-16 | 4.5e-17 |
| `test_qcache_tet[3]` | 1 | 2.5e-14 | 4.4e-16 | 3.1e-16 | 1.2e-17 |
| `test_qcache_tet[4]` | 1 | 3.1e-14 | 6.5e-16 | 3.2e-16 | 8.5e-18 |
| `test_qcache_hex[2]` | 1 | 3.4e-15 | 4.3e-16 | 0.0e+00 | 0.0e+00 |
| `test_qcache_hex[3]` | 1 | 3.2e-14 | 5.6e-16 | 0.0e+00 | 0.0e+00 |
| `test_qcache_curved_p3` | - | (residual only, `jac_check=F`) | | | |
| `test_qcache_curved_p2` | 1 | 7.0e-13 | 1.7e-15 | 4.1e-16 | 9.2e-18 |
| `test_qcache_vbc[norm]` | 1 | 1.0e-14 | 3.7e-16 | 3.5e-16 | 1.7e-17 |
| `test_qcache_vbc[tang]` | 1 | 6.7e-15 | 3.7e-16 | 2.2e-16 | 2.0e-17 |
| `test_qcache_revolved[tet]` | 1 | 3.3e-14 | 4.6e-16 | 2.7e-16 | 8.1e-18 |
| `test_qcache_revolved[hex]` | 1 | 4.2e-14 | 3.8e-16 | 2.6e-16 | 2.5e-18 |
| `test_qcache_pmg` | 2 | 3.0e-14 | 4.2e-16 | 3.3e-16 | 1.9e-17 |
| `test_qcache_pmg_p4` | 3 | 2.3e-14 | 7.0e-16 | 4.1e-16 | 2.3e-17 |
| `test_qcache_pmg_coarse_quad` | 2 | 2.0e-14 | 4.4e-16 | 3.5e-16 | 3.8e-17 |
| `test_qcache_mpi2` | 1 | 6.5e-14 | 4.0e-16 | 3.4e-16 | 7.0e-18 |
| **max** | | 7.0e-13 | 1.7e-15 | 5.5e-16 | 7.5e-17 |

## Validation

- **Regression** (`results/phase4_regression.log`, new defaults: residual and Jacobian on the cache, position
  map), from `builds/build_release/tests` with `OMP_NUM_THREADS=2`: `test_alfven.py`, `test_sound.py` with
  `-k "test_nl_r1_p2 or test_lin_r1_p2 or test_hex_nl_r1_p2 or (test_sound and test_nl_r1_p3) or visc or gyro or
  hall or qcache"` (incl. 2-rank NK/MFNK and linear runs, sound p3 MG), `test_MUG_visc.py`, `test_MUG_gyro.py`,
  `test_MUG_hall.py`, `test_MUG_qcache.py`: **77 passed, 20 skipped** (PETSc and excluded parameter combinations),
  0 failed. The known p-MG failures (alfven `test_nl_r1_p3`/`test_lin_r1_p3`) are outside this subset. This ran
  on commit `695c937`, before the experimental coarse quadrature was added; the full `test_MUG_qcache.py` was rerun
  on the final build (35 passed).
- **Runtime checks**: Debug build (`-fcheck=all,no-array-temps`) of `xmhd_bench`, cache-path check runs of
  `pmg_p4` (3 levels), `floors_solid`, `revolved[tet]`, `vbc[norm]`, `hex[2]`: no runtime errors, results as
  in release.
- **Valgrind** memcheck (release build, 1 thread; p-MG p3→p2, Hall, hyper terms, gyroviscosity, shaping, solid
  cells, floor sources; residual and Jacobian checks): 0 errors (`results/phase4_valgrind.txt`). No leak record
  points into the new code (the largest, 0.74 MB, is the driver's mass-matrix graph).

## Performance A/B

Machine: Intel Xeon (Cooper Lake) cloud VM, 4 vCPUs (shared host, ±15-30% run-to-run noise). Ubuntu 24.04 system
libraries (`oft-build` section 5): GCC 13 `-O2` + LTO, OpenMPI 4.1, SuperLU, METIS, and OpenBLAS 0.3.26. The
`run_bench.py` default matrix is the same as in Phases 0-3, with ILU block-Jacobi. Script `phase4_ab.sh`, tables from
`phase4_tables.py`, configurations interleaved, 2 repetitions each, medians shown:

| config | executable | residual | Jacobian | assembly |
|---|---|---|---|---|
| `p3` | Phase 3 branch (`b6f2c8b` + dump hook) | cache (`xmhd_qcache=2`) | original | search (`mloc%ind` within a cell) |
| `q0pm0` | this branch | direct | new direct path | search |
| `q0pm1` | this branch | direct | new direct path | position map |
| `q2pm0` | this branch | cache | cache path | search |
| `q2pm1` | this branch | cache | cache path | position map (**default**) |
| `q2pm1h` | this branch | cache | cache path, hand-written GEMM loop | position map |

The 1-thread runs use profile mode with `--prof-reps 2,1,0`: two Jacobian builds per run, after the first build
in `xmhd_run`, which also builds the caches and position maps (their setup is reported separately). The 4-thread
runs use the full default matrix (profile, 2 MPI ranks, NK/MFNK run cases).

**BLAS threading.** `dgemm` is called from every OpenMP thread. Ubuntu's default OpenBLAS is the pthreads build,
which threads each of those small GEMMs itself. With it the new build did not scale with OpenMP threads: tet p4
took 2.33 s at 1 thread and 2.21 s at 4 threads; hex p3 took 3.06 s and 2.10 s. With `OPENBLAS_NUM_THREADS=1` it
took 0.71 s and 0.94 s (`results/p4blas_*.json`). OFT's own `build_libs.py` builds OpenBLAS with `USE_OPENMP=1`,
which runs single-threaded inside OpenMP regions, so the 4-thread runs below use Ubuntu's OpenMP variant
(`libopenblas0-openmp`, selected with update-alternatives). The first, pthreads 4-thread runs are kept as
`results/p4ab_t4pthreads_*`. The 1-thread runs were not affected: the small-kernel path is single-threaded and
the profile shows no BLAS threads; spot checks with the OpenMP variant agree to within 8%.

### 1 thread: Jacobian build time per call (prof_reps=2,1,0)

| case | level | p3 [s] | q0pm0 [s] | q0pm1 [s] | q2pm0 [s] | q2pm1 [s] | q2pm1h [s] | p3/q2pm1 | q2pm0/q2pm1 (posmap) | p3/q2pm0 (cache+GEMM) | q0pm1/q2pm1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `tet_p2` | fine | 3.924 | 2.685 | 1.344 | 2.595 | 1.306 | 2.069 | 3.00 | 1.99 | 1.51 | 1.03 |
| `tet_p3` | fine | 23.105 | 10.407 | 5.202 | 10.258 | 5.141 | 11.535 | 4.49 | 2.00 | 2.25 | 1.01 |
| `tet_p3_mg` | fine | 23.412 | 10.300 | 5.176 | 10.280 | 5.104 | 11.302 | 4.59 | 2.01 | 2.28 | 1.01 |
| `tet_p3_mg` | coarse | 7.651 | 4.176 | 3.016 | 3.773 | 2.558 | 4.414 | 2.99 | 1.47 | 2.03 | 1.18 |
| `tet_p4` | fine | 12.451 | 4.592 | 2.368 | 4.372 | 2.333 | 5.976 | 5.34 | 1.87 | 2.85 | 1.01 |
| `hex_p3` | fine | 24.064 | 4.351 | 3.661 | 3.980 | 3.056 | 12.333 | 7.88 | 1.30 | 6.05 | 1.20 |
| `curved_p3_mg` | fine | 31.189 | 14.470 | 7.738 | 14.375 | 7.417 | 15.646 | 4.21 | 1.94 | 2.17 | 1.04 |
| `curved_p3_mg` | coarse | 10.651 | 6.440 | 4.698 | 5.439 | 3.562 | 6.278 | 2.99 | 1.53 | 1.96 | 1.32 |

Repetitions: p3: 2/2/2/2/2/2, q0pm0: 2/2/2/2/2/2, q0pm1: 2/2/2/2/2/2, q2pm0: 2/2/2/2/2/2, q2pm1: 2/2/2/2/2/2, q2pm1h: 2/2/2/2/2/2

- **Overall** (`p3/q2pm1`): fine-level Jacobian builds are **3.0x (tet p2), 4.5x (tet p3), 5.3x (tet p4), 7.9x
  (hex p3) and 4.2x (curved p3) faster**, and coarse p-levels 3.0x, at 1 thread.
- **Position map** (`q2pm0/q2pm1`): 1.9-2.0x on tets, 1.3x on hexes, 1.5x on coarse levels. Searching the column
  lists was half of the remaining build time. Setup is one search pass per level and run (0.2-4 s at 1 thread,
  below).
- **Cache + batched GEMM** (`p3/q2pm0`, both searching): 2.2-2.9x on tets, 6.1x on hexes, 2.0x on coarse levels.
  Most of the gain is the batched contraction, which the new direct path shares: `q0pm1/q2pm1` is 1.0 on fine
  tet levels. The cache path is 1.2x faster than the direct path on hexes and 1.2-1.3x on coarse levels, where
  the direct path evaluates the fine-level basis at every qp to reconstruct the fields.
- **BLAS vs hand-written loop** (`q2pm1h/q2pm1`): `dgemm` (OpenBLAS AVX-512 small-matrix kernels) is 1.6x (tet p2)
  to 4.1x (hex p3) faster than the cache-blocked Fortran loop built with generic `-O2`. `dgemm` is the default.

### 4 threads: Jacobian build time per call (default matrix, OpenMP OpenBLAS)

| case | level | p3 [s] | q2pm1 [s] | p3/q2pm1 |
|---|---|---|---|---|
| `tet_p2` | fine | 1.116 | 0.422 | 2.65 |
| `tet_p3` | fine | 6.223 | 1.768 | 3.52 |
| `tet_p3_mg` | fine | 6.480 | 1.584 | 4.09 |
| `tet_p3_mg` | coarse | 2.107 | 0.740 | 2.85 |
| `tet_p4` | fine | 3.341 | 0.745 | 4.48 |
| `tet_p4_mg` | fine | 3.383 | 0.712 | 4.75 |
| `tet_p4_mg` | coarse | 1.727 | 0.493 | 3.51 |
| `hex_p3` | fine | 6.416 | 0.882 | 7.28 |
| `curved_p3_mg` | fine | 8.807 | 2.152 | 4.09 |
| `curved_p3_mg` | coarse | 3.031 | 1.058 | 2.86 |
| `tet_p3_mg_omp1` | fine | 22.608 | 4.966 | 4.55 |
| `tet_p3_mg_omp1` | coarse | 7.294 | 2.451 | 2.98 |
| `tet_p3_mg_mpiref` | fine | 0.859 | 0.202 | 4.26 |
| `tet_p3_mg_mpiref` | coarse | 0.277 | 0.095 | 2.92 |
| `tet_p3_mg_mpi2` | fine | 0.810 | 0.195 | 4.16 |
| `tet_p3_mg_mpi2` | coarse | 0.312 | 0.090 | 3.48 |

Repetitions: p3: 2/2/2/2/2/2/2/2/2/2, q2pm1: 2/2/2/2/2/2/2/2/2/2

### 4 threads, run mode (5 time steps; NK and MFNK)

| case | cfg | step time (excl. 1st) [s] | run time [s] | Jacobian builds total [s] | builds | NL residual total [s] | lits | nlits | per-step its = p3 | final B, V, Ti rel. diff vs p3 |
|---|---|---|---|---|---|---|---|---|---|---|
| `tet_p2_run` | p3 | 1.869 | 12.987 | 1.079 | 7 | 0.248 | 56 | 18 | True | 8.4e-12, 3.8e-12, 1.8e-12 |
| `tet_p2_run` | q0pm1 | 1.718 | 12.406 | 0.455 | 7 | 0.331 | 56 | 18 | True | 2.3e-12, 1.3e-12, 5.9e-13 |
| `tet_p2_run` | q2pm1 | 1.708 | 11.951 | 0.405 | 7 | 0.220 | 56 | 18 | True | 9.1e-12, 3.4e-12, 1.8e-12 |
| `tet_p3_mg_run` | p3 | 14.712 | 93.785 | 6.925 | 6 | 0.754 | 217 | 16 | True | 4.8e-10, 3.5e-11, 1.2e-11 |
| `tet_p3_mg_run` | q0pm1 | 13.184 | 86.167 | 1.984 | 6 | 0.961 | 217 | 16 | True | 9.2e-10, 8.1e-11, 2.2e-11 |
| `tet_p3_mg_run` | q2pm1 | 13.458 | 87.766 | 1.938 | 6 | 0.563 | 217 | 16 | True | 4.5e-10, 3.3e-11, 1.2e-11 |
| `tet_p3_mg_mfnk_run` | p3 | 13.447 | 86.270 | 5.685 | 5 | 5.766 | 164 | 10 | True | 9.5e-12, 2.7e-12, 8.5e-13 |
| `tet_p3_mg_mfnk_run` | q0pm1 | 13.013 | 83.822 | 2.028 | 5 | 8.135 | 164 | 10 | True | 1.1e-09, 7.2e-11, 2.3e-11 |
| `tet_p3_mg_mfnk_run` | q2pm1 | 12.027 | 79.550 | 1.701 | 5 | 4.751 | 164 | 10 | True | 4.2e-11, 2.7e-12, 8.4e-14 |

The last column compares the median of each configuration with the first `p3` repetition, so the `p3` rows show the
spread between `p3`'s own two repetitions (OpenMP reduction order). `q0pm1` has one repetition.

- **Iteration counts are identical at every step** for NK (tet p2, tet p3 p-MG) and MFNK (tet p3 p-MG). This holds
  for the original Jacobian (`p3`), the new direct path (`q0pm1`) and the cache path (`q2pm1`).
- Final field norms differ by ≤1.1e-9, the same as same-mode repeats at 4 threads, which differ by up to 1.1e-9
  (Phase 3; OpenMP reduction/atomic order).
- Time per run in Jacobian builds drops 2.7-3.6x. The step time drops only 9-11% (tet p3 p-MG NK: 14.7 → 13.5 s
  per step): the builds were already a minor part of these runs (8% → 2.5% of the step at tet p3 p-MG). The GMRES
  solve with block-Jacobi ILU / p-MG smoothing, including the lazily computed ILU factorizations, dominates (see
  next steps). `XMHD_PROF preupdate` does not capture the factorizations because they happen inside the solve.

## Memory

### Memory (1 thread; identical at 4 threads except RSS noise, see results/p4ab_t4_*)

| case | cells | qcache [MB] | coarse-level caches [MB] | posmap (all levels) [MB] | posmap KB/cell (finest) | posmap setup [s] | RSS p3 [MB] | RSS q2pm0 [MB] | RSS q2pm1 [MB] | matvec_norm rel q2pm1 vs p3 |
|---|---|---|---|---|---|---|---|---|---|---|
| `tet_p2` | 12288 | 4.9 | 0.0 | 85 | 6.9 | 0.78 | 792 | 792 | 874 | 0.0e+00 |
| `tet_p3` | 12288 | 7.6 | 0.0 | 312 | 25.4 | 2.20 | 3489 | 3489 | 3787 | 0.0e+00 |
| `tet_p3_mg` | 12288 | 7.6 | 2.6 | 398 | 25.4 | 2.93 | 4195 | 4198 | 4577 | 0.0e+00 |
| `tet_p4` | 1536 | 3.8 | 0.0 | 115 | 75.1 | 0.77 | 1557 | 1558 | 1668 | 0.0e+00 |
| `hex_p3` | 128 | 44.1 | 0.0 | 35 | 274.9 | 0.21 | 723 | 784 | 817 | 7.8e-16 |
| `curved_p3_mg` | 16305 | 26.2 | 3.3 | 528 | 25.4 | 4.05 | 5285 | 5288 | 5791 | 0.0e+00 |

- The **position map** needs 4 bytes per local entry per (row space, column space) pair. That is
  `4*sum(nce_s*nce_t)` bytes per cell: 6.9 KB/cell for tets p2, 25 KB/cell (p3) and 75 KB/cell (p4), and 275 KB/cell
  for hexes p3. Summed over p-MG levels this is 9-19% of the Jacobian's CSR storage (12 bytes per non-zero;
  tet p3: 312 MB vs 2.2 GB). RSS grows by ~8-13% (tet p3: 3.49 → 3.79 GB). The memory is independent of the thread count (per rank).
  `xmhd_jac_posmap_maxmb` skips levels above a limit, and `xmhd_jac_posmap=0` disables the map.
- **Caches**: the finest-level cache (residual and Jacobian) is unchanged from Phase 3 (1.3-44 MB). Coarse p-level
  caches add 0.5-3.3 MB, since they share the finest geometry.
- Per-thread work arrays (`E`, `W`, `C`) are a few MB per thread at p=3 (eg. Lagrange `W` is ~1000 x 105 doubles).

## Experimental coarse-level quadrature (`xmhd_coarse_quad=T`, "cq") vs default, 4 threads

| case | coarse build default [s] | coarse build cq [s] | ratio | fine build default [s] | fine build cq [s] | step time default [s] | step time cq [s] | lits default | lits cq | nlits default | nlits cq |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `tet_p3_mg` | 0.730 | 0.391 | 1.87 | 1.576 | 1.576 | - | - | - | - | - | - |
| `tet_p4_mg` | 0.490 | 0.258 | 1.90 | 0.757 | 0.728 | - | - | - | - | - | - |
| `curved_p3_mg` | 1.008 | 0.595 | 1.69 | 2.103 | 2.267 | - | - | - | - | - | - |
| `tet_p3_mg_run` | 0.102 | 0.064 | 1.59 | - | - | 13.025 | 13.593 | 217 | 217 | 16 | 16 |

`tet_p3_mg_run` q2pm1 per-step (lits, nlits): [(42, 3), (48, 4), (38, 3), (42, 3), (47, 3)]

`tet_p3_mg_run` q2pm1 per-step (lits, nlits): [(42, 3), (48, 4), (38, 3), (42, 3), (47, 3)]

`tet_p3_mg_run` cq per-step (lits, nlits): [(42, 3), (48, 4), (38, 3), (42, 3), (47, 3)]

`tet_p3_mg_run` cq per-step (lits, nlits): [(42, 3), (48, 4), (38, 3), (42, 3), (47, 3)]

- Coarse-level builds are 1.5-2.3x faster with the lower quadrature order (fine levels unchanged, within noise).
- On `tet_p3_mg_run` the linear and nonlinear iteration counts are identical at every step with and without the
  switch: the p-MG preconditioner with ILU smoothing is insensitive to the under-integration of the coarse
  operators here. This is one configuration. MG convergence should be checked on harder problems (strong
  advection, large `dt`, anisotropic conduction) before the switch is considered for a default. It stays
  experimental and off.

## Remaining hotspots

`perf` of profile mode with the final build (1 thread, `--prof-reps 3,0,0`, `runs/p4perf/*/perf_report.txt`). Shares
of the samples attributed to the Jacobian build (the rest of these runs is setup: projections, graphs, the first
build):

| part | tet p3 (5.1 s/build) | hex p3 (3.0 s/build) | notes |
|---|---|---|---|
| scatter with the position map (`xmhd_jac_add_native`) | 27% | 5% | ~34k entries per cell at tet p3, random access into the 27 MB value array (memory bound) |
| trial-function contraction (`xmhd_jac_wform`) | 26% | 19% | vectorized 4-6 term updates over `nce` per used row of `D` |
| `dgemm` | 17% | 57% | OpenBLAS small-matrix kernels on tets; large-matrix kernels + packing on hexes (64/144/125 Lagrange/H(Curl)/Grad(H^1) DOFs per cell) |
| pointwise kernel (`xmhd_jac_qp`) | 8% | 2% | ≤43x43 `D` per qp, zeroed per qp |
| scaling/BCs (`xmhd_jac_finish_cell`) | 6% | <1% | |
| field reconstruction (`fem_qcache_interp_cell`) | 6% | 6% | |
| basis rebuild (`fem_qcache_eval_cell`) | 5% | 11% | |
| `Jac%zero` | 4% | 1% | |

Tuning in this phase: `W` is now written once per element instead of read-modify-written once per non-zero of
`D` (2.3 → ~1 s per build at tet p3), and the scatter uses row-major positions, row by row (1.2 → ~1.0 s). The
build went from 6.5 to 5.1 s. Further CPU options:
- compact `D` to its static non-zero pattern, so `wform` loops over non-zeros only;
- skip zero rows (Dirichlet BC rows) in the scatter;
- reorder cells for locality of the matrix rows they touch;
- `-march`-level vectorization of the Fortran loops. The generic `-O2` build is SSE2 only, while OpenBLAS
  dispatches AVX-512 at runtime.

## Recommendations

### Defaults

- Keep `xmhd_qcache=2` (residual and Jacobian on the cache) and `xmhd_jac_posmap=1`. Mode 0 stays the reference
  implementation (used by the equivalence checks), and mode 1 (FULL) remains for the residual only.
- Document the BLAS requirement for threaded runs: an OpenMP-aware BLAS (OFT's `build_libs.py` default) or
  `OPENBLAS_NUM_THREADS=1`. The `oft-build` skill's system-library recipe installs the pthreads OpenBLAS. It
  should select `libopenblas0-openmp` (or set the variable) so that sandbox timings match OFT's own builds.
- Keep `xmhd_coarse_quad` off. It is an opt-in experiment (results above).

### Time stepping (now dominated by the linear solve and preconditioner)

With the Jacobian 3-8x cheaper, the NK/MFNK step time is dominated by GMRES iterations with block-Jacobi ILU / p-MG
smoothing and by preconditioner factorization updates. Suggested next steps:
1. **Reuse the symbolic factorization** of the ILU/LU blocks between preconditioner updates. The sparsity pattern
   is fixed, so only the numeric phase needs repeating.
2. **Retune the update frequencies** (`xmhd_prefreq`, `nl_update`). Jacobian rebuilds are now cheap relative to
   factorizations, so rebuilding J more often while refactoring less often (or vice versa) may reduce iteration
   counts at little cost.
3. **MFNK**: each GMRES iteration costs one residual evaluation (Phase 3: 1.6x faster). A cheaper assembled
   Jacobian also makes "assembled J for the matvec, matrix-free only for the nonlinear residual" competitive,
   and that should be re-measured.
4. **p-MG**: coarse levels are now cheap to rebuild. Deeper p hierarchies (p=4→3→2→1), or more smoothing on coarse
   levels, should be re-evaluated for iteration counts. The experimental coarse quadrature mostly changes coarse
   build time (see above).
5. **Matvec/smoother kernels**: native CSR matvec and ILU triangular solves are bandwidth bound. A block CSR
   (BSR) layout over the field blocks would cut index traffic. It would also make the position map smaller.

### GPU porting

The Jacobian cell loop now has the structure of a GPU kernel:
- per cell: gathers (DOF tables), `fem_qcache_interp_cell` and `fem_qcache_eval_cell` (cell-independent reference
  tables, which fit in shared memory up to p=3, plus per-qp 3x3 geometry), a pointwise kernel with no allocation
  or type-bound calls, two batched GEMMs per row space, and a scatter with precomputed positions. On a GPU the
  scatter is atomic adds with no search, which is the main reason to precompute positions.
- **Batched GEMM**: the shapes are fixed per configuration (eg. tet p3: Lagrange `20 x 35 x ~600` (value
  components) and `20 x 105 x ~240` (gradients); H(Curl) `45 x 105 x ~400`). Batch over cells with
  `cublasDgemmStridedBatched`, or one cell per thread block with `W` in shared memory.
- **`D` is sparse with a static pattern**: ≤43x43 per qp, of which roughly a quarter is non-zero for the
  default two-temperature model (estimate from the block structure). For GPU, generate the pattern
  once from the plan and the settings (or record it at the first qp), and store only the non-zeros per qp. This cuts
  register/shared-memory pressure and lets `wform` loop over non-zeros only. The same compaction would help on CPUs.
- Constants and settings: put the module settings read by `xmhd_jac_qp`/`xmhd_resid_qp` into one constant struct.
  Replace the `res_profile` procedure pointer by a tabulated profile, and precompute shaping factors per (cell, qp)
  (as recommended in Phase 3).
- Boundary conditions (`xmhd_jac_finish_cell`: Dirichlet rows, `vbc` rotations, revolved-periodic packing) are
  per-cell and branchy, but only touch boundary/periodic DOFs. On a GPU, apply them as a separate pass over flagged
  rows of the assembled matrix, or over flagged cells only.
- Solid cells use a separate block plan. Bin cells by kind (and by orientation signature for H(Curl) tets) to avoid
  divergence.

### Not done

- PETSc matrices keep `mat_add_local` (no position map), as specified.
- No canonical-frame (`fem_qcache_mat_from_ref`) variant: the basis rebuild is only 5-11% of the build after the
  batched contraction, and the canonical frame would add per-block transforms of every local matrix (and does not
  apply to H(Curl) faces on tets).
- `xmhd_massmatrix_apply` (linear runs) is still on direct evaluation, as in Phase 3.

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
  with random vectors agree to ≤1.8e-15. The CSR position map vs searching: ≤5.4e-16 (assembly order only).
- Run mode: linear and nonlinear iteration counts are identical at every step for NK (tet p2, tet p3 p-MG)
  and MFNK (tet p3 p-MG), for the original, new direct and new cache Jacobians (see below).
- Speed: see the A/B section. Jacobian build on tets p3, 1 thread: 22.5 s → 5.1 s.
- New options: `xmhd_jac_posmap` (0/1, default 1), `xmhd_jac_posmap_maxmb` (memory cap per level and rank) and
  the experimental, default-off `xmhd_coarse_quad`. Memory: see the memory section.

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
cache-blocked hand-written loop by 2x (see the A/B section), so `dgemm` is the default (`xmhd_jac_blas`, a
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
   canonical frame does not apply. Measured: basis evaluation is now <2% of the build, so the canonical tables plus
   `fem_qcache_mat_from_ref` would gain nothing;
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
- Built on the first Jacobian build of a level (one search pass, 0.2-2 s, reported in `XMHD_PROF jac`), and
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
the Jacobian of every level, built at the check state, as CSR (`xmhd_jac_dump`). `jac_compare.py` runs the 33
configurations of `test_MUG_qcache.py`, plus a 3-level p-MG case (`pmg_p4`), with any executable, and compares
dumps entrywise. The reference is a build of the unmodified Phase 3 branch plus the dump hook only, so its
Jacobian is the original code. The runs used 1 thread because the projected initial state is thread
nondeterministic: CG reductions make separate 2-thread runs differ by up to 1e-6 in row-relative Jacobian entries,
so cross-run comparisons need 1 thread. In-process comparisons are not affected. Largest difference relative to the row maximum (`results/phase4_jdump_*.txt`):

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

EQUIV_TABLE

## Validation

- **Regression** (`results/phase4_regression.log`, new defaults: residual and Jacobian on the cache, position
  map), from `builds/build_release/tests` with `OMP_NUM_THREADS=2`: `test_alfven.py`, `test_sound.py` with
  `-k "test_nl_r1_p2 or test_lin_r1_p2 or test_hex_nl_r1_p2 or (test_sound and test_nl_r1_p3) or visc or gyro or
  hall or qcache"` (incl. 2-rank NK/MFNK and linear runs, sound p3 MG), `test_MUG_visc.py`, `test_MUG_gyro.py`,
  `test_MUG_hall.py`, `test_MUG_qcache.py`: **77 passed, 20 skipped** (PETSc and excluded parameter combinations),
  0 failed. The known p-MG failures (alfven `test_nl_r1_p3`/`test_lin_r1_p3`) are outside this subset.
- **Runtime checks**: Debug build (`-fcheck=all,no-array-temps`) of `xmhd_bench`, cache-path check runs of
  `pmg_p4` (3 levels), `floors_solid`, `revolved[tet]`, `vbc[norm]`, `hex[2]`: no runtime errors, results as
  in release.
- **Valgrind** memcheck (release build, 1 thread; p-MG p3→p2, Hall, hyper terms, gyroviscosity, shaping, solid
  cells, floor sources; residual and Jacobian checks): 0 errors (`results/phase4_valgrind.txt`). No leak record
  points into the new code (the largest, 0.74 MB, is the driver's mass-matrix graph).

PERF_SECTIONS

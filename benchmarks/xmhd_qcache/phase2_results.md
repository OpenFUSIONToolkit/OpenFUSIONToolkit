# Phase 2: quadrature-point cache module (xMHD quadrature-caching project)

Branch `claude/xmhd-qcache-phase2`, starting from `claude/xmhd-qcache-phase1` (`b10a39c`). Phase 2 adds a
reusable module in `src/fem` with precomputed basis and geometry data at quadrature points. It has two
storage modes behind one interface, plus regression tests and a kernel benchmark. The xMHD physics loops are
unchanged; that port is Phase 3. The only change to existing code is that `xmhd_setup_cell_shared` now calls the
moved, identical flagging routine `fem_cell_shared_flags`.

## Summary

- **New module `fem_quad_cache`** (`src/fem/fem_quad_cache.F90`):
  - `oft_qcache_geom`: per-mesh geometry (`gop`, `cgop`, `det` once per straight/affine cell, per qp otherwise),
    per-cell edge/face orientation states, optional physical qp coordinates and thread-partition flags;
  - `oft_fem_qcache`: Lagrange/H^1/H(Curl) basis data in **FULL** or **REFERENCE** mode;
  - allocation-free kernels per (cell, qp) and per cell (all qp at once).
- **Validation:** both modes match the direct evaluators to at most 3.4e-14 relative, and to ~1e-15 in most
  checks. Covered: every (cell, qp) on straight/periodic/curved tets p=1-4, affine and curved hexes p=1-3,
  2 MPI ranks, the p-1 levels at the finest quadrature, field reconstruction, test-function accumulation and local
  matrices. All 20 new tests pass. Existing `src/tests/fem` tests and the MUG subset pass (see Regression tests).
- **Orientation structure (key check):** every DOF's function depends only on its own entity's orientation, for
  every basis, order and cell type tested. Each orientation state is either an exact **signed permutation** of the
  canonical (state 1) functions or not reachable by any transform within the same space:
  - hexes: signed permutations for all three spaces;
  - tets: signed permutations for Lagrange (all p), H^1 order <= 3 and all edges. H^1 faces of order >= 4 mix with
    the face's edges/vertices and the interior bubbles. H(Curl) faces (p >= 2) are **not in the span of the
    cell's H(Curl) functions at all**: they differ by gradients of H^1 functions, which OFT keeps in the separate
    grad(H^1) space.

  The REFERENCE tables are therefore indexed by (entity, orientation state). States without a transform use their own
  table block, and states with one can be applied as a per-cell coefficient permutation.
- **Kernel A/B (per pass over all cells and qp; ratios to the current direct evaluation, 1 / 4 threads):**
  - REFERENCE mode with the **cell-batched** kernels (`refcel`) is fastest wherever fields are reconstructed:
    reconstruction takes **0.31-0.38x** the direct time on tets p3/p4/curved (0.51/0.62x at p2) and 0.19/0.31x
    on hex p3. Reconstruction + accumulation (residual-like) takes 0.63-0.67x on tets p3/p4/curved
    (0.88/1.07x at p2) and 0.46/0.57x on hex p3.
  - FULL mode is memory-bandwidth bound on tets: p3 streams 100 KB/cell (1.2 GB per pass), and reconstruction is
    1.0-1.2x the direct time. On hex p3 it is the fastest for basis-only lookups (0.16/0.20x) and gives 0.38/0.47x
    for reconstruction, at 4.7 MB/cell.
  - Basis-only rebuild in REFERENCE mode costs about the same as the hand-optimized direct tet evaluators
    (0.8-1.3x). It is 0.38x on hexes.
- **Memory:** FULL needs 21 / 100 / 305 / 4747 KB per cell (tet p2/p3/p4, hex p3), as Phase 0 estimated.
  REFERENCE needs 0.1-1.4 MB of tables on tets and 22 MB on hex p3, plus 200 B per straight cell of geometry
  (~5 KB per curved cell).
- **Recommendation for Phase 3:** use REFERENCE mode with the cell-batched kernels as the default. Restructure the
  residual into gather -> `fem_qcache_interp_cell` -> pointwise physics over all qp -> `fem_qcache_accum_cell` ->
  scatter. Keep FULL as an option for hexes, where memory allows.

## Design

### Objects and lifecycle

```fortran
USE fem_quad_cache
TYPE(oft_qcache_geom), TARGET :: geom
TYPE(oft_fem_qcache) :: qlag, qh1, qhc
CALL geom%setup(mesh, quad, store_pts=shaping_active, store_shared=.TRUE.)  ! once per mesh + quadrature
CALL qlag%setup(lag_fe, geom, qcache_mode_ref)    ! values + gradients (default for Lagrange)
CALL qh1%setup(h1_fe, geom, qcache_mode_ref)      ! gradients (default for H^1)
CALL qhc%setup(hcurl_fe, geom, qcache_mode_ref)   ! values + curls (default for H(Curl))
...
WRITE(*,*) geom%mem(), qlag%mem(), qlag%setup_time  ! memory [bytes] and setup time
CALL qlag%delete; CALL geom%delete                   ! explicit teardown (idempotent)
```

- `setup` calls `delete` first, so re-setup is safe. `delete` is safe on objects that were never set up. `ready()`
  reports whether a cache is usable. Setup aborts with a clear message for unsupported FE types, a mesh mismatch,
  or a geometry cache that is not set up. Setup must be called outside OpenMP parallel regions; the kernels are
  thread-safe.
- The quadrature rule is an argument, so a coarse p-level FE space can be cached at the **finest** quadrature, as
  xMHD's coarse builds use. This is tested (`*_pm1` checks).
- No module or threadprivate state. Template evaluation resets the standard evaluators' threadprivate
  `oriented_cell` before and after use.

### Layout: `(dof, comp, qp, cell)` and `(row, lcomp, qp)`

The DOF index is fastest in every table. Field reconstruction (`sum_j c_j phi_j`) and the test-function loop
(`res_j += phi_j . flux`) are then both unit-stride over DOFs, which are the long loops (10-125). The component
index (1 or 3) is short and is unrolled into separate accumulators. All data for one (cell, qp) is one contiguous
`nce x ncomp` block (one `nce x ncomp*nqp` block per cell), so the same layout is a GEMM operand over quadrature
points, either for local-matrix assembly or for the cell-batched kernels. `(dof, qp, comp, cell)` would split a
qp's components across `nce*nqp` strides, which is worse for the per-qp kernels and no better for GEMM. The
prototype's `(dof, qp, cell, comp)` layout is avoided.

### Geometry (`oft_qcache_geom`)

- `goff(cell)` and `gstride(cell)` (0/1) index the geometry: entry `g = goff(cell) + (qp-1)*gstride(cell)`.
  - Straight cells (`cell_is_curved` false) store one entry, as xMHD does.
  - Curved tets and **all hexes** are tested at setup. If `gop` and `det` at every qp agree with qp 1 to 1e-14
    relative, one entry is stored, otherwise one per qp. The periodic cube hexes are detected as affine, the sphere
    hex27 and curved-cylinder cells are not (3180 of 16305 cells per qp on `cyl_tilt`).
- **Folded logical frame:** OFT's tet Jacobian sets `gop(:,1) = -(gop(:,2)+gop(:,3)+gop(:,4))` (barycentric
  coordinates), so only 3 columns are independent. The cache stores `gop(:,2:4)` and the cross products
  `cgop(:,4:6)` for tets, and `gop(:,1:3)` and `cgop(:,(/1,2,4/))` for hexes. Gradients, H(Curl) values and curls
  then all have exactly 3 logical components. This is 25% (gradients) and 50% (curls) fewer flops than the 4/6
  components the evaluators use, and the logical representation becomes unique, which the transform analysis
  below relies on. `fem_qcache_geom_get` returns `gop(3,4)`/`cgop(3,6)` and `det*w` in the evaluators' format.
- `orient(ne+nf, cell)` holds the edge states (1: `lce>0`, 2: `lce<0`) and the face states
  (`lcfo>0 -> lcfo`, `lcfo<0 -> face_np-lcfo`, i.e. 1..6 for triangles and 1..8 for quads).
- Optional `pt(3,nqp,nc)` (physical qp coordinates, eg. for xMHD shaping functions; enable when
  `xmhd_shaping_active()`) and `cell_shared(nc)` (the Phase 1 partition-boundary flags, now computed by
  `fem_cell_shared_flags` in `src/fem`; `xmhd_setup_cell_shared` calls it).

### FULL mode

`q(iq)%full(nce, npc, nqp, nc)` holds physical values/gradients/curls, built with the standard evaluators and the
cached geometry. Lookups are a copy, or direct use of the `full(:,:,qp,cell)` slice.

### REFERENCE mode

- **Tables:** `q(iq)%ref(nrows, 3, nqp)` (1 component for scalar values) has one block of rows per
  (local entity, orientation state). Rows `1:nce` hold state 1 of every entity in local DOF order (the
  "canonical table"); the other states follow. Vertices and the interior have one state, edges 2, faces 6 (tri)
  or 8 (quad).
- **Generation:** the tables come from the **standard evaluators** (with unit `gop`/`cgop`, as the evaluators are
  linear in them) on an internal 6-8 cell orientation-template mesh. It is allocated with `MOLD=` of the real mesh
  and shares its `cell_ed`/`cell_fc`. Template cell `s` has every edge in state `min(s,2)` and every face in
  state `s`. The FE object is shallow-copied with the mesh pointer swapped, so nothing is modified and no
  basis evaluator needs a new option. Hex is handled by the same code; the prototype's `cell_orient` flag silently
  did not handle it. Tables are therefore exact copies of what the evaluators compute, for any basis/order/cell
  type they support.
- **Orientation transforms:** setup fits `phi^s_j = sum_r M_jr phi^1_r` for every (entity, state) (LAPACK
  `dgels`). It tries the regressors in order: the entity's own canonical functions; plus those of its
  vertices, edges and the interior; finally all canonical functions of the cell. Fits within 1e-10 of a signed
  permutation of the own block are snapped to exact `xf_perm`/`xf_sgn`. Remaining good fits are kept as dense
  blocks, and failed fits are marked kind 3 (use the state's own block). The accepted fits have residuals of
  1e-15 or less.
- **Using the tables:** there are three equivalent paths, all exact.
  1. *Gather:* `fem_qcache_cell_rows(cache, cell, rows)` once per cell, then the per-qp kernels with `rows`.
     Works for every space.
  2. *Canonical frame, per qp:* `fem_qcache_coefs_to_ref` (coefficients, once per cell), the per-qp kernels
     without `rows` (contiguous), then `fem_qcache_res_from_ref` / `fem_qcache_mat_from_ref` (results, once per
     cell). Only when `fem_qcache_canonical_ok(cache)`, i.e. no kind 3 states.
  3. *Cell-batched* (`fem_qcache_interp_cell`, `fem_qcache_accum_cell`): coefficients are mapped to the canonical
     frame once (signed permutations), and the canonical table is contracted for **all qp at once** (a
     `nce x 3nqp` by `nce x nfld` product). Kind 3 entities get a small correction from their own state block
     (`nd x 3nqp`, `lda = nrows`, no copy). The per-qp geometry transform is applied last, and accumulation is the
     transpose. There are no per-DOF gathers, it works for every space, and it is the fastest path.

### Kernel API (all allocation-free, no `SELECT TYPE`/type-bound calls, thread-safe)

| Routine | Purpose |
|---|---|
| `fem_qcache_geom_get(geom,cell,qp,gop,detw[,cgop])` | `gop(3,4)`, `det*w`, `cgop(3,6)` (evaluator format) |
| `fem_qcache_cell_rows(q,cell,rows)` | table row of each local DOF (gather path) |
| `fem_qcache_eval(q,iq,cell,qp,out[,rows])` | physical basis `out(nce,npc)` at one qp |
| `fem_qcache_eval_cell(q,iq,cell,out[,rows])` | physical basis `out(nce,npc,nqp)` for a cell (GEMM operand) |
| `fem_qcache_interp(q,iq,cell,qp,coefs,nfld,out[,rows])` | `out(npc,nfld) = sum_j coefs(j,f) phi_j` at one qp |
| `fem_qcache_accum(q,iq,cell,qp,flux,res[,rows])` | `res(j) += phi_j . flux` at one qp |
| `fem_qcache_interp_cell(q,iq,cell,coefs,nfld,out,work)` | `out(npc,nqp,nfld)` for all qp (cell frame in/out) |
| `fem_qcache_accum_cell(q,iq,cell,flux,nfld,res,work)` | `res(nce,nfld) += sum_qp phi . flux(npc,nqp,nfld)` |
| `fem_qcache_coefs_to_ref`, `fem_qcache_res_from_ref`, `fem_qcache_mat_from_ref` | canonical-frame maps (no-op in FULL mode) |
| `fem_qcache_work_size(q,nfld)`, `fem_qcache_canonical_ok(q)`, `q%ready()`, `q%mem()`, `geom%mem()` | queries |

`iq` is `qcache_val`, `qcache_grad` or `qcache_curl`. For H(Curl), `qcache_val` is the vector value. Fluxes
passed to the accumulation kernels must include `det*w`.

## Validation

`src/tests/fem/test_fem_qcache.{F90,py}`. For each case, Lagrange (order p; values+gradients), H^1 (order p+1;
values+gradients) and H(Curl) (order p; values+curls) are set up on the finest level, and also their p-1 levels at
the finest quadrature. Three paths are checked: FULL, REFERENCE-gather (`ref`) and REFERENCE-canonical (`refc`,
falling back to gather where no transform exists). Each is compared with the direct evaluators **at every
(cell, qp)**, using the Jacobian evaluated at every qp as the reference geometry. The checks are:
- basis values/gradients/curls (`eval`);
- reconstruction of 2 fields from deterministic coefficients (`interp`, `icell`);
- test-function accumulation over all qp (`accum`, `acell`);
- a local matrix `sum_qp phi_i . phi_j` on the first 64 cells, including `fem_qcache_mat_from_ref`;
- `det*w` and physical qp coordinates.

Errors are max-normalized (by the max basis magnitude, or by `sum |c phi|` for reconstructions); the tolerance is
1e-12.

| test | mesh | max error (all checks) |
|---|---|---|
| `test_tet[1-4]` | cube, straight tets (1536 cells) | <= 1e-15 |
| `test_tet_periodic[2,3]` | triply periodic cube | <= 1e-15 |
| `test_tet_mpi[2,3]` (2 ranks) | distributed cube | <= 1e-15 |
| `test_hex[1-3]` | cube hexes (affine, collapsed to 1 geometry entry/cell) | 3e-15 |
| `test_hex_mpi` (2 ranks) | distributed cube hexes | <= 3e-15 |
| `test_curved_tet[1-4]`, `test_curved_tet_mpi` | `sphere_tet10_test` (grid_order=2, curved boundary cells) | 2.5e-15 |
| `test_curved_hex[1-3]` | `sphere_hex27_test` (grid_order=2, all cells per-qp geometry) | 1.4e-15 |

All 20 tests pass (`OMP_NUM_THREADS=2`, about 45 s). Every edge state (2) and every face state (6 tri / 8 quad)
occurs on these meshes (eg. tet cube face states: 1221/799/1020/1017/1256/831; curved hex: all 8), so each
table block is exercised against real cells. In the benchmark runs below, the reconstruction/residual
checksums of all six variants agree with direct evaluation to <= 3.4e-14 relative.

### Orientation structure (minimal representation)

Number of (entity, state) pairs with state > 1 by transform kind (identity / signed permutation / dense / none),
from the test runs:

| space | tet p=1 | tet p=2 | tet p=3 | tet p=4 | hex p=1-3 |
|---|---|---|---|---|---|
| Lagrange (order p) | - | 6/0/0/0 | 20/6/0/0 | 0/26/0/0 | all identity/perm |
| H^1 (order p+1) | 6/0/0/0 | 20/6/0/0 | 0/10/16/0 (faces: 16 regressors) | 0/10/16/0 (25 regressors) | all identity/perm |
| H(Curl) (order p) | 0/6/0/0 | 0/10/0/16 | 0/10/0/16 | 0/10/0/16 | all identity/perm |

- Edges and Lagrange faces are signed permutations (edge flips reverse the order/sign of the modes).
- H^1 tet faces of order >= 4 (3+ face modes) are linear combinations of the canonical face, edge, vertex
  **and interior** functions (`nd x 16` and `nd x 25` blocks). The interior bubbles vanish on every face, so the
  extensions of rotated face modes into the cell differ by bubbles.
- H(Curl) tet faces (p >= 2) are not in the span of any canonical H(Curl) functions of the cell. A diagnostic in
  the test fits them against the canonical H(Curl) functions **plus the gradients of the canonical H^1 (order
  p+1) functions**. That fit is exact (1e-15 for values, 5e-15 to 2e-14 for curls with the same H(Curl)
  coefficients, p=2-4). So OFT's gradient-free H(Curl) space is not rotation-invariant on its own, but
  H(Curl) + grad(H^1) is. This is harmless for conformity (neighbors share the global orientation). It means a
  per-cell transform for B fields would have to couple the two spaces. The cache uses the per-state table blocks
  instead (kind 3), which is exact.
- An earlier version fitted the redundant 4/6-component (unfolded) representation. That made the fit too strict
  (spurious failures for H^1 order >= 4 and H(Curl)+grad). The folded frame removes this.

## Kernel benchmark

Machine: Intel Xeon @ 2.30 GHz cloud VM (4 vCPUs, shared host), Ubuntu 24.04, GCC/gfortran 13 `-O2` with LTO,
OpenMPI 4.1, system OpenBLAS/SuperLU/ARPACK/HDF5/METIS (`oft-build` skill, section 5), `OFT_BUILD_TESTS=ON`.
Raw results are in `results/phase2_kernels_t4.json` and `results/phase2_kernels_t1.json`.

`src/tests/fem/fem_qcache_bench.F90` (built with the tests) and `benchmarks/xmhd_qcache/run_qcache_bench.py`. These
use the Phase 0 harness meshes and orders, with the xMHD spaces on the finest level and the finest H(Curl)
quadrature. Each timing is one pass over all cells and qp, wall clock, `!$omp do schedule(static)` over cells.
Variants and kernels are interleaved over 7 (4 threads) or 5 (1 thread) repetitions. Tables show the median
(min) in ms and median/direct.

Kernels:
- **basis**: physical Lagrange values+gradients, H^1 gradients, H(Curl) values+curls into local arrays;
- **recon**: basis + reconstruction of the xMHD fields from local coefficients, approximating
  `xmhd_interp_apply`: 6 Lagrange fields (V, N, Ti, Te) with gradients, B from grad(H^1) + H(Curl), J from
  curl(H(Curl)); this includes the per-cell coefficient gather;
- **resid**: recon + test-function accumulation of all fields (`res_j += phi_j f + grad phi_j . G`, etc.), the
  structure of the residual without the physics.

Variants:
- `direct`: the current code path (`mesh%jacobian` once per straight cell, per qp for curved cells/hexes; `oft_*_eval_all`);
- `full` / `ref` / `refc`: FULL, REFERENCE-gather, REFERENCE-canonical-where-available, with the per-qp kernels;
- `fullc` / `refcel`: FULL and REFERENCE with the cell-batched kernels (basis via `fem_qcache_eval_cell` with rows).

#### tet_p2 (12288 cells, 14 qp/cell, 4 thread(s); nce lag/h1/hcurl = 10/20/14)

Time per pass over all cells and quadrature points, median (min) of 7 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 9.16 (8.31) | 1.00 | 21.17 (19.90) | 1.00 | 32.31 (26.10) | 1.00 |
| full | 7.18 (6.72) | 0.78 | 23.74 (22.26) | 1.12 | 37.01 (35.29) | 1.15 |
| ref | 12.10 (11.94) | 1.32 | 14.16 (13.24) | 0.67 | 32.32 (31.40) | 1.00 |
| refc | 12.92 (12.08) | 1.41 | 22.23 (21.80) | 1.05 | 64.19 (40.06) | 1.99 |
| fullc | 11.14 (9.15) | 1.22 | 33.50 (24.83) | 1.58 | 54.60 (35.67) | 1.69 |
| refcel | 20.34 (11.05) | 2.22 | 13.02 (11.24) | 0.62 | 34.50 (25.24) | 1.07 |

Max relative checksum difference vs direct (recon, resid): 5.5e-15

#### tet_p3 (12288 cells, 35 qp/cell, 4 thread(s); nce lag/h1/hcurl = 20/35/29)

Time per pass over all cells and quadrature points, median (min) of 7 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 47.65 (39.86) | 1.00 | 93.34 (88.64) | 1.00 | 131.29 (117.15) | 1.00 |
| full | 37.22 (33.22) | 0.78 | 107.11 (102.03) | 1.15 | 153.94 (151.24) | 1.17 |
| ref | 43.68 (43.13) | 0.92 | 50.21 (49.42) | 0.54 | 132.97 (114.75) | 1.01 |
| refc | 49.84 (47.66) | 1.05 | 83.57 (81.51) | 0.90 | 154.20 (135.42) | 1.17 |
| fullc | 46.56 (43.87) | 0.98 | 111.06 (108.90) | 1.19 | 177.98 (166.04) | 1.36 |
| refcel | 50.85 (44.10) | 1.07 | 35.28 (33.18) | 0.38 | 85.44 (78.82) | 0.65 |

Max relative checksum difference vs direct (recon, resid): 3.0e-15

#### tet_p4 (1536 cells, 61 qp/cell, 4 thread(s); nce lag/h1/hcurl = 35/56/53)

Time per pass over all cells and quadrature points, median (min) of 7 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 18.22 (17.78) | 1.00 | 34.42 (33.72) | 1.00 | 47.47 (44.33) | 1.00 |
| full | 15.69 (14.99) | 0.86 | 41.18 (37.09) | 1.20 | 62.95 (57.10) | 1.33 |
| ref | 17.02 (15.87) | 0.93 | 18.39 (17.80) | 0.53 | 48.18 (44.21) | 1.02 |
| refc | 17.91 (17.60) | 0.98 | 26.75 (26.33) | 0.78 | 45.17 (44.38) | 0.95 |
| fullc | 19.89 (19.41) | 1.09 | 41.84 (40.31) | 1.22 | 74.33 (68.41) | 1.57 |
| refcel | 17.94 (16.86) | 0.99 | 11.28 (11.06) | 0.33 | 30.74 (30.30) | 0.65 |

Max relative checksum difference vs direct (recon, resid): 2.3e-15

#### hex_p3 (128 cells, 512 qp/cell, 4 thread(s); nce lag/h1/hcurl = 64/125/88)

Time per pass over all cells and quadrature points, median (min) of 7 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 152.31 (119.91) | 1.00 | 146.98 (141.50) | 1.00 | 188.83 (158.02) | 1.00 |
| full | 30.62 (21.24) | 0.20 | 68.36 (58.52) | 0.47 | 80.60 (78.66) | 0.43 |
| ref | 59.42 (51.68) | 0.39 | 47.36 (39.54) | 0.32 | 73.36 (66.41) | 0.39 |
| refc | 49.79 (43.08) | 0.33 | 61.52 (45.94) | 0.42 | 76.86 (65.11) | 0.41 |
| fullc | 27.66 (26.43) | 0.18 | 61.87 (56.22) | 0.42 | 104.73 (98.84) | 0.55 |
| refcel | 62.33 (59.14) | 0.41 | 45.09 (34.80) | 0.31 | 107.11 (95.59) | 0.57 |

Max relative checksum difference vs direct (recon, resid): 3.4e-14

#### curved_p3 (16305 cells, 35 qp/cell, 4 thread(s); nce lag/h1/hcurl = 20/35/29)

Time per pass over all cells and quadrature points, median (min) of 7 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 69.68 (68.51) | 1.00 | 130.91 (126.98) | 1.00 | 172.80 (170.79) | 1.00 |
| full | 50.43 (48.79) | 0.72 | 140.01 (133.90) | 1.07 | 208.09 (200.29) | 1.20 |
| ref | 60.70 (57.33) | 0.87 | 69.57 (68.75) | 0.53 | 174.82 (158.01) | 1.01 |
| refc | 64.22 (63.82) | 0.92 | 116.53 (110.34) | 0.89 | 186.92 (183.14) | 1.08 |
| fullc | 67.30 (61.82) | 0.97 | 163.53 (152.59) | 1.25 | 242.46 (226.15) | 1.40 |
| refcel | 61.95 (57.23) | 0.89 | 50.05 (49.57) | 0.38 | 114.53 (110.89) | 0.66 |

Max relative checksum difference vs direct (recon, resid): 1.3e-14

### Single thread

#### tet_p2 (12288 cells, 14 qp/cell, 1 thread(s); nce lag/h1/hcurl = 10/20/14)

Time per pass over all cells and quadrature points, median (min) of 5 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 33.18 (31.66) | 1.00 | 77.74 (73.92) | 1.00 | 97.56 (95.91) | 1.00 |
| full | 28.34 (25.88) | 0.85 | 89.33 (85.44) | 1.15 | 137.71 (134.25) | 1.41 |
| ref | 40.79 (39.73) | 1.23 | 51.25 (49.97) | 0.66 | 118.41 (116.34) | 1.21 |
| refc | 44.39 (41.94) | 1.34 | 80.89 (80.49) | 1.04 | 144.43 (138.83) | 1.48 |
| fullc | 36.01 (35.47) | 1.09 | 97.53 (93.38) | 1.25 | 144.33 (136.30) | 1.48 |
| refcel | 43.65 (41.15) | 1.32 | 39.67 (39.38) | 0.51 | 85.70 (84.15) | 0.88 |

Max relative checksum difference vs direct (recon, resid): 6.9e-15

#### tet_p3 (12288 cells, 35 qp/cell, 1 thread(s); nce lag/h1/hcurl = 20/35/29)

Time per pass over all cells and quadrature points, median (min) of 5 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 152.97 (147.81) | 1.00 | 333.41 (322.83) | 1.00 | 459.14 (441.48) | 1.00 |
| full | 121.09 (118.77) | 0.79 | 382.58 (371.71) | 1.15 | 586.94 (581.22) | 1.28 |
| ref | 168.22 (165.87) | 1.10 | 188.11 (184.81) | 0.56 | 438.35 (432.62) | 0.95 |
| refc | 183.68 (176.98) | 1.20 | 318.10 (311.66) | 0.95 | 534.90 (521.17) | 1.17 |
| fullc | 168.68 (167.32) | 1.10 | 417.56 (408.02) | 1.25 | 627.71 (615.20) | 1.37 |
| refcel | 173.07 (167.35) | 1.13 | 127.22 (124.83) | 0.38 | 302.67 (297.96) | 0.66 |

Max relative checksum difference vs direct (recon, resid): 2.4e-15

#### tet_p4 (1536 cells, 61 qp/cell, 1 thread(s); nce lag/h1/hcurl = 35/56/53)

Time per pass over all cells and quadrature points, median (min) of 5 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 69.92 (67.40) | 1.00 | 133.67 (130.13) | 1.00 | 172.05 (169.13) | 1.00 |
| full | 52.74 (50.81) | 0.75 | 137.55 (134.15) | 1.03 | 219.11 (211.25) | 1.27 |
| ref | 62.68 (61.09) | 0.90 | 72.02 (69.62) | 0.54 | 171.72 (169.49) | 1.00 |
| refc | 68.59 (66.43) | 0.98 | 103.27 (101.33) | 0.77 | 174.52 (170.81) | 1.01 |
| fullc | 66.58 (66.08) | 0.95 | 147.66 (144.31) | 1.10 | 244.91 (239.84) | 1.42 |
| refcel | 68.14 (66.03) | 0.97 | 42.07 (40.03) | 0.31 | 115.18 (111.76) | 0.67 |

Max relative checksum difference vs direct (recon, resid): 2.8e-15

#### hex_p3 (128 cells, 512 qp/cell, 1 thread(s); nce lag/h1/hcurl = 64/125/88)

Time per pass over all cells and quadrature points, median (min) of 5 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 458.51 (445.15) | 1.00 | 555.43 (540.15) | 1.00 | 603.39 (593.10) | 1.00 |
| full | 75.65 (73.74) | 0.16 | 212.22 (210.84) | 0.38 | 268.18 (262.44) | 0.44 |
| ref | 174.92 (170.82) | 0.38 | 133.47 (132.63) | 0.24 | 238.93 (237.20) | 0.40 |
| refc | 143.81 (140.56) | 0.31 | 166.74 (161.48) | 0.30 | 237.82 (235.98) | 0.39 |
| fullc | 84.77 (83.32) | 0.18 | 202.02 (197.41) | 0.36 | 321.58 (317.16) | 0.53 |
| refcel | 187.23 (184.75) | 0.41 | 107.31 (106.13) | 0.19 | 275.09 (269.94) | 0.46 |

Max relative checksum difference vs direct (recon, resid): 1.1e-14

#### curved_p3 (16305 cells, 35 qp/cell, 1 thread(s); nce lag/h1/hcurl = 20/35/29)

Time per pass over all cells and quadrature points, median (min) of 5 interleaved repetitions, ms. Ratio = median / direct median.

| variant | basis | ratio | recon | ratio | resid | ratio |
|---|---|---|---|---|---|---|
| direct | 277.68 (267.01) | 1.00 | 520.21 (496.09) | 1.00 | 690.13 (675.21) | 1.00 |
| full | 164.78 (161.11) | 0.59 | 526.52 (523.44) | 1.01 | 785.13 (770.77) | 1.14 |
| ref | 220.72 (217.85) | 0.79 | 257.90 (252.89) | 0.50 | 600.85 (592.89) | 0.87 |
| refc | 254.55 (250.35) | 0.92 | 431.96 (421.49) | 0.83 | 718.87 (702.26) | 1.04 |
| fullc | 226.77 (223.45) | 0.82 | 583.14 (579.47) | 1.12 | 869.98 (855.97) | 1.26 |
| refcel | 230.82 (224.58) | 0.83 | 188.54 (185.11) | 0.36 | 431.78 (426.13) | 0.63 |

Max relative checksum difference vs direct (recon, resid): 5.0e-15

### Setup time and memory

| case | cells | qp | geom setup [s] | geom [MB] | FULL setup lag/h1/hcurl [s] | FULL [MB] | FULL [KB/cell] | REF setup lag/h1/hcurl [s] | REF tables [MB] | max RSS [MB] |
|---|---|---|---|---|---|---|---|---|---|---|
| tet_p2 | 12288 | 14 | 0.0006 | 2.46 | 0.014/0.019/0.024 | 253.2 | 20.6 | 0.0002/0.0034/0.0006 | 0.09 | 312 |
| tet_p3 | 12288 | 35 | 0.0004 | 2.46 | 0.051/0.072/0.115 | 1235.2 | 100.5 | 0.0006/0.0038/0.0024 | 0.45 | 1346 |
| tet_p4 | 1536 | 61 | 0.0002 | 0.31 | 0.021/0.022/0.042 | 469.2 | 305.5 | 0.0041/0.0019/0.0128 | 1.40 | 557 |
| hex_p3 | 128 | 512 | 0.0075 | 0.05 | 0.058/0.077/0.094 | 607.7 | 4747.3 | 0.0315/0.0356/0.0589 | 22.38 | 665 |
| curved_p3 | 16305 | 35 | 0.0285 | 19.70 | 0.079/0.115/0.180 | 1639.0 | 100.5 | 0.0007/0.0035/0.0030 | 0.45 | 1764 |

- FULL setup costs about 1.5-5 direct basis passes (it is one direct evaluation, plus the stores). REFERENCE
  setup (template evaluation plus the transform fits) takes 1-15 ms for tets and ~0.13 s for hex p3; most of it
  is the least-squares fits. Geometry setup takes <1 ms for straight meshes and 30 ms for the curved cylinder.
- Geometry costs 200 B per straight/affine cell (`gop`+`cgop` 18 doubles, `det`, offsets, 10-18 orientation
  ints) and ~5.3 KB per curved cell (per-qp entries at 35 qp).
- REFERENCE memory is dominated by geometry for large tet meshes; the tables are 0.1-1.4 MB (tets) and 22 MB
  (hex p3, 512 qp).

### Interpretation

1. **Batching is what matters.** The per-qp REFERENCE gather kernels (`ref`) already reconstruct fields 1.5-2x
   faster than direct evaluation, because the physical basis is never formed (logical contraction, then one
   3x3 transform). They give nothing for basis-only or the residual on tets: OFT's tet evaluators are
   hand-unrolled and cheap. Contracting the cell's coefficients with the canonical table for all qp at once
   (`refcel`) also removes per-qp call overhead and per-DOF gathers. Reconstruction then runs 2.6-3.2x faster on
   tets p3/p4/curved and 3-5x faster on hex p3. The canonical per-qp variant (`refc`) gains nothing over the gather;
   the canonical frame only pays off in the batched form.
2. **FULL mode is bandwidth bound on tets.** Basis lookups are a copy at memory bandwidth. Reconstruction from
   FULL tables reads 100 KB/cell at p3 and is no faster than evaluating. It wins only where evaluation is
   expensive (hex p3: basis 0.16-0.20x direct) and costs 4.7 MB/cell there.
3. **The residual-like kernel** gains less (0.63-0.67x on tets p3/p4/curved) than reconstruction, because
   accumulation into 6+2 residual vectors is a large share of it in every variant. Phase 3's pointwise physics
   adds the same absolute cost to every variant. At p2 (10-20 DOFs, 14 qp) the per-cell work is too small for
   batching to pay off (0.88-1.07x), though reconstruction alone is still 0.5-0.6x.
5. **Noise.** Repeat 4-thread runs on this VM differ by up to ~30% (eg. `tet_p2` `refcel` basis 20 ms median,
   11 ms min). The 1-thread tables agree with the 4-thread medians to ~10% in ratio and are the more reliable.
4. **Hexes** gain the most in every mode (no per-qp `mesh%jacobian` for affine cells, and table lookups replace the
   tensor-product evaluation). Phase 0 found basis evaluation is 63% of the hex residual.

## Regression tests

From `builds/build_release/tests` with `OMP_NUM_THREADS=2`:

@@REGRESSION@@

## Notes for Phase 3

### Default mode and API

- **Default: REFERENCE mode with the cell-batched kernels** (`fem_qcache_interp_cell`/`fem_qcache_accum_cell`).
  Memory is negligible and setup is milliseconds. It is exact for every space, order and cell type, and it is the
  fastest variant for reconstruction and the residual. Offer FULL as an option (eg. `xmhd_qcache=2`) for hex
  meshes, where basis-only lookups are 5-6x faster than direct, if memory allows (4.7 MB/cell at p3).
- **Setup:** in `xmhd_setup_ops`/`xmhd_run`, once per mesh, with explicit teardown:
  `geom%setup(mesh, xmhd_ML_hcurl%levels(nlevels)%fe%quad, store_pts=xmhd_shaping_active(), store_shared=.TRUE.)`,
  then one `oft_fem_qcache` per space **and per FE level** (the coarse p-levels at the fine quadrature: 3 small
  caches per level). Do not store them only in `xmhd_run` (prototype problem 4). Tie them to the xMHD ops object,
  rebuild if the mesh changes, and free them in the ops teardown.
- **Residual (`xmhd_errmatrix_apply`)**, per cell:
  1. gather local coefficients (as `xmhd_interp_apply` does): `lf(nce_l,5[+Te,N2])`, `bgrad(nce_h)`,
     `bcurl(nce_c)[,J2]`;
  2. `fem_qcache_interp_cell` for: Lagrange `qcache_val` and `qcache_grad` (all Lagrange fields at once, `nfld` =
     5-7), H^1 `qcache_grad` (B part), H(Curl) `qcache_val` (B part) and `qcache_curl` (J). Outputs are
     `(3,nqp,nfld)` arrays; `fem_qcache_geom_get` (or `geom%det`, `geom%wts`) gives `det*w`, and `geom%pt` gives
     the shaping positions;
  3. pointwise physics over all qp of the cell, filling flux arrays `(npc,nqp,nfld)` including `det*w`
     (vectorizable over qp);
  4. `fem_qcache_accum_cell` for each (space, quantity), then the Phase 1 scatter (`geom%cell_shared` for
     atomics).

  Workspace: one `work(fem_qcache_work_size(q, nfld_max))` per thread.
- **Jacobian (`xmhd_build_ops`)**: get physical basis blocks per cell with `fem_qcache_eval_cell(q, iq, cell, B,
  rows)` (`B(nce,3,nqp)`, exact for all spaces), then form the local matrices
  `K_ij = sum_q sum_cd B_i,c(q) D_cd(q) B_j,d(q)` as batched GEMMs over qp: scale `B` by `D*detw` per qp, then
  one `(nce x 3nqp) * (3nqp x nce)` DGEMM per block. Where `fem_qcache_canonical_ok` holds for both spaces
  (Lagrange, H^1 order <= 3, everything on hexes), the GEMM can use the shared canonical table plus
  `fem_qcache_mat_from_ref`. The (much larger) scatter cost from Phase 0/1 is separate. Precomputed CSR positions
  remain the main lever there.
- **Coarse p-levels:** cache each coarse level's spaces with the **fine** `geom` (same quadrature, as
  `xmhd_build_ops` uses). The tests check this (`*_pm1`). REFERENCE tables for p2 at the p3 quadrature are about
  0.1 MB, so this is essentially free. It replaces all coarse-level `oft_*_eval_all` calls. Whether to lower the
  coarse quadrature order (Phase 0 conclusion 2) is an independent decision; it would need a second `geom`.

### Known bugs in the prototype residual (`xmhd_mg_tests`), to avoid when porting

1. **Velocity gradient transposed:** the prototype built `dV(j,:) = grad V_j`, but `xmhd_interp_unpack` uses
   `dV(:,k) = grad V_k`, i.e. `val(6:14)=RESHAPE(dlf0(:,1:3),(/9/))` with `dlf0(:,k)` = gradient of field `k`. With
   `fem_qcache_interp_cell(q, qcache_grad, ...)` on the 3 velocity fields, `out(:,m,k)` is grad V_k: pack it the
   same way.
2. **Missing scaling:** the prototype dropped `n2_scale` on N2/dN2 and `j2_scale` on J2/J2c. Apply them exactly
   as `xmhd_interp_unpack` does.

### Other notes

- The per-qp kernels (`fem_qcache_eval/interp/accum`) exist for incremental porting (drop-in for the current
  per-qp structure). They are about as fast as direct evaluation on tets, so the gain only comes with the
  per-cell restructuring.
- Don't use the canonical per-qp path for H(Curl) on tets (p >= 2): `fem_qcache_canonical_ok` is false there,
  and `fem_qcache_coefs_to_ref` aborts. The cell-batched kernels handle it automatically.
- GPU mapping: per cell, the REFERENCE kernels are one small dense contraction against a cell-independent table
  that fits in shared memory (0.1-1.4 MB for tets), plus a 3x3 transform per qp and signed-permutation index
  maps. Geometry is 19 doubles per cell (straight) and per-cell orientation codes. Replace the `cell_shared`
  atomics with coloring if needed.
- Prototype problems 1-4 from the brief are addressed: layout (DOF-contiguous); no `cell_orient` option, so hex
  is handled by the template mesh; no dense `cell_trans` (signed permutations, or small per-entity blocks
  applied once per cell); explicit setup/teardown/memory reporting and safe behavior when not set up.

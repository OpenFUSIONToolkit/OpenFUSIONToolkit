# Phase 3: xMHD nonlinear residual on the quadrature-point cache

Branch `claude/xmhd-qcache-phase3`, starting from `claude/xmhd-qcache-phase2` (`5328aab`). Phase 3 ports the
xMHD nonlinear residual (`xmhd_errmatrix_apply`) to the Phase 2 cache (`src/fem/fem_quad_cache.F90`). It checks
that the new path gives the same result as the old one to round-off, and measures the speedup.

**Summary**
- New option `xmhd_qcache` in `&xmhd_options`: 0 = direct basis evaluation (default, unchanged), 1 = cache in
  FULL mode, 2 = cache in REFERENCE mode. There is no fallback path: every configuration in the test matrix is
  supported by both cache modes.
- The pointwise physics is one routine, `xmhd_resid_qp`, shared by the direct and the cache path. The direct path
  still gives **bit-identical** checksums to the Phase 2 build (1 thread).
- Equivalence: 32 configurations (physics options, tets p2-p4, hexes p2/p3, curved tets, BCs, revolved
  periodicity, p-multigrid, 2 MPI ranks) × 3 residual evaluations (`dt>0`, `dt>0` with `back_step`, `dt<0`
  diagnostics) × 3 mode pairs. The largest difference is 3e-14 relative to the residual of the same
  equation, and diagnostics agree to 1e-15 (see below). NK and MFNK runs give identical iteration counts.
- Speed (REFERENCE vs direct, NL residual): see the tables. The new path is never slower in the measured
  cases.
- One Phase 2 kernel change: transposed reference tables, so that accumulation is a unit-stride dot product
  (accumulation 1.7x faster in the residual).

## Design

### Shared pointwise physics

`xmhd_resid_qp(u0, np_prev, v_prev, jp_prev, pt, shp, det, dt, he, eta_cell, solid, back_step, neg_src, nflag,
diag_vals, neg_vols, fl_v, fl_g, fh_v, fh_c, fg_g)` holds everything that used to be inside the quadrature loop of
`xmhd_errmatrix_apply` after the field reconstruction. That covers floors, transport coefficients (Braginskii,
shaping, `eta_temp`/`res_profile`), diagnostics, the solid-region branch, B/J2/Bg, V (JxB, gravity, advection,
pressure, `kin`/`iso` viscosity, `nu_par`, `nu_gyro`, upwinding), n/n2, T and Te. The code is the original,
statement for statement. It maps the unpacked fields (`xmhd_loc_values`, from `xmhd_interp_unpack`) to
test-function fluxes that already include `det*w`:

| residual block | test-function contraction |
|---|---|
| Lagrange field `k` (1-3 V, 4 n, 5 T, 6 Te, 7 n2) | `res_j += phi_j*fl_v(k) + grad(phi_j).fl_g(:,k)` |
| H(Curl) field `k` (1 B, 2 J2) | `res_j += phi_j.fh_v(:,k) + curl(phi_j).fh_c(:,k)` |
| Grad(H^1) part of B | `res_j += grad(phi_j).fg_g` |

`xmhd_resid_masks(dt, act_fluid, act_solid)` gives the active fluxes for fluid and solid cells. These match the
original `IF` structure (eg. J2 only for `dt>0`; solid cells only B). The previous-step reconstructions are passed
in unscaled, and the kernel applies `den_scale`/`vel_scale` exactly where the old code did. Shaping factors are
passed in as `shp(4)`, which is 1 for functions that are not set (multiplying by 1 is exact), so the kernel has
no procedure-pointer calls for shaping. It still calls `res_profile` if `eta_temp>0`.

**Trade-off.** Both loops call the kernel once per quadrature point, and the direct loop contracts the returned
fluxes with its basis arrays. Two contractions change form at round-off level: the H(Curl) mass term goes from
`(phi.B)*det` to `phi.(B*det)` (same for J2). This did not change the checksums of the direct path. With one
thread, `nlop_norm`/`matvec_norm` are identical to all 16 printed digits against a build of `5328aab` on
`tet_p2`, `hex_p3` and `curved_p3_mg`. The extra call costs nothing measurable. With 1 thread the refactored
direct path is 1.04-1.28x *faster* than before, probably because the shaping/diagnostic branches are hoisted.
That is within VM noise.

### Cache path (`xmhd_resid_qcache`)

Per cell, with per-thread work arrays allocated once per call:
1. **gather** the local coefficients from the interpolator's local arrays (`full_interp%lf_loc`, `Te_loc`,
   `N2_loc`, `bcurl_loc`, `bgrad_loc`, `J2_loc`, which already have the revolved-periodic velocity unpacking)
   and the previous-step fields, using precomputed `dofs_*(nce,nc)` tables;
2. **reconstruct** with `fem_qcache_interp_cell`: Lagrange values for all Lagrange fields at once (V, n, T,
   [Te], [n2], [n_prev], [V_prev]) and gradients for the current ones; Grad(H^1) gradients; H(Curl) values
   (B, [J2]) and curls (B, [J2], [B_prev]);
3. **pack** each qp into the 40-value array exactly as `xmhd_interp_apply` does (in particular
   `vals(6+3(k-1):8+3(k-1)) = grad(V_k)`, so `dV(:,k)=grad(V_k)`). Then `xmhd_interp_unpack`, which applies
   `vel_scale`, `den_scale`, `n2_scale` and `j2_scale`, then `xmhd_resid_qp`;
4. **accumulate** the active fluxes with `fem_qcache_accum_cell` (Lagrange values and gradients for up to 7
   fields at once; H(Curl) values and curls; Grad(H^1) gradients);
5. **scatter** with the Phase 1 `xmhd_scatter_add` (atomics only for `geom%cell_shared` cells), with the same
   scale factors as the direct path.

BCs, the revolved-periodic repack, `restore_local` and the diagnostics reduction are shared by both paths.
`det*w` comes from `geom%det*geom%wts`. Physical points come from `geom%pt` when shaping is active (stored only
then). For the `dt<0` diagnostics without shaping they are computed per cell with `mesh%log2phys`, once per step.

Prototype bugs from `xmhd_mg_tests` (`xmhd_errmatrix_capply`) and how they are avoided:
- (a) velocity gradient orientation: the same packing and unpack routine as the direct path is used, so there is
  no separate convention to get wrong. Advection, `iso` viscosity, upwinding and viscous heating all match to 1e-15;
- (b, c) `n2_scale`, `j2_scale`, `vel_scale`, `den_scale` are applied by the shared `xmhd_interp_unpack`
  and kernel; the hyper-resistivity/diffusivity cases (`hyper`, `hyper_centered`, `test_qcache_tet/hex`) match to
  round-off.

### Selection and lifecycle

- `xmhd_qcache` (0/1/2, default 0) is read with `&xmhd_options` and documented in the options table of
  `xmhd_read_settings`. Invalid values abort. No configuration needs a fallback.
- The cache (`xmhd_qcache_data`: one `oft_qcache_geom`, three `oft_fem_qcache`, DOF tables, element sizes) is
  stored per FE level in `xmhd_ops%qc`. It is built for the finest level in `xmhd_run` right after
  `xmhd_alloc_ops` (with `store_pts=xmhd_shaping_active()` and `store_shared=.TRUE.`).
- `xmhd_errmatrix_apply` checks `xmhd_qcache_valid(oft_xmhd_ops%qc)` on every call: the mode, FE space and mesh
  pointers, and that points are stored if shaping is now active. It rebuilds the cache if the check fails. So
  every caller of `oft_xmhd_ops%A%apply` is safe: `xmhd_run` (diagnostics, NK), MFNK (`oft_mf_matrix` applies `A`
  directly), profile mode, the check mode, and any external caller. `xmhd_lin_run` does not apply `A`.
- Teardown: `xmhd_qcache_free` (all levels) at the end of `xmhd_run` and `xmhd_lin_run`, before `oft_finalize`
  in profile mode, and in the check mode.
- Setup time and memory are printed under `pm`/`debug>=1` and in `xmhd_prof_report`
  (`XMHD_PROF qcache mode,setup_s,mem_mb`, parsed by `run_bench.py`).

### Phase 2 kernel change

The first profile of the cache path (tet p3, 1 thread) had `fem_qcache_accum_cell` at 55% of the residual time,
twice `fem_qcache_interp_cell`. The canonical contraction in the accumulation was written as rank-1 updates
(`rw(:,f) += tab(:,kq)*lf(kq,f)`: two loads and a store per FMA). The kind-3 entity corrections (H(Curl) tet faces
at p>=2) used inner loops of 3-6 DOFs. REFERENCE mode now also stores the transposed table
`q(iq)%reft(nlc*nqp,nrows)` (0.1-1.4 MB for tets). Accumulation is then a unit-stride dot product over the
quadrature points, and the kind-3 corrections become unit-stride updates/dot products of length `nlc*nqp`.
Effect on the residual (tet p3, 1 thread): 0.52 s → 0.39 s per apply; hex p3: 0.37 s → 0.22 s. All 20
`test_fem_qcache` tests pass.

## Equivalence

`xmhd_run(..., qcache_check=.TRUE.)` (`xmhd_qcache_check`, driven by `xmhd_bench` `mode='qcheck'` and
`src/tests/physics/MUG/test_MUG_qcache.py`) evaluates `A%apply(u)` on the same projected nonlinear state (guide
field + ABC perturbation, sheared flow, modulated n/Ti/Te) with `xmhd_qcache` = 0, 2 and 1, for:
1. `dt>0`, no previous step;
2. `dt>0` with `self%up` = a copy of `u` with a smooth ±5% perturbation of every DOF (`back_step`: previous-step n,
   V (upwinding) and J (Hall));
3. `dt<0` (the diagnostics call: `diag_vals` compared too).

For each mode pair it reports, per field block, the maximum absolute difference and the block's maximum, the
relative difference of all 7 diagnostic values, and the maximum difference of the floor flags (`neg_flag`).
Test-only hooks: `qcheck_solid` marks every fourth cell solid, `qcheck_negsrc` sets a floor source pattern. The
bench options `den_floor`/`temp_floor` (set within the field range, so floors activate) and `g_accel` set the
corresponding module variables, `shaping_mix` uses different shaping functions for diff/cond, and `res_prof` sets
a `res_profile`.

**Metric.** Each block's difference is taken relative to the largest entry of the same equation: B = {Bc, Bg},
V = {Vx, Vy, Vz}, others by themselves. The Grad(H^1) block of B on its own is a cancellation residue
(`div(B)~0`). Its entries are 1e-4 (tets) to 1e-9 (hexes) of the individual quadrature contributions, so its
self-relative difference (shown separately) is meaningless at the 1e-12 level. The test requires ≤1e-12 per
equation, ≤1e-13 relative to the largest residual entry, ≤1e-12 for diagnostics and floor flags. An
implementation error gives O(1) differences.

{{QCHECK_TABLE}}

All 32 cases pass (`test_MUG_qcache.py`, 32 passed in ~105 s with 2 threads). The full log is in
`results/phase3_qcheck.log`.

Coverage, as required:
- one/two temperature (`single_temp_hall`, `brag_single_temp`); Hall on/off; upwinding on/off (`no_upwind_iso`);
- `visc_type` kin/iso, `nu_par` on/off (`kin_gyro_nopar`), `nu_gyro` (`iso_gyro`, `kin_gyro_nopar`, tet/hex);
- Braginskii (one/two temperature), thermal equilibration, `xmhd_diss_centered`;
- n2/j2 hyper terms, with `dt>0` and `dt<0` (diagnostics skip them) and with centering;
- all four shaping functions, with different functions per coefficient, plus `eta_temp` with and without `res_profile`;
- floors and floor sources, gravity, solid cells (`floors_solid`, `curved_p3`, `mpi2`);
- `vbc='norm'`/`'tang'` (`xmhd_vbcdir`) with Neumann n/T on a non-periodic box;
- revolved periodicity (reflected 30° tet and hex wedges from `test_vlag_periodic.py`, with `vbc='norm'`);
- straight tets p2-p4, curved (quadratic) tets p3 (`cyl_Marklin.h5`), hexes p2/p3, a p-MG hierarchy (finest
  level checked), and 2 MPI ranks (distributed mesh, `nlevels=2`);
- `adv_den`/`adv_temp`/`advec`/`jcb` off (`no_advec_den`, `no_temp`).

Valgrind (memcheck) on a cache-path run (tet p3, Hall, hyper terms, shaping, solid cells) reports 0 errors.

### Run mode (NK and MFNK)

From the A/B runs (`run_bench.py` run cases, 5 time steps, 4 threads):

{{RUN_TABLE}}

### Regression tests with the cache path

Temporary build with the default flipped to `xmhd_qcache=2` (one-line source edit, not committed; reverted and
rebuilt afterwards), from `builds/build_release/tests` with `OMP_NUM_THREADS=2`:

| suite | result |
|---|---|
| `test_alfven.py`, `test_sound.py` with `-k "test_nl_r1_p2 or test_lin_r1_p2 or test_hex_nl_r1_p2 or (test_sound and test_nl_r1_p3)"` (incl. 2-rank NK/MFNK variants and sound p3 MG), `test_MUG_visc.py`, `test_MUG_gyro.py`, `test_MUG_hall.py` | 43 passed, 20 skipped (8 PETSc, 12 excluded parameter combinations), 0 failed |

A breakpoint on `xmhd_qcache_setup` confirmed that the test binaries built the cache from `xmhd_run`. With
the default (0) build: `test_MUG_qcache.py` 32 passed and `fem/test_fem_qcache.py` 20 passed. The known p-MG
failures (alfven `test_nl_r1_p3`/`test_lin_r1_p3`) are not in this subset.

## Performance A/B

Machine: Intel Xeon cloud VM, 4 vCPUs (shared host, ±15-30% run-to-run noise). Ubuntu 24.04 system libraries
(`oft-build` section 5): GCC 13 `-O2` + LTO, OpenMPI 4.1, OpenBLAS, SuperLU, METIS. `run_bench.py` default
matrix (same as Phases 0-2), ILU block-Jacobi. Script `phase3_ab.sh`, tables from `phase3_tables.py`:
- 1 thread: profile cases with `--prof-reps 1,10,0` (one Jacobian build, 10 NL applies, no matvecs);
  `xmhd_qcache=0,2` interleaved, 2 repetitions each, plus 1 run of `xmhd_qcache=1`;
- 4 threads: full default matrix (profile, MPI and run/MFNK cases), `xmhd_qcache=0,2` interleaved, 2
  repetitions each, plus 1 run of `xmhd_qcache=1`.

Medians over repetitions. "NL apply" is profile mode's average time per `A%apply`. "loop/apply" is the Phase 1
cell-loop timer divided by the number of applies (it includes the diagnostics and checksum calls). The cache is
built once in `xmhd_run` and is not part of either.

{{AB_TABLES}}

## Remaining hotspots (tet p3, 1 thread, `perf`)

`perf record` of `xmhd_bench` profile mode with `--prof-reps 1,60,0` (62 NL applies), self time attributed per
function and divided by the number of applies:

| | direct (`xmhd_qcache=0`), 0.68 s/apply | REFERENCE (`xmhd_qcache=2`), 0.39 s/apply |
|---|---|---|
| field reconstruction | `xmhd_interp_apply` (+ `RESHAPE`) ~0.25 s (37%) | `fem_qcache_interp_cell` ~0.14 s (36%) |
| basis evaluation | `oft_*_eval_all` ~0.21 s (31%) | none (tables) |
| test-function accumulation + scatter | cell loop body ~0.14 s (20%) | `fem_qcache_accum_cell` ~0.16 s (41%); gather/pack/scatter ~0.02 s (5%) |
| pointwise physics | `xmhd_resid_qp` + `cross_product` ~0.04 s (6%) | `xmhd_resid_qp` + `cross_product` + unpack ~0.06 s (15%) |

In the cache path, ~75% of the time is two small dense contractions per (cell, quantity): coefficients ×
reference table, and table × fluxes. Shapes are `nce x (3*nqp)` by `nfld`, eg. 20 x 105 x 6 for Lagrange p3,
and 45 x 105 x 2 for H(Curl) p3 plus the kind-3 face blocks. Counting flops, this is about 60k FMAs per cell
at p3, which runs at roughly 2 GFLOP/s per core with the `-O2` (SSE2) build. On CPUs the next gains would come
from `-march`-level vectorization (AVX2/AVX-512 4-8x wider than the SSE2 the generic build targets), batching
several cells into one GEMM, or skipping zero fluxes. On GPUs these are exactly the batched-GEMM shapes that map
well (one cell per thread block, table in shared memory).

## Recommendations

### Default

**Flip the default to `xmhd_qcache=2` (REFERENCE).** {{FLIP_EVIDENCE}} Keep 0 for debugging and as a reference
implementation, since the equivalence check compares against it. Don't make FULL (1) the default: it stores
`nce*3*nqp` doubles per cell and space (≈100 KB/cell for tets p3, ≈4.7 MB/cell for hexes p3, see the memory
column), and it is not faster than REFERENCE in these cases. Keep it as an option for small hex meshes.

If the default is flipped:
- the MUG regression subset already passes with 2 (above);
- users with custom shaping functions get physical points stored (24 B/qp/cell), as needed;
- for the cases where the residual dominates the run time (MFNK, many NL iterations), the gain carries over
  directly: each matrix-free GMRES iteration is one residual apply.

### Next phase: Jacobian build (`xmhd_build_ops`)

The Jacobian build is now the largest single cost in every profile (20-27% of the residual-focused profile runs
above, and most of the time in run mode). Recommended structure, reusing this phase's pieces:
1. **Same per-cell skeleton as `xmhd_resid_qcache`**: gather, `fem_qcache_interp_cell` for the fields (the
   build reconstructs the same fields through `xmhd_interp_apply` today), and one shared pointwise kernel
   (`xmhd_jac_qp`, the analogue of `xmhd_resid_qp`). It returns the per-qp coefficient tensors `D_ab(q)` for
   every active block pair (value/gradient/curl of the test and trial functions), with `det*w` included.
2. **Batched local matrices via GEMM over qp**: get the physical basis per cell with
   `fem_qcache_eval_cell(q, iq, cell, B, rows)` (`B(nce,3,nqp)`, exact for all spaces). Form
   `K_ab = B_a (D_ab B_b)^T` as one `(nce_a x 3nqp) x (3nqp x nce_b)` GEMM per block pair, after scaling `B_b` by
   `D_ab` per qp. Where `fem_qcache_canonical_ok` holds for both spaces (Lagrange; H^1 order<=3; all hex spaces),
   use the canonical table plus `fem_qcache_mat_from_ref` and skip the per-cell basis build. The residual
   profile shows these small contractions run at ~2 GFLOP/s/core in the `-O2` SSE2 build, so a batched or
   `-march`-tuned GEMM matters.
3. **Precomputed CSR positions**: for each cell and block pair store the positions of the `nce_a x nce_b` local
   entries in the native matrix's value array (int32), computed once per sparsity pattern. The scatter then
   needs no `search4`/`mat_add_values` (7-10% of the profile even in residual-dominated runs, and the dominant
   build cost in Phase 0). Memory is `4*sum(nce_a*nce_b)` bytes per cell (≈60 KB/cell for tets p3, all blocks), so
   offer it for the fine level only, or per block.
4. **Coarse p-levels**: cache each coarse level's spaces on the fine `geom`, as Phase 2 tested (`*_pm1`).
   Decide separately whether to lower the coarse quadrature order (Phase 0 conclusion 2), which needs a second
   `geom` and is a numerics change (MG convergence must be rechecked).

### GPU porting

The cache path is close to a GPU kernel already:
- per cell: indexed gathers, two small dense contractions per space against **cell-independent** tables
  (Lagrange p3: 17 KB, H(Curl) p3 values+curls ~75 KB plus kind-3 face blocks; these fit in shared memory up to
  p=3), per-qp 3x3 geometry transforms (19 doubles per straight cell), the pointwise kernel, and atomics for
  shared cells (or a coloring);
- `xmhd_resid_qp` has no allocation, no type-bound calls and only reads module settings. Before porting,
  (i) put the settings into one constant struct, (ii) replace the `res_profile` procedure pointer with a
  tabulated or built-in profile, (iii) precompute the shaping factors per (cell, qp) at setup (4 doubles/qp,
  replacing `geom%pt` and the user-function calls), and (iv) provide device versions of `brag_*`,
  `elec_ion_therm_rate`, `cross_product`, `magnitude`;
- the orientation handling (`qcache_to_canon`/`from_canon`: signed permutations, small dense blocks, kind-3
  table blocks) is per cell and branchy. On a GPU, sort or bin cells by orientation signature, or apply the
  permutations as gathers inside the kernel;
- the `dt<0` diagnostics reduction and floor flags are per-cell reductions (`neg_flag` max, `diag_vals`
  sums).

### Not done in this phase

- `xmhd_massmatrix_apply` was not ported: `xmhd_lin_run` applies it once or twice per time step, and linear
  runs are dominated by the assembled-Jacobian solve, so the gain would be negligible. The same pattern (a
  shared pointwise kernel plus `xmhd_resid_qcache`'s skeleton) applies if it is needed.
- `xmhd_build_ops` (the Jacobian) is unchanged, as planned.


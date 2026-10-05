Toroidal and parallel current densities in TokaMaker and IMAS     {#doc_tokamaker_current_conventions}
================

TokaMaker's `jphi-*` profiles and IMAS (FUSE, IMAS.jl, FRESCO) use different, equally valid
definitions of "toroidal current density". This note derives the exact conversions between them
and from the flux-surface-averaged parallel current ⟨J·B⟩, as used by the bootstrap solver
(`solve_bootstrap` / `solve_with_bootstrap`).

## Setting
Axisymmetric equilibrium, **B** = F(ψ)∇φ + ∇ψ×∇φ, with ψ the poloidal flux per radian. The sign of
ψ is chosen so that the Grad–Shafranov equation reads Δ*ψ = −μ0R²p′ − FF′ (′ = d/dψ), so that

    (A1)  j_φ = R p′ + FF′/(μ0 R)

which is TokaMaker's convention (`util.get_jphi_from_GS`). COCOS 11 (IMAS) uses ψ per full turn
and a different sign: that changes only the factors 2π and ±1 in what follows, never the geometry.
Flux-surface average: ⟨A⟩ = ∮A dl/B_p ⁄ ∮dl/B_p, so V′ = dV/dψ = 2π∮dl/B_p. F, F′ and p′ are
flux functions.

## Parallel current
Ampère's law gives **J** = (F′/μ0)∇ψ×∇φ + j_φ R∇φ, so J·B = (F/R) j_φ + (F′/μ0) B_p². With (A1):

    (A2)  J·B = F p′ + F′B²/μ0,   hence   ⟨J·B⟩ = F p′ + F′⟨B²⟩/μ0

## The two toroidal conventions
Averaging (A1):

    (A3)  TokaMaker:  J_TM   ≡ ⟨j_φ⟩           = ⟨R⟩p′ + ⟨1/R⟩FF′/μ0
    (A4)  IMAS:       J_IMAS ≡ ⟨j_φ/R⟩/⟨1/R⟩  = (p′ + ⟨1/R²⟩FF′/μ0)/⟨1/R⟩

TokaMaker's jphi → FF′ map (`jphi_update`, `jphi_bs_update`) inverts (A3) exactly. IMAS.jl
`J_tor` and FRESCO `PressureJt` use (A4). On the FF′ part the two differ by ⟨1/R²⟩/⟨1/R⟩² = 1 + O(ε²),
largest at the edge.

## IMAS ↔ TokaMaker
Eliminating FF′/μ0 between (A3) and (A4):

    (A5)  J_TM   = ⟨R⟩p′ + ⟨1/R⟩ (J_IMAS⟨1/R⟩ − p′)/⟨1/R²⟩
          J_IMAS = [p′ + ⟨1/R²⟩ (J_TM − ⟨R⟩p′)/⟨1/R⟩] / ⟨1/R⟩

## Parallel → toroidal
From (A2), FF′/μ0 = F(⟨J·B⟩ − Fp′)/⟨B²⟩. Substituting into (A4):

    (A6)  J_IMAS = F⟨1/R²⟩⟨J·B⟩/(⟨B²⟩⟨1/R⟩) + p′ (1 − F²⟨1/R²⟩/⟨B²⟩)/⟨1/R⟩

which is IMAS.jl `JparB_2_JtoR` (with `includes_bootstrap=true`; COCOS 11 adds −2π to p′), a check on
the derivation. Substituting into (A3) instead:

    (A7)  J_TM = F⟨1/R⟩⟨J·B⟩/⟨B²⟩ + p′G,     G = ⟨R⟩ − F²⟨1/R⟩/⟨B²⟩

## Components
(A6) and (A7) are exact for the total current with the equilibrium's own F and averages, and at
fixed geometry they are linear in ⟨J·B⟩. So for a split ⟨J·B⟩ = Σ_k⟨J·B⟩_k (ohmic, bootstrap,
current drive, …) each component converts with the first term alone, and the single pressure term
is assigned to one component; following IMAS.jl (`includes_bootstrap=true`), to the bootstrap.
Building J_TM this way guarantees

    (A8)  the equilibrium's own ⟨J·B⟩ (= F p′ + F′⟨B²⟩/μ0) = Σ_k ⟨J·B⟩_k

The pressure term is the pressure-driven, non-field-aligned current (diamagnetic plus
Pfirsch–Schlüter, the latter with ⟨j_PS·B⟩ = 0). Neoclassical bootstrap and current-drive models
give the field-aligned part as ⟨J·B⟩ (Hirshman & Sigmar 1981; Sauter et al. 1999; Redl et al. 2021),
which is why it converts with the first term only. For B_p → 0, G → ⟨R⟩ − ⟨1/R⟩/⟨1/R²⟩ ≥ 0.

Earlier bootstrap conversions used ⟨J·B⟩⟨R⟩/F, then ⟨J·B⟩/⟨|B|⟩. For B_p → 0 these exceed the
first term of (A7) by ⟨R⟩⟨1/R²⟩/⟨1/R⟩ and ⟨1/R²⟩/⟨1/R⟩² respectively, and both omit p′G.

## Plasma current
The cross-section element between neighbouring surfaces is dA = dl dψ/|∇ψ| = dl dψ/(R B_p), so for
any flux function f

    (A9a) ∫ f dA = ∫dψ f ∮dl/(R B_p) = ∫dψ (V′/2π) f ⟨1/R⟩

and for the toroidal current

    (A9b) I_p = ∫ j_φ dA = ∫dψ ∮(j_φ/R) dl/B_p = ∫dψ (V′/2π) ⟨j_φ/R⟩

By (A4), ⟨j_φ/R⟩ = ⟨1/R⟩J_IMAS, so comparing with (A9a)

    (A9c) I_p = ∫ J_IMAS dA      (exact; J_IMAS from J_TM by (A5))
    (A9d) I_p = ∫dψ (V′/2π) [J_TM⟨1/R²⟩/⟨1/R⟩ + p′(1 − ⟨R⟩⟨1/R²⟩/⟨1/R⟩)]

TokaMaker evaluates (A9c) as `gs_flux_int(J_IMAS(J_TM))`, which is (A9a) on the plasma region.
Two earlier errors in this measure have been removed:
- the integrand was J_TM/(⟨R⟩⟨1/R⟩), i.e. ∫dψ (V′/2π) J_TM/⟨R⟩, which equals (A9b) only if
  ⟨j_φ/R⟩ = ⟨j_φ⟩/⟨R⟩;
- `gs_flux_int` integrated over the whole plasma region of the mesh, crediting every point outside
  the LCFS with the profile's LCFS value (the profile interpolator returns the LCFS end there). For a
  profile with a finite edge value this overstated I_p by several percent (ITER bootstrap case:
  5–7 %); ∫ 1 dA was the region's area, not the plasma's.

`jphi_bs_update` used to absorb both errors by rescaling its I_p target by
`gs_itor_nl`/`gs_flux_int` on the total. With the exact measure that factor converged to 1.0004
(1-D quadrature error) and has been removed, so the SWB inductive scale α is set by (A9c) alone;
converged I_p sits within a few 1e-4 of the target. `jphi_update` normalises its first iterate by
(A9c) and then follows the FEM current (`gs_itor_nl`).

## Where these are used
- TokaMaker `jphi_bs_update` / `calculate_bootstrap` (Fortran) and `solve_with_bootstrap(use_python_solve=True)`:
  the bootstrap enters `jphi_total` by (A7) with p′G; `boot_profs['jdotb_bs_raw']` holds the Redl ⟨J·B⟩.
- TokaMaker `jphi_update` / `jphi_bs_update` I_p normalisation, and the α root of
  `solve_with_bootstrap(use_python_solve=True)` (`compute_flux_integral` of J_IMAS): (A9c).
- bouquet: currents are stored as J_TM; FUSE/IMAS currents are read and written through (A5)–(A7).

## References
- J. P. Freidberg, *Ideal MHD* (Cambridge University Press, 2014): Grad–Shafranov equilibrium; (A1), (A2).
- J. Wesson, *Tokamaks*, 4th ed. (Oxford University Press, 2011): equilibrium and flux-surface averages.
- S. P. Hirshman and D. J. Sigmar, Nucl. Fusion **21**, 1079 (1981): neoclassical currents as ⟨J·B⟩;
  Pfirsch–Schlüter current with zero ⟨J·B⟩.
- O. Sauter, C. Angioni and Y. R. Lin-Liu, Phys. Plasmas **6**, 2834 (1999); A. Redl et al.,
  Phys. Plasmas **28**, 022502 (2021): bootstrap current as ⟨j_∥B⟩.
- O. Sauter and S. Yu. Medvedev, Comput. Phys. Commun. **184**, 293 (2013): COCOS (IMAS signs and 2π).
- IMAS.jl `src/physics/currents.jl`, `JparB_2_JtoR`:
  `<J⋅B> = (1/f)(<B²>/<1/R²>)(<Jt/R> + dp/dψ(1 − f²<1/R²>/<B²>))`, i.e. (A6).

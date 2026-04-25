module Bootstrap

# Port of bootstrap.py. Pure-Julia ports of:
#   - Hmode_profiles
#   - parameterize_edge_jBS
#   - calculate_ln_lambda
#   - redl_bootstrap (Redl et al. PoP 28, 022502 (2021))
#
# `analyze_bootstrap_edge_spike`, `solve_jphi`, `find_optimal_scale`, and
# `solve_with_bootstrap` are stubbed for incremental implementation — they
# couple to the Grad-Shafranov solver and to scipy's curve_fit / find_peaks
# which need careful rebuilding on Julia primitives. Each stub clearly
# states what's missing.

using Distributions: SkewNormal, pdf
using Roots: find_zero, Brent

export Hmode_profiles, parameterize_edge_jBS, calculate_ln_lambda,
       redl_bootstrap, get_jphi_from_GS

const _MU0 = 4π * 1e-7
const _EC = 1.602176634e-19

# ----------------------------------------------------------------------------

"""
    Hmode_profiles(; edge=0.08, ped=0.4, core=2.5, rgrid=201, expin=1.5, expout=1.5,
                     widthp=0.04, xphalf=nothing) -> Vector{Float64}

H-mode profile (tanh pedestal + power-law core), copied from OMFIT.
Returns the profile sampled on `range(0, 1; length=rgrid)`.
"""
function Hmode_profiles(; edge::Real=0.08, ped::Real=0.4, core::Real=2.5,
                          rgrid::Integer=201, expin::Real=1.5, expout::Real=1.5,
                          widthp::Real=0.04,
                          xphalf::Union{Nothing,Real}=nothing)
    w_E1 = 0.5 * widthp
    xph = xphalf === nothing ? 1.0 - w_E1 : Float64(xphalf)
    xped = xph - w_E1
    pconst = 1.0 - tanh((1.0 - xph) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)
    coretanh = 0.5 * a_t * (1.0 - tanh(-xph / w_E1) - pconst) + edge
    xpsi = collect(range(0.0, 1.0; length=rgrid))
    val = @. 0.5 * a_t * (1.0 - tanh((xpsi - xph) / w_E1) - pconst) + edge
    @inbounds for i in eachindex(xpsi)
        xtoped_i = xpsi[i] / xped
        if xtoped_i^expin < 1.0
            val[i] += (core - coretanh) * (1.0 - xtoped_i^expin)^expout
        end
    end
    return val
end

# ----------------------------------------------------------------------------

function _generate_baseline_prof(psi::AbstractVector, amp::Real, center::Real,
                                  width::Real, offset::Real, sk::Real,
                                  y_sep::Real, blend_width::Real, tail_alpha::Real)
    raw_shape(x) = pdf(SkewNormal(center, width, sk), x)
    x_grid = collect(range(max(0.0, center - 3 * width),
                            min(1.0, center + 3 * width); length=10000))
    y_grid = raw_shape.(x_grid)
    idx_peak = argmax(y_grid)
    val_peak_raw = y_grid[idx_peak]
    x_peak = x_grid[idx_peak]

    k_smooth = (amp / width) * (blend_width / 4.0)
    k_smooth < 1e-5 && (k_smooth = 1e-5)

    n = length(psi)
    profile_left = similar(psi, Float64)
    stitch_height = 0.0
    if amp > offset
        diff = amp - offset
        if diff > 1e-10
            argument = 1.0 - exp((offset - amp) / k_smooth)
            argument <= 0 && (argument = 1e-16)
            internal_amp = amp + k_smooth * log(argument)
        else
            internal_amp = amp
        end
        @inbounds for i in 1:n
            spike_i = raw_shape(psi[i]) / val_peak_raw * internal_amp
            # logaddexp(a, b) = max(a,b) + log(1 + exp(-|a-b|))
            a = offset / k_smooth; b = spike_i / k_smooth
            m = max(a, b); profile_left[i] = k_smooth * (m + log1p(exp(-abs(a - b))))
        end
        stitch_height = amp
    else
        @inbounds for i in 1:n
            profile_left[i] = offset
        end
        stitch_height = offset
    end

    dist_to_edge = max(1.0 - x_peak, 1e-5)

    if y_sep < stitch_height
        target_cos = (y_sep / stitch_height)^(1.0 / tail_alpha)
        target_cos = clamp(target_cos, -1.0, 1.0)
        ω = acos(target_cos) / dist_to_edge
        profile_right = similar(psi, Float64)
        @inbounds for i in 1:n
            base_cos = max(cos(ω * (psi[i] - x_peak)), 0.0)
            profile_right[i] = stitch_height * base_cos^tail_alpha
        end
    else
        target_cosh = (y_sep / stitch_height)^(1.0 / tail_alpha)
        ω = acosh(target_cosh) / dist_to_edge
        profile_right = @. stitch_height * cosh(ω * (psi - x_peak))^tail_alpha
    end

    out = similar(psi, Float64)
    @inbounds for i in 1:n
        out[i] = psi[i] <= x_peak ? profile_left[i] : profile_right[i]
    end
    if y_sep >= 0
        @inbounds for i in 1:n
            out[i] = max(out[i], 0.0)
        end
    end
    return out
end

"""
    parameterize_edge_jBS(psi, amp, center, width, offset, sk;
                          y_sep=0.0, blend_width=0.03, tail_alpha=1.5)

Parameterized edge bootstrap profile (skew-normal spike blended with a
powered-cosine tail). Inverts the offset via Brent's method to match the
user-specified `offset` value.
"""
function parameterize_edge_jBS(psi::AbstractVector, amp::Real, center::Real,
                                width::Real, offset::Real, sk::Real;
                                y_sep::Real=0.0, blend_width::Real=0.03,
                                tail_alpha::Real=1.5)
    if offset == 0.0
        return _generate_baseline_prof(psi, amp, center, width, -1e10, sk,
                                        y_sep, blend_width, tail_alpha)
    end
    function obj(α)
        prof = _generate_baseline_prof(psi, amp, center, width, α * offset, sk,
                                        y_sep, blend_width, tail_alpha)
        return prof[1] - offset
    end
    α_opt = find_zero(obj, (-1e10, 10 * amp), Brent(); rtol=1e-6)
    return _generate_baseline_prof(psi, amp, center, width, α_opt * offset, sk,
                                    y_sep, blend_width, tail_alpha)
end

# ----------------------------------------------------------------------------

"""
    calculate_ln_lambda(Te, Ti, ne, ni; Zeff=1.0, electron_lnLambda_model="sauter",
                        ion_lnLambda_model="unity") -> (ln_λ_e, ln_λ_ii)

Coulomb logarithms. `Te`, `Ti` in eV; `ne`, `ni` in m^-3. Both arrays/scalars.
"""
function calculate_ln_lambda(Te::AbstractArray, Ti::AbstractArray,
                              ne::AbstractArray, ni::AbstractArray;
                              Zeff::Union{Real,AbstractArray}=1.0,
                              electron_lnLambda_model::AbstractString="sauter",
                              ion_lnLambda_model::AbstractString="unity")
    ln_λ_e = if electron_lnLambda_model == "NRL"
        @. 23.5 - log(sqrt(ne / 1e6) * Te^(-5.0/4.0)) -
            sqrt(1e-5 + (log(Te) - 2)^2 / 16.0)
    else
        @. 31.3 - log(sqrt(ne) / Te)
    end
    Z_lnLam = if ion_lnLambda_model == "Zavg"
        @. max(ne / ni, 1.0)
    else
        1.0
    end
    ln_λ_ii = @. 30.0 - log(Z_lnLam^3 * sqrt(ni) / Ti^1.5)
    return max.(ln_λ_e, 10.0), max.(ln_λ_ii, 10.0)
end

# ----------------------------------------------------------------------------

function _F31_poly(X, Zeff)
    dZ = @. Zeff^1.2 - 0.71
    return @. (1 + 0.15 / dZ) * X - (0.22 / dZ) * X^2 +
              (0.01 / dZ) * X^3 + (0.06 / dZ) * X^4
end

"""
    redl_bootstrap(; psi_N, Te, Ti, ne, ni, pe, pi, Zeff, R, q, eps, fT, I_psi,
                     dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi=nothing,
                     dp_dpsi=nothing, ln_lambda_e=17.0, ln_lambda_ii=17.0,
                     use_legacy_L34=false, use_sign_q=false,
                     ion_collisionality_model="Zeff", formula_form="jB",
                     Zeff_override=nothing, nu_e_star_override=nothing,
                     nu_i_star_override=nothing) -> (j_BS, coeffs::Dict)

Bootstrap current via the Redl et al., Phys. Plasmas 28, 022502 (2021)
formulas. Returns the bootstrap current density (A/m^2) and a dictionary of
coefficients `L31`, `L32`, `L34`, `alpha`, etc.
"""
function redl_bootstrap(; psi_N::AbstractArray, Te::AbstractArray, Ti::AbstractArray,
                          ne::AbstractArray, ni::AbstractArray,
                          pe::AbstractArray, pi::AbstractArray,
                          Zeff::AbstractArray, R::AbstractArray, q::AbstractArray,
                          eps::AbstractArray, fT::AbstractArray, I_psi::AbstractArray,
                          dT_e_dpsi::AbstractArray, dT_i_dpsi::AbstractArray,
                          dn_e_dpsi::AbstractArray,
                          dn_i_dpsi::Union{Nothing,AbstractArray}=nothing,
                          dp_dpsi::Union{Nothing,AbstractArray}=nothing,
                          ln_lambda_e::Union{Real,AbstractArray}=17.0,
                          ln_lambda_ii::Union{Real,AbstractArray}=17.0,
                          use_legacy_L34::Bool=false, use_sign_q::Bool=false,
                          ion_collisionality_model::AbstractString="Zeff",
                          formula_form::AbstractString="jB",
                          Zeff_override::Union{Nothing,Real,AbstractArray}=nothing,
                          nu_e_star_override::Union{Nothing,AbstractArray}=nothing,
                          nu_i_star_override::Union{Nothing,AbstractArray}=nothing)
    Zeff_eff = Zeff_override === nothing ? Zeff :
        (Zeff_override isa AbstractArray ? Zeff_override : fill(Float64(Zeff_override), length(Zeff)))
    formula_form == "jboot1" && dn_i_dpsi === nothing &&
        error("dn_i_dpsi required for jboot1 form")

    p_total = pe .+ pi
    R_pe = pe ./ p_total

    nu_e_star = if nu_e_star_override !== nothing
        nu_e_star_override
    else
        @. 6.921e-18 * q * R * ne * Zeff_eff * ln_lambda_e / (Te^2 * eps^1.5)
    end

    nu_i_star = if nu_i_star_override !== nothing
        nu_i_star_override
    elseif ion_collisionality_model == "Koh"
        @. 4.90e-18 * q * R * ni * Zeff_eff * ln_lambda_ii / (Ti^2 * eps^1.5)
    elseif ion_collisionality_model == "unity"
        @. 4.90e-18 * q * R * ni * 1.0 * ln_lambda_ii / (Ti^2 * eps^1.5)
    else
        @. 4.90e-18 * q * R * ni * (Zeff_eff^4) * ln_lambda_ii / (Ti^2 * eps^1.5)
    end

    # --- L31 ---
    Zclip = max.(Zeff_eff .- 1.0, 0.0)
    ft_31_d1 = @. (0.67 * (1 - 0.7 * fT) * sqrt(nu_e_star)) / (0.56 + 0.44 * Zeff_eff)
    ft_31_d2 = @. ((0.52 + 0.086 * sqrt(nu_e_star)) * (1 + 0.87 * fT) * nu_e_star) /
                  (1 + 1.13 * sqrt(Zclip))
    X31 = @. fT / (1 + ft_31_d1 + ft_31_d2)
    L31 = _F31_poly(X31, Zeff_eff)

    # --- L34 ---
    L34 = if use_legacy_L34
        f33_d1 = @. 0.25 * (1 - 0.7 * fT) * sqrt(nu_e_star) * (1 + 0.45 * sqrt(Zclip))
        f33_d2 = @. (0.61 * (1 - 0.41 * fT) * nu_e_star) / sqrt(Zeff_eff)
        f33teff = @. fT / (1 + f33_d1 + f33_d2)
        _F31_poly(f33teff, Zeff_eff)
    else
        L31
    end

    # --- L32 ---
    dee_2 = @. (0.23 * (1 - 0.96 * fT) * sqrt(nu_e_star)) / sqrt(Zeff_eff)
    dee_3 = @. (0.13 * (1 - 0.38 * fT) * nu_e_star / (Zeff_eff^2)) *
              (sqrt(1 + 2 * sqrt(Zclip)) +
               fT^2 * sqrt((0.075 + 0.25 * (Zeff_eff - 1)^2) * nu_e_star))
    X32_ee = @. fT / (1 + dee_2 + dee_3)
    F32_ee = @. (0.1 + 0.6 * Zeff_eff) /
                (Zeff_eff * (0.77 + 0.63 * (1 + (Zeff_eff - 1)^1.1))) *
                (X32_ee - X32_ee^4) +
                0.7 / (1 + 0.2 * Zeff_eff) *
                (X32_ee^2 - X32_ee^4 - 1.2 * (X32_ee^3 - X32_ee^4)) +
                1.3 / (1 + 0.5 * Zeff_eff) * X32_ee^4

    dei_2 = @. (0.87 * (1 + 0.39 * fT) * sqrt(nu_e_star)) / (1 + 2.95 * (Zeff_eff - 1)^2)
    dei_3 = @. 1.53 * (1 - 0.37 * fT) * nu_e_star * (2 + 0.375 * (Zeff_eff - 1))
    X32_ei = @. fT / (1 + dei_2 + dei_3)
    F32_ei = @. -(0.4 + 1.93 * Zeff_eff) / (Zeff_eff * (0.8 + 0.6 * Zeff_eff)) *
                  (X32_ei - X32_ei^4) +
                5.5 / (1.5 + 2 * Zeff_eff) *
                  (X32_ei^2 - X32_ei^4 - 0.8 * (X32_ei^3 - X32_ei^4)) -
                1.3 / (1 + 0.5 * Zeff_eff) * X32_ei^4
    L32 = F32_ee .+ F32_ei

    # --- alpha ---
    alpha0 = @. -(0.62 + 0.055 * (Zeff_eff - 1)) /
                  (0.53 + 0.17 * (Zeff_eff - 1)) *
                  (1 - fT) /
                  (1 - (0.31 - 0.065 * (Zeff_eff - 1)) * fT - 0.25 * fT^2)
    alpha = @. ((alpha0 + 0.7 * Zeff_eff * sqrt(fT) * sqrt(nu_i_star)) /
                (1 + 0.18 * sqrt(nu_i_star)) - 0.002 * nu_i_star^2 * fT^6) /
              (1 + 0.004 * nu_i_star^2 * fT^6)

    # --- Assemble ---
    j_bootstrap, term_decomp = if formula_form == "jboot1"
        dp = dp_dpsi === nothing ?
            (@. (ne * dT_e_dpsi + Te * dn_e_dpsi + ni * dT_i_dpsi + Ti * dn_i_dpsi) * _EC) :
            dp_dpsi
        bra1 = @. L31 * dp / pe
        bra2 = @. L32 * dT_e_dpsi / Te
        bra3 = @. L34 * alpha * (1 - R_pe) / R_pe * dT_i_dpsi / Ti
        j = @. -I_psi * pe * (bra1 + bra2 + bra3)
        td = Dict{String,Any}(
            "bra1" => (@. -I_psi * pe * bra1),
            "bra2" => (@. -I_psi * pe * bra2),
            "bra3" => (@. -I_psi * pe * bra3),
        )
        (j, td)
    else
        term_n = @. p_total * L31 * (1.0 / ne) * dn_e_dpsi
        term_Te = @. pe * (L31 + L32) * (1.0 / Te) * dT_e_dpsi
        term_Ti = @. pi * (L31 + L34 * alpha) * (1.0 / Ti) * dT_i_dpsi
        j = @. -I_psi * (term_n + term_Te + term_Ti)
        td = Dict{String,Any}(
            "term_n" => (@. -I_psi * term_n),
            "term_Te" => (@. -I_psi * term_Te),
            "term_Ti" => (@. -I_psi * term_Ti),
        )
        (j, td)
    end

    if use_sign_q
        sgn = sign.(q)
        j_bootstrap = j_bootstrap .* sgn
        for k in keys(term_decomp)
            term_decomp[k] = term_decomp[k] .* sgn
        end
    end

    coeffs = Dict{String,Any}(
        "L31" => L31, "L32" => L32, "L34" => L34,
        "alpha" => alpha, "alpha0" => alpha0,
        "nu_e_star" => nu_e_star, "nu_i_star" => nu_i_star,
        "Zeff_used" => Zeff_eff,
        "R_pe" => R_pe,
        "formula_form" => formula_form,
        "F32_ee" => F32_ee, "F32_ei" => F32_ei,
    )
    merge!(coeffs, term_decomp)
    return j_bootstrap, coeffs
end

"""
    get_jphi_from_GS(ffprime, pprime, R_avg, one_over_R_avg)

Compute toroidal current density from Grad-Shafranov flux functions.
"""
get_jphi_from_GS(ffprime, pprime, R_avg, one_over_R_avg) =
    @. ffprime * (one_over_R_avg / _MU0) + R_avg * pprime

# ----------------------------------------------------------------------------
# TODO: solve_jphi, find_optimal_scale, analyze_bootstrap_edge_spike,
# solve_with_bootstrap. These compose the GS solver and require the full
# Phase 2 + 5 method set; left as a follow-up since the pure bootstrap
# building blocks above are already self-contained and testable.

end # module

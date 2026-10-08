module Bootstrap

# Port of bootstrap.py. Pure-Julia ports of:
#   - Hmode_profiles
#   - parameterize_edge_jBS
#   - calculate_ln_lambda
#   - redl_bootstrap (Redl et al. PoP 28, 022502 (2021))
#
# `analyze_bootstrap_edge_spike`, `solve_jphi`, `find_optimal_scale`, and
# `solve_with_bootstrap` couple to the Grad-Shafranov solver via the callback
# Refs wired by the top-level TokaMaker module. `solve_with_bootstrap` uses the
# native Redl model; the OMFIT-Sauter path and the curve-fit `parameterized_spike`
# are not ported (they raise).

using Distributions: SkewNormal, pdf
using Roots: find_zero, Brent

export Hmode_profiles, parameterize_edge_jBS, calculate_ln_lambda,
       redl_bootstrap, get_jphi_from_GS, solve_with_bootstrap

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
# Composer functions that combine bootstrap calculations with the GS solver.

"""
    solve_jphi(gs, ffp_prof, pp_prof, Ip_target, pax_target)

Set the toroidal-current profile (`jphi-linterp`) and pressure profile, then
solve for equilibrium. Mirrors `bootstrap.solve_jphi`.
"""
function solve_jphi(gs, ffp_prof::AbstractDict, pp_prof::AbstractDict,
                    Ip_target::Real, pax_target::Real)
    ffp_prof["type"] = "jphi-linterp"
    # Late binding to avoid Bootstrap → Core dependency at module load
    Base.invokelatest(set_targets!_anon[], gs, Ip_target, pax_target)
    Base.invokelatest(set_profiles!_anon[], gs, ffp_prof, pp_prof)
    Base.invokelatest(solve!_anon[], gs)
    return nothing
end

# Stub callbacks resolved at runtime so this module doesn't need to know the
# Tokamaker type at compile time. The TokaMaker top-level module wires these.
const set_targets!_anon = Ref{Function}()
const set_profiles!_anon = Ref{Function}()
const solve!_anon = Ref{Function}()
const get_profiles_anon = Ref{Function}()
const get_q_anon = Ref{Function}()
const get_stats_anon = Ref{Function}()
const compute_flux_integral_anon = Ref{Function}()
const sauter_fc_anon = Ref{Function}()
const psi_bounds_anon = Ref{Function}()

"""
    find_optimal_scale(gs, psi_N, pressure, ffp_prof, pp_prof, j_inductive,
                       Ip_target; psi_pad=1e-3, spike_prof=nothing,
                       find_j0=true, scale_j0=1.0, tolerance=1e-2,
                       max_iter=5, verbose=true)

Secant-method search for a scale factor `α` such that the input current
profile `α * j_inductive + spike_prof` produces the requested `j_0` (when
`find_j0=true`) or `I_p` (when `find_j0=false`) after a GS solve. Mirrors
`bootstrap.find_optimal_scale`.
"""
function find_optimal_scale(gs, psi_N::AbstractVector, pressure::AbstractVector,
                            ffp_prof::AbstractDict, pp_prof::AbstractDict,
                            j_inductive::AbstractVector, Ip_target::Real;
                            psi_pad::Real=1e-3,
                            spike_prof::Union{Nothing,AbstractVector}=nothing,
                            find_j0::Bool=true, scale_j0::Real=1.0,
                            tolerance::Real=1e-2, max_iter::Integer=5,
                            verbose::Bool=true)
    n_psi = length(psi_N)
    spike = spike_prof === nothing ? zeros(Float64, length(j_inductive)) :
            Vector{Float64}(spike_prof)
    j_ind = Vector{Float64}(j_inductive)
    last_jphi = zeros(Float64, n_psi)

    function get_j0_error(scale_val)
        matched = scale_val .* j_ind .+ spike
        ffp_prof["type"] = "jphi-linterp"; ffp_prof["y"] = matched
        solve_jphi(gs, ffp_prof, pp_prof, Ip_target, pressure[1])
        eq = gs.equilibrium
        profs = get_profiles_anon[](eq; npsi=n_psi, psi_pad=psi_pad)
        qres = get_q_anon[](eq; npsi=n_psi, psi_pad=psi_pad)
        out_jphi = get_jphi_from_GS(profs.F .* profs.Fp, profs.Pp,
                                     qres.ravgs[1, :], qres.ravgs[2, :])
        last_jphi = out_jphi
        input_j0 = scale_val * j_ind[1] + spike[1]
        output_j0 = out_jphi[1]
        diff = input_j0 - output_j0
        rel = abs(diff) / abs(output_j0 + eps(Float64))
        verbose && @info "  scale=$(round(scale_val; digits=4))  j0_in=$(round(input_j0; sigdigits=4))  j0_out=$(round(output_j0; sigdigits=4))  rel=$(round(rel; sigdigits=3))"
        return diff, rel
    end

    function get_Ip_error(scale_val)
        ffp_prof["type"] = "jphi-linterp"
        ffp_prof["y"] = scale_j0 .* j_ind .+ spike
        solve_jphi(gs, ffp_prof, pp_prof, Ip_target * scale_val, pressure[1])
        stats = get_stats_anon[](gs.equilibrium; lcfs_pad=psi_pad)
        out_Ip = stats["Ip"]
        diff = out_Ip - Ip_target
        rel = abs(diff) / max(abs(Ip_target), eps(Float64))
        verbose && @info "  scale=$(round(scale_val; digits=4))  Ip_out=$(round(out_Ip / 1e6; sigdigits=4)) MA  rel=$(round(rel; sigdigits=3))"
        return diff, rel
    end

    err_fun = find_j0 ? get_j0_error : get_Ip_error

    p0 = 1.0
    err0, rel0 = err_fun(p0)
    rel0 < tolerance && return p0, last_jphi
    p1 = if find_j0
        err0 < 0 ? 1.2 : 0.8
    else
        err0 < 0 ? 1.01 : 0.99
    end
    err1, rel1 = err_fun(p1)
    rel1 < tolerance && return p1, last_jphi

    for i in 1:max_iter
        abs(err1 - err0) < 1e-9 && break
        p_new = p1 - err1 * (p1 - p0) / (err1 - err0)
        p_new = clamp(p_new, 0.1, 5.0)
        err_new, rel_new = err_fun(p_new)
        rel_new < tolerance && return p_new, last_jphi
        p0, err0 = p1, err1
        p1, err1 = p_new, err_new
    end
    verbose && @info "find_optimal_scale: max_iter reached, returning best last value"
    return p1, last_jphi
end

# ----------------------------------------------------------------------------
# Peak detection — minimal subset of scipy.signal.find_peaks for our use.
#
# We only need: positions of local maxima (with optional height threshold).
# Returns indices into `y` (1-based).

"Find local-max indices in `y` where each peak satisfies y[i] > height."
function find_peaks(y::AbstractVector; height::Real=-Inf)
    n = length(y)
    n < 3 && return Int[]
    peaks = Int[]
    @inbounds for i in 2:n-1
        if y[i] > y[i-1] && y[i] >= y[i+1] && y[i] > height
            push!(peaks, i)
        end
    end
    return peaks
end

"Approximate full-width-at-half-max around peak index `idx`. Returns
`(width, left_idx, right_idx)` with left/right indices clamped to the array."
function peak_widths(y::AbstractVector, idx::Int; rel_height::Real=0.5)
    n = length(y)
    target = y[idx] * rel_height
    li = idx
    while li > 1 && y[li] > target
        li -= 1
    end
    ri = idx
    while ri < n && y[ri] > target
        ri += 1
    end
    return ri - li, li, ri
end

"""
    analyze_bootstrap_edge_spike(psi_N, j_BS) -> Dict

Detect the bootstrap edge spike, compute a smoothly-blended core+spike
profile via cubic Hermite interpolation. Mirrors a simplified version of
`bootstrap.analyze_bootstrap_edge_spike` — skips the optional curve_fit
parameterized-spike output (which uses scipy curve_fit and isn't strictly
required for `solve_with_bootstrap`).
"""
function analyze_bootstrap_edge_spike(psi_N::AbstractVector, j_BS::AbstractVector)
    edge_mask = psi_N .>= 0.7
    psi_edge = psi_N[edge_mask]
    j_edge = j_BS[edge_mask]
    edge_peaks = find_peaks(j_edge; height=0.0)
    if isempty(edge_peaks)
        @info "No clear bootstrap edge peak found"
        return nothing
    end
    # Pick peak closest to psi_N=1
    idx_peak_local = edge_peaks[argmax(psi_edge[edge_peaks])]
    peak_psi = psi_edge[idx_peak_local]
    peak_height = j_edge[idx_peak_local]
    # Find core minimum between psi=0.5 and the peak
    core_mask = (psi_N .> 0.5) .& (psi_N .< peak_psi)
    j_core = j_BS[core_mask]; psi_core = psi_N[core_mask]
    isempty(j_core) && return nothing
    lmin_arg = argmin(j_core)
    lmin_jBS = j_core[lmin_arg]
    fit_mask = psi_N .>= psi_core[lmin_arg]
    masked_spike = fill(lmin_jBS, length(psi_N))
    masked_spike[fit_mask] .= j_BS[fit_mask]
    # Smoothing window via cubic Hermite blend
    jBS_min_loc = psi_core[lmin_arg]
    dist_to_peak = peak_psi - jBS_min_loc
    blend_width = min(0.5 * dist_to_peak, 0.2)
    x_start = jBS_min_loc - blend_width / 2.0
    x_end = jBS_min_loc + blend_width / 2.0
    idx_start = argmin(abs.(psi_N .- x_start))
    idx_end = argmin(abs.(psi_N .- x_end))
    idx_start = max(1, idx_start)
    idx_end = min(length(psi_N) - 1, idx_end)
    y_start = masked_spike[idx_start]; dy_start = 0.0
    y_end = masked_spike[idx_end]
    dy_end = (masked_spike[idx_end + 1] - masked_spike[idx_end - 1]) /
             (psi_N[idx_end + 1] - psi_N[idx_end - 1])
    x_window = psi_N[idx_start:idx_end]
    t = (x_window .- x_window[1]) ./ (x_window[end] - x_window[1] + eps(Float64))
    h00 = @. 2 * t^3 - 3 * t^2 + 1
    h10 = @. t^3 - 2 * t^2 + t
    h01 = @. -2 * t^3 + 3 * t^2
    h11 = @. t^3 - t^2
    dx = x_window[end] - x_window[1]
    masked_spike[idx_start:idx_end] .= h00 .* y_start .+ h10 .* dx .* dy_start .+
                                        h01 .* y_end .+ h11 .* dx .* dy_end
    fwhm_idx = peak_widths(j_edge, idx_peak_local; rel_height=0.5)
    sigma_init = fwhm_idx[1] > 0 ?
        (psi_edge[min(fwhm_idx[3], length(psi_edge))] -
         psi_edge[max(fwhm_idx[2], 1)]) / 2.355 : 0.01
    return Dict{String,Any}(
        "sigma" => sigma_init,
        "background" => lmin_jBS,
        "peak_psi" => peak_psi,
        "peak_height" => peak_height,
        "masked_spike" => masked_spike,
    )
end

# ----------------------------------------------------------------------------
# Self-consistent bootstrap + Grad-Shafranov solve.

"`numpy.gradient` with unit spacing: 2nd-order central in the interior, 1st-order at the edges."
function _np_gradient(f::AbstractVector)
    n = length(f)
    g = zeros(Float64, n)
    n == 1 && return g
    g[1] = f[2] - f[1]
    g[n] = f[n] - f[n-1]
    @inbounds for i in 2:n-1
        g[i] = (f[i+1] - f[i-1]) / 2
    end
    return g
end

_nan_to_zero(x::AbstractVector) = map(v -> isnan(v) ? 0.0 : v, x)

# Match NumPy: a negative base raised to a fractional power yields NaN (which is
# later nan_to_num'd) rather than throwing, as Julia's `^` would.
_pownan(x::Real, p::Real) = x < 0 ? NaN : float(x)^p

"""
    solve_with_bootstrap(gs; ne, Te, ni, Ti, Zeff, Ip_target, inductive_jphi,
                         Zis=[1.0], scale_jBS=1.0, isolate_edge_jBS=false,
                         psi_pad=1e-3, iterations=3, parameterize_jBS=false,
                         use_OMFIT_sauter=false, verbose=true) -> Dict

Self-consistently compute the bootstrap current from H-mode profiles and iterate
the Grad-Shafranov equilibrium. Port of `bootstrap.solve_with_bootstrap` using
the native Redl (2021) model.

Profiles `ne` [m^-3], `Te`/`Ti` [eV], `ni` [m^-3], `Zeff` are sampled uniformly
in normalized psi on `[0,1]`. `inductive_jphi` (the inductive toroidal current
profile) is required.

Returns a `Dict` with `"total_j_phi"`, `"j_BS"`, `"j_inductive"`,
`"isolated_j_BS"`, `"scale_j0"`, `"scale_Ip"`.

Not supported (vs Python): `use_OMFIT_sauter` (needs the `omfit_classes` Python
package) and `parameterize_jBS` (the curve-fit `parameterized_spike` is not
ported); both raise. Diagnostic plotting is omitted.
"""
function solve_with_bootstrap(gs; ne, Te, ni, Ti, Zeff, Ip_target,
                              inductive_jphi=nothing, Zis=[1.0], scale_jBS::Real=1.0,
                              isolate_edge_jBS::Bool=false, psi_pad::Real=1e-3,
                              iterations::Integer=3, parameterize_jBS::Bool=false,
                              use_OMFIT_sauter::Bool=false, verbose::Bool=true)
    use_OMFIT_sauter && error("use_OMFIT_sauter is not supported in the Julia port (requires the omfit_classes Python package); use the native Redl model.")
    inductive_jphi === nothing && error("inductive_jphi must be specified.")

    ne = Vector{Float64}(ne); Te = Vector{Float64}(Te)
    ni = Vector{Float64}(ni); Ti = Vector{Float64}(Ti)
    Zeff = Vector{Float64}(Zeff)
    ind_jphi = Vector{Float64}(inductive_jphi)

    pressure = (_EC .* ne .* Te) .+ (_EC .* ni .* Ti)   # [Pa]
    n_psi = length(pressure)
    psi_N = collect(range(0.0, 1.0; length=n_psi))

    function calculate_profiles_and_bootstrap(include_jBS::Bool)
        eq = gs.equilibrium
        gp = get_profiles_anon[](eq; npsi=n_psi, psi_pad=psi_pad)
        f = gp.F
        sf = sauter_fc_anon[](eq; npsi=n_psi, psi_pad=psi_pad)
        ft = 1.0 .- sf.fc
        # Negative inverse-aspect-ratio can occur at degenerate edge samples;
        # NumPy turns the downstream sqrt/^1.5 into NaN (later zeroed) rather than
        # erroring, so mirror that by mapping negatives to NaN.
        eps = map(e -> e < 0 ? NaN : e, sf.r_avgs[3, :] ./ sf.r_avgs[1, :])
        qr = get_q_anon[](eq; npsi=n_psi, psi_pad=psi_pad)
        qvals = qr.q
        R_avg = qr.ravgs[1, :]

        pb = psi_bounds_anon[](gs)
        psi_range = pb[2] - pb[1]
        d_psi_eff = _np_gradient(psi_N) .* psi_range
        d_psi_eff[d_psi_eff .== 0] .= 1e-9
        pprime_local = _np_gradient(pressure) ./ d_psi_eff

        j_BS_final = zeros(Float64, n_psi)
        if include_jBS
            dn_e_dpsi = _np_gradient(ne) ./ d_psi_eff
            dT_e_dpsi = _np_gradient(Te) ./ d_psi_eff
            dn_i_dpsi = _np_gradient(ni) ./ d_psi_eff
            dT_i_dpsi = _np_gradient(Ti) ./ d_psi_eff
            ln_le, ln_lii = calculate_ln_lambda(Te, Ti, ne, ni; Zeff=Zeff,
                electron_lnLambda_model="NRL", ion_lnLambda_model="Zavg")
            Zdom = 1.0
            Zavg = ne ./ ni
            Zion = _pownan.(Zdom^2 .* Zavg .* Zeff, 0.25)
            eps15 = _pownan.(eps, 1.5)
            nu_i_star = 4.90e-18 .* abs.(qvals) .* R_avg .* ni .* Zion .^ 4 .* ln_lii ./
                        (Ti .^ 2 .* eps15)
            nu_e_star = 6.921e-18 .* abs.(qvals) .* R_avg .* ne .* Zeff .* ln_le ./
                        (Te .^ 2 .* eps15)
            # A negative collisionality (e.g. NRL ln_Λ turning negative at some
            # samples) yields NaN under NumPy's sqrt downstream; mirror that
            # rather than letting Julia's sqrt throw.
            nu_i_star = map(v -> v < 0 ? NaN : v, nu_i_star)
            nu_e_star = map(v -> v < 0 ? NaN : v, nu_e_star)
            j_BS_neo, _ = redl_bootstrap(; psi_N=psi_N, Te=Te, Ti=Ti, ne=ne, ni=ni,
                pe=_EC .* (ne .* Te), pi=_EC .* (ni .* Ti), Zeff=Zeff, R=R_avg,
                q=qvals, eps=eps, fT=ft, I_psi=f,
                dT_e_dpsi=dT_e_dpsi, dT_i_dpsi=dT_i_dpsi,
                dn_e_dpsi=dn_e_dpsi, dn_i_dpsi=dn_i_dpsi,
                ln_lambda_e=ln_le, ln_lambda_ii=ln_lii,
                nu_e_star_override=nu_e_star, nu_i_star_override=nu_i_star,
                use_legacy_L34=false, use_sign_q=true, formula_form="jboot1")
            j_BS_final = _nan_to_zero(j_BS_neo .* (R_avg ./ f))
        end

        spike_prof = zeros(Float64, n_psi)
        if include_jBS
            if isolate_edge_jBS
                parameterize_jBS && error("parameterize_jBS=true is not supported (parameterized_spike not ported; see PORT_PLAN.md).")
                res = analyze_bootstrap_edge_spike(psi_N, j_BS_final)
                masked = res === nothing ? zeros(Float64, n_psi) : Vector{Float64}(res["masked_spike"])
                spike_prof = masked .* scale_jBS
            else
                spike_prof = j_BS_final .* scale_jBS
            end
        end

        # Scale the inductive current so total Ip matches the target.
        obj(alpha) = compute_flux_integral_anon[](gs, psi_N, alpha .* ind_jphi .+ spike_prof) - Ip_target
        alpha_opt = try
            find_zero(obj, (1.0, 10.0 * Ip_target), Brent(); rtol=1e-6)
        catch
            verbose && @warn "Root scalar failed to bracket. Defaulting to alpha=1.0"
            1.0
        end
        matched_j_inductive = alpha_opt .* ind_jphi
        matched_jphi = matched_j_inductive .+ spike_prof

        pp_dict = Dict{String,Any}("type" => "linterp", "x" => copy(psi_N),
                                   "y" => pprime_local ./ pprime_local[1])
        ffp_dict = Dict{String,Any}("type" => "jphi-linterp", "x" => copy(psi_N),
                                    "y" => _nan_to_zero(matched_jphi))
        return pp_dict, ffp_dict, j_BS_final, matched_j_inductive, spike_prof
    end

    set_targets!_anon[](gs, Ip_target, pressure[1])

    verbose && @info ">>> Matching input core j_phi with G-S solution"
    pp_prof, ffp_prof, j_bs_curr, matched_j_inductive, spike_prof =
        calculate_profiles_and_bootstrap(true)
    pp_prof["y"][end] = 0.0
    ffp_prof["type"] = "jphi-linterp"
    ffp_prof["y"] = matched_j_inductive .+ spike_prof
    solve_jphi(gs, ffp_prof, pp_prof, Ip_target, pressure[1])

    pp_prof, ffp_prof, j_bs_curr, matched_j_inductive, spike_prof =
        calculate_profiles_and_bootstrap(true)
    pp_prof["y"][end] = 0.0

    verbose && @info ">>> Finding optimal j_phi scale factor"
    final_scale_j0, _ = find_optimal_scale(gs, psi_N, pressure, ffp_prof, pp_prof,
        matched_j_inductive, Ip_target; psi_pad=psi_pad, spike_prof=spike_prof,
        find_j0=true, verbose=verbose)
    verbose && @info ">>> Finding optimal Ip scale factor"
    final_scale_Ip, _ = find_optimal_scale(gs, psi_N, pressure, ffp_prof, pp_prof,
        matched_j_inductive, Ip_target; psi_pad=psi_pad, spike_prof=spike_prof,
        find_j0=false, scale_j0=final_scale_j0, tolerance=0.001, verbose=verbose)

    verbose && @info ">>> Iterating on H-mode equilibrium solution"
    tmp_jphi = _nan_to_zero(final_scale_j0 .* matched_j_inductive .+ spike_prof)
    for _ in 1:iterations
        pp_prof, ffp_prof, j_bs_curr, matched_j_inductive, spike_prof =
            calculate_profiles_and_bootstrap(true)
        pp_prof["y"][end] = 0.0
        ffp_prof["type"] = "jphi-linterp"
        ffp_prof["y"] = final_scale_j0 .* matched_j_inductive .+ spike_prof
        solve_jphi(gs, ffp_prof, pp_prof, Ip_target * final_scale_Ip, pressure[1])
        gp = get_profiles_anon[](gs.equilibrium; npsi=n_psi, psi_pad=psi_pad)
        qr = get_q_anon[](gs.equilibrium; npsi=n_psi, psi_pad=psi_pad)
        tmp_jphi = get_jphi_from_GS(gp.F .* gp.Fp, gp.Pp, qr.ravgs[1, :], qr.ravgs[2, :])
    end

    return Dict{String,Any}(
        "total_j_phi" => tmp_jphi,
        "j_BS" => j_bs_curr,
        "j_inductive" => matched_j_inductive,
        "isolated_j_BS" => spike_prof,
        "scale_j0" => final_scale_j0,
        "scale_Ip" => final_scale_Ip,
    )
end

end # module

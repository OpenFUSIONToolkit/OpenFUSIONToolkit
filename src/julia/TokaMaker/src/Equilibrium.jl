module EquilibriumModule

using ..LibPath: OFT_PATH_SLEN
using ..CInterface
using ..SettingsModule

export TokaMakerEquilibrium

mutable struct TokaMakerEquilibrium
    tmaker_ptr::Ptr{Cvoid}      # parent TokaMaker pointer (for context, not owned)
    eq_ptr::Ptr{Cvoid}
    psi_convention::Int
    F0::Float64
    np::Int
    ncoils::Int
    # Pointer caches populated lazily via tokamaker_get_refs. These point at
    # the Fortran-owned scalars/arrays so writes mutate solver state directly.
    ffp_scale_ptr::Ptr{Float64}
    p_scale_ptr::Ptr{Float64}
    o_point_ptr::Ptr{Float64}
    lim_point_ptr::Ptr{Float64}
    x_points_ptr::Ptr{Float64}
    psi_bounds_ptr::Ptr{Float64}
    diverted_ptr::Ptr{UInt8}
    has_plasma_ptr::Ptr{UInt8}
    refs_loaded::Bool
    finalized::Bool

    function TokaMakerEquilibrium(tmaker_ptr::Ptr{Cvoid}, eq_ptr::Ptr{Cvoid};
                                   psi_convention::Int=0, F0::Float64=0.0,
                                   np::Int=0, ncoils::Int=0)
        eq = new(tmaker_ptr, eq_ptr, psi_convention, F0, np, ncoils,
                 C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
                 false, false)
        finalizer(_destroy_equil!, eq)
        return eq
    end
end

"Populate the Fortran-side reference pointers from `tokamaker_get_refs`."
function load_refs!(eq::TokaMakerEquilibrium)
    eq.refs_loaded && return eq
    o_pt = Ref{Ptr{Float64}}(C_NULL)
    lim_pt = Ref{Ptr{Float64}}(C_NULL)
    x_pts = Ref{Ptr{Float64}}(C_NULL)
    div = Ref{Ptr{UInt8}}(C_NULL)
    pb = Ref{Ptr{Float64}}(C_NULL)
    ffps = Ref{Ptr{Float64}}(C_NULL)
    ps = Ref{Ptr{Float64}}(C_NULL)
    has_pl = Ref{Ptr{UInt8}}(C_NULL)
    buf = errbuf()
    c_tokamaker_get_refs(eq.eq_ptr, o_pt, lim_pt, x_pts, div, pb, ffps, ps, has_pl, buf)
    check_err(buf, "get_refs")
    eq.o_point_ptr = o_pt[]
    eq.lim_point_ptr = lim_pt[]
    eq.x_points_ptr = x_pts[]
    eq.diverted_ptr = div[]
    eq.psi_bounds_ptr = pb[]
    eq.ffp_scale_ptr = ffps[]
    eq.p_scale_ptr = ps[]
    eq.has_plasma_ptr = has_pl[]
    eq.refs_loaded = true
    return eq
end

# Direct accessors that write through the cached Fortran pointers.
function ffp_scale(eq::TokaMakerEquilibrium)
    load_refs!(eq); return unsafe_load(eq.ffp_scale_ptr)
end
function set_ffp_scale!(eq::TokaMakerEquilibrium, value::Real)
    load_refs!(eq); unsafe_store!(eq.ffp_scale_ptr, Float64(value)); return value
end
function p_scale(eq::TokaMakerEquilibrium)
    load_refs!(eq); return unsafe_load(eq.p_scale_ptr)
end
function set_p_scale!(eq::TokaMakerEquilibrium, value::Real)
    load_refs!(eq); unsafe_store!(eq.p_scale_ptr, Float64(value)); return value
end
function o_point(eq::TokaMakerEquilibrium)
    load_refs!(eq)
    return [unsafe_load(eq.o_point_ptr, i) for i in 1:2]
end
function lim_point(eq::TokaMakerEquilibrium)
    load_refs!(eq)
    return [unsafe_load(eq.lim_point_ptr, i) for i in 1:2]
end
function psi_bounds(eq::TokaMakerEquilibrium)
    load_refs!(eq)
    return [unsafe_load(eq.psi_bounds_ptr, i) for i in 1:2]
end
function diverted(eq::TokaMakerEquilibrium)
    load_refs!(eq); return unsafe_load(eq.diverted_ptr) != 0
end
function has_plasma(eq::TokaMakerEquilibrium)
    load_refs!(eq); return unsafe_load(eq.has_plasma_ptr) != 0
end

function _destroy_equil!(eq::TokaMakerEquilibrium)
    eq.finalized && return
    if eq.eq_ptr != C_NULL
        try
            buf = errbuf()
            c_tokamaker_equil_destroy(eq.eq_ptr, buf)
        catch
            # finalizer must not throw
        end
        eq.eq_ptr = C_NULL
    end
    eq.finalized = true
    return nothing
end

"""
    get_psi(eq; normalized=true)

Return the poloidal flux on mesh vertices [np]. When `normalized=true` the
values are scaled to `[0,1]` (psi_axis=0, psi_edge=1) using the spheromak
convention (psi_convention=1) or to `[1,0]` for tokamak (psi_convention=0)
matching `_core.py:abspsi_to_normalized`.
"""
function get_psi(eq::TokaMakerEquilibrium; normalized::Bool=true)
    psi = zeros(Float64, eq.np)
    psi_lim = Ref{Float64}(0.0)
    psi_max = Ref{Float64}(0.0)
    buf = errbuf()
    c_tokamaker_get_psi(eq.eq_ptr, psi, psi_lim, psi_max, buf)
    check_err(buf, "get_psi")
    if normalized
        denom = psi_max[] - psi_lim[]
        if abs(denom) > 0
            @. psi = (psi - psi_lim[]) / denom
        end
    end
    return psi
end

function set_psi!(eq::TokaMakerEquilibrium, psi::Vector{Float64};
                 update_bounds::Bool=false)
    length(psi) == eq.np || error("psi length $(length(psi)) != np=$(eq.np)")
    buf = errbuf()
    c_tokamaker_set_psi(eq.eq_ptr, psi, update_bounds, buf)
    check_err(buf, "set_psi")
    return nothing
end

function get_globals(eq::TokaMakerEquilibrium)
    Itor = Ref{Float64}(0.0)
    centroid = zeros(Float64, 2)
    vol = Ref{Float64}(0.0)
    pvol = Ref{Float64}(0.0)
    dflux = Ref{Float64}(0.0)
    tflux = Ref{Float64}(0.0)
    bp_vol = Ref{Float64}(0.0)
    buf = errbuf()
    c_tokamaker_get_globals(eq.eq_ptr, Itor, centroid, vol, pvol, dflux, tflux, bp_vol, buf)
    check_err(buf, "get_globals")
    return (Ip=Itor[], centroid=centroid, vol=vol[], pvol=pvol[],
            dflux=dflux[], tflux=tflux[], bp_vol=bp_vol[])
end

function get_refs(eq::TokaMakerEquilibrium)
    load_refs!(eq)
    o_p = [unsafe_load(eq.o_point_ptr, i) for i in 1:2]
    lim_p = [unsafe_load(eq.lim_point_ptr, i) for i in 1:2]
    # X-point storage in Fortran is gs_factory%xpts(2, MAX_XPOINTS) where
    # MAX_XPOINTS in grad_shaf.F90 is 20 (see x_loc shape used by Python).
    x_pts = zeros(Float64, 2, 20)
    for j in 1:20, i in 1:2
        x_pts[i, j] = unsafe_load(eq.x_points_ptr, (j-1)*2 + i)
    end
    pb = [unsafe_load(eq.psi_bounds_ptr, i) for i in 1:2]
    return (o_point=o_p, lim_point=lim_p, x_points=x_pts,
            diverted=unsafe_load(eq.diverted_ptr) != 0,
            psi_bounds=pb,
            ffp_scale=unsafe_load(eq.ffp_scale_ptr),
            p_scale=unsafe_load(eq.p_scale_ptr),
            has_plasma=unsafe_load(eq.has_plasma_ptr) != 0)
end

"""
    get_xpoints(eq) -> (Matrix{Float64}, Bool)

Returns `(x_points, diverted)`. `x_points` is shaped `[n_active, 2]` with rows
`(R, Z)` to match Python's `get_xpoints` convention (`x_points[i, :]` is the
i-th candidate X-point). Sentinel zero rows are dropped.
"""
function get_xpoints(eq::TokaMakerEquilibrium)
    refs = get_refs(eq)
    xpts = refs.x_points  # (2, 20) Fortran-side storage
    # Python sentinel: first slot has R<0 if no X-points found, otherwise the
    # active list ends at the first slot with R<0.
    if xpts[1, 1] < 0.0
        return zeros(Float64, 0, 2), refs.diverted
    end
    cutoff = size(xpts, 2)
    for j in 1:size(xpts, 2)
        if xpts[1, j] < 0.0
            cutoff = j - 1
            break
        end
    end
    out = Matrix{Float64}(undef, cutoff, 2)
    for i in 1:cutoff
        out[i, 1] = xpts[1, i]
        out[i, 2] = xpts[2, i]
    end
    return out, refs.diverted
end

# ----------------------------------------------------------------------------
# Inductance and stats

"""
    calc_inductance(eq, ncoils) -> (Lp, M_p_c::Vector{Float64})

Plasma self-inductance and per-coil mutual inductance. Backed by
`tokamaker_get_plasma_Lmat`, which returns a flat vector of length
`ncoils + 1`: first `ncoils` entries are mutual M_{p,coil_i}, last entry is
plasma self-inductance.
"""
function calc_inductance(eq::TokaMakerEquilibrium, ncoils::Integer)
    Lmat = zeros(Float64, ncoils + 1)
    buf = errbuf()
    c_tokamaker_get_plasma_Lmat(eq.eq_ptr, Lmat, buf)
    check_err(buf, "get_plasma_Lmat")
    return Lmat[end], Lmat[1:ncoils]
end

"""
    get_stats(eq; lcfs_pad=nothing, axis_pad=0.02, li_normalization="std",
              geom_type="max", beta_Ip=nothing, free_boundary=true,
              mirror_mode=false) -> Dict{String,Any}

Mirrors `TokaMaker_equilibrium.get_stats` in `_core.py:2357`. Returns a dict
with keys: Ip, Ip_centroid, kappa, kappaU, kappaL, delta, deltaU, deltaL,
R_geo, a_geo, vol, q_0, q_95, P_ax, P_max, W_MHD, beta_pol, dflux, tflux,
l_i, plus beta_tor and beta_n if `F0 > 0`.
"""
function get_stats(eq::TokaMakerEquilibrium;
                   lcfs_pad::Union{Nothing,Real}=nothing,
                   axis_pad::Real=0.02,
                   li_normalization::AbstractString="std",
                   geom_type::AbstractString="max",
                   beta_Ip::Union{Nothing,Real}=nothing,
                   free_boundary::Bool=true,
                   mirror_mode::Bool=false)
    g = get_globals(eq)
    Ip = beta_Ip === nothing ? g.Ip : Float64(beta_Ip)
    p_psi = collect(range(0.0, 1.0; length=100))
    p_psi[1] = 0.001
    profs = get_profiles(eq; psi=p_psi)
    p = profs.P
    if mirror_mode
        return Dict{String,Any}(
            "Ip" => Ip, "Ip_centroid" => g.centroid, "vol" => g.vol,
            "P_ax" => p[1], "P_max" => maximum(p), "W_MHD" => g.pvol * 1.5,
            "dflux" => g.dflux, "tflux" => g.tflux,
        )
    end
    pad = lcfs_pad === nothing ?
          ((diverted(eq) || !free_boundary) ? 0.01 : 0.0) :
          Float64(lcfs_pad)
    qres = get_q(eq; psi=Float64[1.0 - pad, 0.95, axis_pad], compute_geo=true)
    qvals = qres.q
    rbounds = qres.rbounds
    zbounds = qres.zbounds
    # rbounds/zbounds use the Python convention [minmax, coord] (1=min, 2=max;
    # 1=R, 2=Z) thanks to the transpose in get_q above.
    # Diverted-topology X-point handling matches Python lines 2394-2410.
    if diverted(eq)
        x_points, _ = get_xpoints(eq)
        if size(x_points, 1) > 0
            x_active = x_points[end, :]
            # Index translation: Python `x_points[-(i+1)]` for `i in 0..n-2`
            # ↔ Julia `x_points[end - i, :]` for `i in 0..n-2`. Julia 1-based
            # loop variable `j = i + 1` runs 1..n-1.
            n_xp = size(x_points, 1)
            if x_active[2] < zbounds[1, 2]
                zbounds[1, :] = x_active
                for j in 1:n_xp - 1
                    i = j - 1
                    if x_points[end - i, 2] > zbounds[2, 2]
                        zbounds[2, :] = x_points[end - i, :]; break
                    end
                end
            elseif x_active[2] > zbounds[2, 2]
                zbounds[2, :] = x_active
                for j in 1:n_xp - 1
                    i = j - 1
                    if x_points[end - i, 2] < zbounds[1, 2]
                        zbounds[1, :] = x_points[end - i, :]; break
                    end
                end
            end
        end
    end
    op = o_point(eq)
    li = if lowercase(li_normalization) == "std"
        (g.bp_vol / g.vol) / (_MU0 * Ip / qres.dl)^2
    elseif lowercase(li_normalization) == "iter"
        2.0 * g.bp_vol / ((_MU0 * Ip)^2 * op[1])
    else
        error("Invalid li_normalization \"$li_normalization\"")
    end
    R_geo, a_geo = if geom_type == "mid"
        rlcfs = trace_surf(eq, 1.0 - pad)
        rlcfs === nothing && error("trace_surf failed in get_stats")
        rlcfs = rlcfs[rlcfs[:, 1] .< op[1], :]
        iLFS = 1
        iHFS = argmin(abs.(rlcfs[:, 2] .- op[2]))
        ((rlcfs[iLFS, 1] + rlcfs[iHFS, 1]) / 2,
         (rlcfs[iLFS, 1] - rlcfs[iHFS, 1]) / 2)
    elseif geom_type == "max"
        ((rbounds[2, 1] + rbounds[1, 1]) / 2,
         (rbounds[2, 1] - rbounds[1, 1]) / 2)
    else
        error("Invalid geom_type \"$geom_type\"")
    end
    width = rbounds[2, 1] - rbounds[1, 1]
    out = Dict{String,Any}(
        "Ip" => Ip,
        "Ip_centroid" => g.centroid,
        "kappa" => (zbounds[2, 2] - zbounds[1, 2]) / width,
        "kappaU" => (zbounds[2, 2] - op[2]) * 2.0 / width,
        "kappaL" => (op[2] - zbounds[1, 2]) * 2.0 / width,
        "delta" => ((rbounds[2, 1] + rbounds[1, 1]) / 2 -
                    (zbounds[2, 1] + zbounds[1, 1]) / 2) * 2.0 / width,
        "deltaU" => ((rbounds[2, 1] + rbounds[1, 1]) / 2 - zbounds[2, 1]) * 2.0 / width,
        "deltaL" => ((rbounds[2, 1] + rbounds[1, 1]) / 2 - zbounds[1, 1]) * 2.0 / width,
        "R_geo" => R_geo,
        "a_geo" => a_geo,
        "vol" => g.vol,
        "q_0" => qvals[3],
        "q_95" => qvals[2],
        "P_ax" => p[1],
        "P_max" => maximum(p),
        "W_MHD" => g.pvol * 1.5,
        "beta_pol" => 100.0 * (2.0 * g.pvol * _MU0 / g.vol) / (Ip * _MU0 / qres.dl)^2,
        "dflux" => g.dflux,
        "tflux" => g.tflux,
        "l_i" => li,
    )
    if eq.F0 > 0.0
        beta_tor = 100.0 * (2.0 * g.pvol * _MU0 / g.vol) / (eq.F0 / R_geo)^2
        out["beta_tor"] = beta_tor
        out["beta_n"] = beta_tor * a_geo * (eq.F0 / R_geo) / (Ip / 1e6)
    end
    return out
end

# ----------------------------------------------------------------------------
# Profiles, q-profile, Sauter

const _MU0 = 4π * 1e-7

function get_profiles(eq::TokaMakerEquilibrium;
                      psi::Union{Nothing,AbstractVector}=nothing,
                      psi_pad::Real=1e-8, npsi::Integer=50)
    psi_user = psi === nothing ?
        collect(range(Float64(psi_pad), 1.0 - Float64(psi_pad); length=npsi)) :
        Vector{Float64}(psi)
    psi_call = eq.psi_convention == 0 ? (1.0 .- psi_user) : copy(psi_user)
    n = length(psi_call)
    f = zeros(Float64, n); fp = zeros(Float64, n)
    p = zeros(Float64, n); pp = zeros(Float64, n)
    buf = errbuf()
    c_tokamaker_get_profs(eq.eq_ptr, n, psi_call, f, fp, p, pp, buf)
    check_err(buf, "get_profs")
    return (psi=psi_user, F=f, Fp=fp, P=p ./ _MU0, Pp=pp ./ _MU0)
end

"""
    get_q(eq; psi=nothing, psi_pad=0.02, npsi=50, compute_geo=false)

Returns `(psi, q, ravgs, dl, rbounds, zbounds)`. With `compute_geo=false`,
`dl/rbounds/zbounds` are returned as `nothing`. `rbounds[:,1]` corresponds to
the *last* sample; the LCFS bounding box for `get_stats` is obtained by
passing `[1.0-lcfs_pad, 0.95, axis_pad]` and reading the first column.
"""
function get_q(eq::TokaMakerEquilibrium;
               psi::Union{Nothing,AbstractVector}=nothing,
               psi_pad::Real=0.02, npsi::Integer=50, compute_geo::Bool=false)
    psi_user = psi === nothing ?
        collect(range(Float64(psi_pad), 1.0 - Float64(psi_pad); length=npsi)) :
        Vector{Float64}(psi)
    psi_call = eq.psi_convention == 0 ? (1.0 .- psi_user) : copy(psi_user)
    n = length(psi_call)
    qvals = zeros(Float64, n)
    ravgs = zeros(Float64, 3, n)
    dl = Ref{Float64}(compute_geo ? 1.0 : -1.0)
    rbounds = zeros(Float64, 2, 2)
    zbounds = zeros(Float64, 2, 2)
    buf = errbuf()
    c_tokamaker_get_q(eq.eq_ptr, n, psi_call, qvals, ravgs, dl, rbounds, zbounds, buf)
    check_err(buf, "get_q")
    psi_save = psi_user
    # Transpose rbounds/zbounds so the Julia view uses Python's [i,j] indexing
    # convention: [minmax_index, coord_index]. Index 1 = min, 2 = max.
    if compute_geo
        return (psi=psi_save, q=qvals, ravgs=ravgs, dl=dl[],
                rbounds=Matrix{Float64}(rbounds'),
                zbounds=Matrix{Float64}(zbounds'))
    else
        return (psi=psi_save, q=qvals, ravgs=ravgs,
                dl=nothing, rbounds=nothing, zbounds=nothing)
    end
end

function trace_surf(eq::TokaMakerEquilibrium, psi::Real)
    psi_call = eq.psi_convention == 0 ? 1.0 - Float64(psi) : Float64(psi)
    pts = Ref{Ptr{Float64}}(C_NULL)
    npts = Ref{Int32}(0)
    buf = errbuf()
    c_tokamaker_trace_surf(eq.tmaker_ptr, psi_call, pts, npts, buf)
    check_err(buf, "trace_surf")
    npts[] <= 0 && return nothing
    raw = unsafe_wrap(Array, pts[], (2, Int(npts[])); own=false)
    return Matrix{Float64}(raw') |> copy
end

function calc_sauter_fc(eq::TokaMakerEquilibrium;
                       psi::Union{Nothing,AbstractVector}=nothing,
                       psi_pad::Real=0.02, npsi::Integer=50)
    psi_in = psi === nothing ?
        collect(range(Float64(psi_pad), 1.0 - Float64(psi_pad); length=npsi)) :
        Vector{Float64}(psi)
    n = length(psi_in)
    fc = zeros(Float64, n)
    r_avgs = zeros(Float64, 3, n)
    modb_avgs = zeros(Float64, 2, n)
    buf = errbuf()
    c_tokamaker_sauter_fc(eq.eq_ptr, n, psi_in, fc, r_avgs, modb_avgs, buf)
    check_err(buf, "sauter_fc")
    return (psi=psi_in, fc=fc, r_avgs=r_avgs, modb_avgs=modb_avgs)
end

function calc_loopvoltage(eq::TokaMakerEquilibrium)
    v = Ref{Float64}(0.0)
    buf = errbuf()
    c_tokamaker_gs_calc_vloop(eq.eq_ptr, v, buf)
    check_err(buf, "calc_vloop")
    return v[]
end

function calc_delstar_curr(eq::TokaMakerEquilibrium, psi::Vector{Float64})
    out = Vector{Float64}(psi)
    buf = errbuf()
    c_tokamaker_get_dels_curr(eq.eq_ptr, out, buf)
    check_err(buf, "get_dels_curr")
    return out
end

function calc_jtor_plasma(eq::TokaMakerEquilibrium, np::Integer)
    j = zeros(Float64, np)
    buf = errbuf()
    c_tokamaker_get_jtor(eq.eq_ptr, j, buf)
    check_err(buf, "get_jtor")
    return j
end

# ----------------------------------------------------------------------------
# I/O

const _DISABLED = -1e99

function save_eqdsk(eq::TokaMakerEquilibrium, filename::AbstractString;
                    nr::Integer=65, nz::Integer=65,
                    rbounds::Union{Nothing,AbstractVector}=nothing,
                    zbounds::Union{Nothing,AbstractVector}=nothing,
                    run_info::AbstractString="TokaMaker.jl",
                    psi_pad::Real=0.001, rcentr::Real=_DISABLED,
                    truncate_eq::Bool=false,
                    limiter_file::AbstractString="none",
                    lcfs_pressure::Real=0.0, cocos::Integer=7)
    rb = rbounds === nothing ? Float64[_DISABLED, _DISABLED] : Vector{Float64}(rbounds)
    zb = zbounds === nothing ? Float64[_DISABLED, _DISABLED] : Vector{Float64}(zbounds)
    buf = errbuf()
    c_tokamaker_save_eqdsk(eq.eq_ptr, padpath(filename), nr, nz, rb, zb,
                           padpath(run_info), Float64(psi_pad), Float64(rcentr),
                           truncate_eq, padpath(limiter_file), Float64(lcfs_pressure),
                           Int32(cocos), buf)
    check_err(buf, "save_eqdsk")
    return filename
end

function save_ifile(eq::TokaMakerEquilibrium, filename::AbstractString;
                    npsi::Integer=65, ntheta::Integer=65,
                    psi_pad::Real=0.01, lcfs_pressure::Real=0.0,
                    pack_lcfs::Bool=false, single_prec::Bool=false)
    buf = errbuf()
    c_tokamaker_save_ifile(eq.eq_ptr, padpath(filename), npsi, ntheta,
                           Float64(psi_pad), Float64(lcfs_pressure),
                           pack_lcfs, single_prec, buf)
    check_err(buf, "save_ifile")
    return filename
end

function save_mug(eq::TokaMakerEquilibrium, filename::AbstractString)
    buf = errbuf()
    c_tokamaker_save_mug(eq.eq_ptr, padpath(filename), buf)
    check_err(buf, "save_mug")
    return filename
end

end # module

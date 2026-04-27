module CoreModule

using OrderedCollections: OrderedDict
using ..LibPath: OFT_PATH_SLEN
using ..SettingsModule
using ..CInterface
using ..OFTEnvModule: OFTEnv
using ..EquilibriumModule
using ..Profiles: write_profile_file

export Tokamaker, setup_mesh!, setup_regions!, setup!, init_psi!, solve!,
       set_profiles!, set_targets!, get_targets, set_isoflux_constraints!,
       set_psi_constraints!, set_saddle_constraints!,
       set_coil_currents!, get_coil_currents,
       vac_solve!, compute_area_integral, compute_flux_integral,
       get_conductor_currents, get_conductor_source,
       reset!, update_settings!, get_psi, set_psi!,
       coil_dict2vec, coil_vec2dict, set_coil_bounds!, set_coil_vsc!, set_coil_reg!

const _MU0 = 4π * 1e-7

mutable struct Tokamaker
    env::OFTEnv
    tmaker_ptr::Ptr{Cvoid}
    mesh_ptr::Ptr{Cvoid}
    equilibrium::Union{TokaMakerEquilibrium,Nothing}
    settings::TokamakerSettings

    cond_dict::OrderedDict{String,Dict{String,Any}}
    vac_dict::OrderedDict{String,Dict{String,Any}}
    coil_dict::OrderedDict{String,Dict{String,Any}}
    coil_sets::OrderedDict{String,Dict{String,Any}}
    coil_set_names::Vector{String}
    virtual_coils::OrderedDict{String,Dict{String,Any}}
    vcoils::OrderedDict{String,Float64}
    dist_coils::OrderedDict{String,Any}

    F0::Float64
    psi_convention::Int

    nregs::Int
    ncoils::Int
    np::Int
    nc::Int
    nvac::Int
    order::Int

    r::Matrix{Float64}            # [np, 3] last col zeros
    lc::Matrix{Int32}             # [nc, 3]
    reg::Vector{Int32}            # [nc]
    Lcoils::Matrix{Float64}       # [ncoils, ncoils]
    lim_contour::Matrix{Float64}
    lim_contours::Vector{Matrix{Float64}}

    targets::OrderedDict{String,Float64}

    finalized::Bool

    function Tokamaker(env::OFTEnv)
        new(env,
            C_NULL,
            C_NULL,
            nothing,
            tokamaker_default_settings(),
            OrderedDict{String,Dict{String,Any}}(),  # cond
            OrderedDict{String,Dict{String,Any}}(),  # vac
            OrderedDict{String,Dict{String,Any}}(),  # coil
            OrderedDict{String,Dict{String,Any}}(),  # coil_sets
            String[],                                  # coil_set_names
            OrderedDict{String,Dict{String,Any}}("#VSC" => Dict{String,Any}("id" => -1, "facs" => Dict{String,Float64}())),
            OrderedDict{String,Float64}(),
            OrderedDict{String,Any}(),
            0.0, 0,
            -1, -1, -1, -1, 0, 0,
            zeros(Float64, 0, 0),
            zeros(Int32, 0, 0),
            Int32[],
            zeros(Float64, 0, 0),
            zeros(Float64, 0, 0),
            Matrix{Float64}[],
            OrderedDict{String,Float64}(),
            false,
            )
    end
end

function _make_finalizer(t::Tokamaker)
    finalizer(_destroy_tmaker!, t)
    return t
end

function _destroy_tmaker!(t::Tokamaker)
    t.finalized && return
    if t.equilibrium !== nothing
        # equilibrium has its own finalizer
        t.equilibrium = nothing
    end
    if t.tmaker_ptr != C_NULL
        try
            buf = errbuf()
            c_tokamaker_destroy(t.tmaker_ptr, buf)
        catch
        end
        t.tmaker_ptr = C_NULL
    end
    t.finalized = true
    return nothing
end

function update_settings!(t::Tokamaker)
    t.tmaker_ptr == C_NULL && return
    buf = errbuf()
    c_tokamaker_set_settings(t.tmaker_ptr, t.settings, buf)
    check_err(buf, "update_settings")
    return nothing
end

function reset!(t::Tokamaker)
    if t.tmaker_ptr != C_NULL
        buf = errbuf()
        c_tokamaker_destroy(t.tmaker_ptr, buf)
        check_err(buf, "reset")
    end
    t.tmaker_ptr = C_NULL
    t.mesh_ptr = C_NULL
    t.equilibrium = nothing
    t.settings = tokamaker_default_settings()
    empty!(t.cond_dict); empty!(t.vac_dict); empty!(t.coil_dict)
    empty!(t.coil_sets); empty!(t.coil_set_names)
    t.virtual_coils = OrderedDict{String,Dict{String,Any}}(
        "#VSC" => Dict{String,Any}("id" => -1, "facs" => Dict{String,Float64}()))
    empty!(t.vcoils); empty!(t.dist_coils)
    t.F0 = 0.0
    t.nregs = -1; t.ncoils = -1; t.np = -1; t.nc = -1; t.nvac = 0; t.order = 0
    t.r = zeros(Float64, 0, 0)
    t.lc = zeros(Int32, 0, 0)
    t.reg = Int32[]
    t.Lcoils = zeros(Float64, 0, 0)
    t.lim_contour = zeros(Float64, 0, 0)
    empty!(t.lim_contours)
    return nothing
end

"""
    setup_mesh!(tm; r, lc, reg=nothing)

Pass a 2D triangular mesh to Tokamaker. `r` is `[np, 2]` (or `[np, ndim]`),
`lc` is `[nc, 3]` 1-based connectivity, `reg` is `[nc]` 1-based region ids.
The Fortran side expects 1-based indices, matching Julia conventions.
"""
function setup_mesh!(t::Tokamaker; r::AbstractMatrix, lc::AbstractMatrix,
                     reg::Union{AbstractVector,Nothing}=nothing,
                     lc_zero_based::Bool=false)
    t.nregs == -1 || error("Mesh already setup; call reset!() first")
    # Accept Python-style user input: r [np, ndim], lc [nc, npc].
    # Fortran wants (ndim, np) and (npc, nc) layouts. In col-major Julia
    # we transpose to get the right memory order.
    # Python's gs_Domain.build_mesh / load_gs_mesh return 0-based lc; auto-detect.
    np_pts = size(r, 1)
    ndim = size(r, 2)
    nc = size(lc, 1)
    npc = size(lc, 2)
    rmat = Matrix{Float64}(transpose(r))      # (ndim, np)
    lc_local = Matrix{Int32}(lc)
    if lc_zero_based || (minimum(lc_local) == 0)
        lc_local .+= Int32(1)
    end
    lcmat = Matrix{Int32}(transpose(lc_local)) # (npc, nc)
    regvec = reg === nothing ? ones(Int32, nc) : Vector{Int32}(reg)
    minimum(regvec) <= 0 && error("Invalid 'reg' array, values must be >= 1")
    np = np_pts
    nregs_out = Ref{Int32}(0)
    mesh_ptr_out = Ref{Ptr{Cvoid}}(C_NULL)
    c_oft_setup_smesh(ndim, np, rmat, npc, nc, lcmat, regvec, nregs_out, mesh_ptr_out)
    t.mesh_ptr = mesh_ptr_out[]
    buf = errbuf()
    tmaker_ptr_out = Ref{Ptr{Cvoid}}(C_NULL)
    c_tokamaker_alloc(tmaker_ptr_out, t.mesh_ptr, buf)
    check_err(buf, "tokamaker_alloc")
    t.tmaker_ptr = tmaker_ptr_out[]
    _make_finalizer(t)
    update_settings!(t)
    t.nregs = Int(nregs_out[])
    return t
end

"""
    setup_regions!(tm; cond_dict=Dict(), coil_dict=Dict())

Configure conducting regions and coils. Each value in `cond_dict` is a Dict
with keys: `reg_id`, `eta` (Ohm-m), optional `noncontinuous`, `allow_xpoints`,
`inner_limiter`, `vac_id`. Each value in `coil_dict` has keys: `reg_id`,
optional `nturns`, `coil_set`, `allow_xpoints`.
"""
function setup_regions!(t::Tokamaker;
                        cond_dict::AbstractDict=OrderedDict{String,Dict{String,Any}}(),
                        coil_dict::AbstractDict=OrderedDict{String,Dict{String,Any}}())
    t.nregs > 0 || error("Mesh not setup; call setup_mesh! first")
    nregs = t.nregs
    xpoint_mask = zeros(Int32, nregs)
    xpoint_mask[1] = Int32(1)
    eta_vals = fill(-2.0, nregs)
    eta_vals[1] = -1.0
    contig_flag = ones(Int32, nregs)

    vac_dict = OrderedDict{String,Dict{String,Any}}()
    cond_local = OrderedDict{String,Dict{String,Any}}()
    for (key, spec) in cond_dict
        spec = Dict{String,Any}(string(k) => v for (k, v) in spec)
        rid = Int(spec["reg_id"])
        if haskey(spec, "vac_id")
            vac_dict[String(key)] = spec
            continue
        end
        eta_vals[rid] = Float64(spec["eta"]) / _MU0
        if get(spec, "noncontinuous", false)
            contig_flag[rid] = Int32(0)
        end
        xpoint_mask[rid] = Int32(get(spec, "allow_xpoints", false) ? 1 : 0)
        if get(spec, "inner_limiter", false)
            contig_flag[rid] = Int32(-1)
        end
        cond_local[String(key)] = spec
    end
    t.vac_dict = vac_dict
    t.cond_dict = cond_local

    coil_local = OrderedDict{String,Dict{String,Any}}()
    coil_sets = OrderedDict{String,Dict{String,Any}}()
    nCoils = 0
    for (key, spec) in coil_dict
        spec = Dict{String,Any}(string(k) => v for (k, v) in spec)
        rid = Int(spec["reg_id"])
        xpoint_mask[rid] = Int32(get(spec, "allow_xpoints", false) ? 1 : 0)
        eta_vals[rid] = -1.0
        coil_set = String(get(spec, "coil_set", String(key)))
        if !haskey(coil_sets, coil_set)
            coil_sets[coil_set] = Dict{String,Any}(
                "id" => nCoils,
                "net_turns" => 0.0,
                "sub_coils" => Vector{Dict{String,Any}}(),
            )
            nCoils += 1
        end
        spec["name"] = String(key)
        push!(coil_sets[coil_set]["sub_coils"], spec)
        coil_sets[coil_set]["net_turns"] += Float64(get(spec, "nturns", 1.0))
        coil_local[String(key)] = spec
    end
    t.coil_dict = coil_local
    t.coil_sets = coil_sets

    t.nvac = 0
    for i in 1:nregs
        if eta_vals[i] < -1.5
            eta_vals[i] = 1e10
            t.nvac += 1
        end
    end

    t.coil_set_names = fill("", nCoils)
    coil_nturns = zeros(Float64, nregs, nCoils)  # column-major: [nregs, nCoils] matches Fortran [ncoils, nregs] in row-major (Python view)
    # Python passes `coil_nturns` as float64[nCoils, nregs] C_CONTIGUOUS, which is row-major.
    # Fortran sees this as col-major [nregs, nCoils]. In Julia (col-major), we allocate [nregs, nCoils]
    # so the underlying memory matches what Fortran expects.
    for (key, info) in coil_sets
        cid = info["id"]::Int
        t.coil_set_names[cid + 1] = key
        for sub in info["sub_coils"]
            rid = Int(sub["reg_id"])
            coil_nturns[rid, cid + 1] = Float64(get(sub, "nturns", 1.0))
        end
    end

    coil_file = padpath("none")
    buf = errbuf()
    c_tokamaker_setup_regions(t.tmaker_ptr, coil_file, eta_vals, contig_flag,
                              xpoint_mask, coil_nturns, nCoils, buf)
    check_err(buf, "setup_regions")
    return t
end

function setup!(t::Tokamaker; order::Integer=2, F0::Real=0.0, full_domain::Bool=false)
    t.np == -1 || error("G-S instance already setup")
    update_settings!(t)
    ncoils_out = Ref{Int32}(0)
    Lmat_out = Ref{Ptr{Float64}}(C_NULL)
    buf = errbuf()
    c_tokamaker_setup(t.tmaker_ptr, Int32(order), full_domain, ncoils_out, Lmat_out, buf)
    check_err(buf, "setup")
    t.order = Int(order)
    t.F0 = Float64(F0)
    t.ncoils = Int(ncoils_out[])
    # Virtual coils slot in immediately after the real coils. Python uses
    # negative indexing in numpy; in Julia we materialize the id explicitly.
    for (i, key) in enumerate(keys(t.virtual_coils))
        t.virtual_coils[key]["id"] = t.ncoils + i - 1
    end
    if t.ncoils > 0 && Lmat_out[] != C_NULL
        t.Lcoils = unsafe_wrap(Array, Lmat_out[], (t.ncoils, t.ncoils); own=false) |> copy
    else
        t.Lcoils = zeros(Float64, max(t.ncoils, 0), max(t.ncoils, 0))
    end
    # Create equilibrium object
    eq_ptr_out = Ref{Ptr{Cvoid}}(C_NULL)
    buf = errbuf()
    c_tokamaker_equil_copy(t.tmaker_ptr, C_NULL, eq_ptr_out, buf)
    check_err(buf, "equil_copy")
    eq = TokaMakerEquilibrium(t.tmaker_ptr, eq_ptr_out[];
                              psi_convention=t.psi_convention,
                              F0=t.F0, np=0, ncoils=t.ncoils)
    buf = errbuf()
    c_tokamaker_equil_set(t.tmaker_ptr, eq.eq_ptr, buf)
    check_err(buf, "equil_set")
    t.equilibrium = eq

    # Get limiter contours
    npts = Ref{Int32}(0); nloops = Ref{Int32}(0)
    rloc = Ref{Ptr{Float64}}(C_NULL); loop_ptr = Ref{Ptr{Int32}}(C_NULL)
    buf = errbuf()
    c_tokamaker_get_limiter(t.tmaker_ptr, npts, rloc, nloops, loop_ptr, buf)
    check_err(buf, "get_limiter")
    if npts[] > 0 && rloc[] != C_NULL
        lim_pts = unsafe_wrap(Array, rloc[], (2, Int(npts[])); own=false) |> copy |> permutedims  # [npts,2]
        loops = unsafe_wrap(Array, loop_ptr[], (Int(nloops[]) + 1,); own=false) |> copy
        contours = Matrix{Float64}[]
        for i in 1:Int(nloops[])
            lo = loops[i]; hi = loops[i+1]
            seg = lim_pts[lo:hi-1, :]
            push!(contours, vcat(seg, lim_pts[lo:lo, :]))
        end
        t.lim_contours = contours
        if !isempty(contours)
            largest = argmax(size(c, 1) for c in contours)
            t.lim_contour = contours[largest]
        end
    end
    # Get plotting mesh
    np_out = Ref{Int32}(0); nc_out = Ref{Int32}(0)
    rloc2 = Ref{Ptr{Float64}}(C_NULL); lcloc = Ref{Ptr{Int32}}(C_NULL); regloc = Ref{Ptr{Int32}}(C_NULL)
    buf = errbuf()
    c_tokamaker_get_mesh(t.tmaker_ptr, np_out, rloc2, nc_out, lcloc, regloc, buf)
    check_err(buf, "get_mesh")
    t.np = Int(np_out[]); t.nc = Int(nc_out[])
    if t.np > 0
        # Fortran returns r as [3, np] in memory; numpy wraps as (np, 3) row-major (= same memory).
        # In Julia (col-major) we wrap as (3, np) and transpose for the user-facing [np, 3].
        rraw = unsafe_wrap(Array, rloc2[], (3, t.np); own=false)
        t.r = Matrix{Float64}(rraw') |> copy
        lcraw = unsafe_wrap(Array, lcloc[], (3, t.nc); own=false)
        t.lc = Matrix{Int32}(lcraw') |> copy
        t.reg = unsafe_wrap(Array, regloc[], (t.nc,); own=false) |> copy
    end
    eq.np = t.np
    return t
end

# ----------------------------------------------------------------------------
# init_psi / solve / get_psi forwarding

function init_psi!(t::Tokamaker; r0::Real=-1.0, z0::Real=0.0, a::Real=0.0,
                   kappa::Real=0.0, delta::Real=0.0,
                   curr_source::Union{Nothing,AbstractVector}=nothing)
    src = curr_source === nothing ? nothing : Vector{Float64}(curr_source)
    if src !== nothing && length(src) != t.np
        error("curr_source length must equal np=$(t.np)")
    end
    buf = errbuf()
    rhs_ptr = src === nothing ? Ptr{Float64}(C_NULL) : pointer(src)
    GC.@preserve src begin
        ccall((:tokamaker_init_psi, _libpath()), Cvoid,
              (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Ptr{Cdouble}, Ptr{UInt8}),
              t.tmaker_ptr, Float64(r0), Float64(z0), Float64(a),
              Float64(kappa), Float64(delta), rhs_ptr, buf)
    end
    check_err(buf, "init_psi")
    return t
end

_libpath() = LibPath.liboftpy[]

function solve!(t::Tokamaker; vacuum::Bool=false)
    buf = errbuf()
    c_tokamaker_solve(t.tmaker_ptr, vacuum, buf)
    check_err(buf, "solve")
    return t.equilibrium
end

"""
    vac_solve!(t; psi=nothing, rhs_source=nothing) -> equilibrium

Solve a vacuum (no-plasma) Grad-Shafranov problem with present coil
currents and an optional RHS current source. Mirrors Python
`TokaMaker.vac_solve`. Use `solve!` instead if isoflux/flux/saddle
constraints are required.
"""
function vac_solve!(t::Tokamaker;
                    psi::Union{Nothing,AbstractVector}=nothing,
                    rhs_source::Union{Nothing,AbstractVector}=nothing)
    psi_vec = if psi === nothing
        zeros(Float64, t.np)
    else
        length(psi) == t.np || error("psi length $(length(psi)) != np=$(t.np)")
        Vector{Float64}(psi)
    end
    rhs_ptr = Ptr{Float64}(C_NULL)
    rhs_vec = Float64[]
    if rhs_source !== nothing
        length(rhs_source) == t.np || error("rhs_source length $(length(rhs_source)) != np=$(t.np)")
        rhs_vec = Vector{Float64}(rhs_source)
        rhs_ptr = pointer(rhs_vec)
    end
    buf = errbuf()
    GC.@preserve psi_vec rhs_vec begin
        c_tokamaker_vac_solve(t.tmaker_ptr, psi_vec, rhs_ptr, buf)
    end
    check_err(buf, "vac_solve")
    return t.equilibrium
end

"""
    compute_area_integral(t, field; reg_mask=-1) -> Float64

Integrate `field` (length `np`) over a region (negative `reg_mask` =
whole mesh). Mirrors Python `TokaMaker.compute_area_integral`.
"""
function compute_area_integral(t::Tokamaker, field::AbstractVector;
                               reg_mask::Integer=-1)
    length(field) == t.np || error("field length $(length(field)) != np=$(t.np)")
    f = Vector{Float64}(field)
    res = Ref{Float64}(0.0)
    buf = errbuf()
    c_tokamaker_area_int(t.tmaker_ptr, f, Int32(reg_mask), res, buf)
    check_err(buf, "compute_area_integral")
    return res[]
end

"""
    get_conductor_currents(t, psi; cell_centered=false, include_Vcoils=false)
        -> (mask::BitVector, currents::Vector{Float64})

For a given `psi` (length `np`), return the toroidal current density in
the conducting regions. Mirrors Python `TokaMaker.get_conductor_currents`
which delegates to the equilibrium snapshot.

* `mask` is per-cell (length `nc`); true where the cell belongs to a
  resistive conductor (or vcoil if `include_Vcoils=true`).
* If `cell_centered=true`, `currents` is per-cell — vertex values
  averaged over each triangle. Otherwise it's the per-vertex current
  density returned by `calc_delstar_curr`.
"""
function get_conductor_currents(t::Tokamaker, psi::AbstractVector;
                                cell_centered::Bool=false,
                                include_Vcoils::Bool=false)
    eq = t.equilibrium
    eq === nothing && error("get_conductor_currents: no equilibrium available")
    length(psi) == t.np || error("psi length $(length(psi)) != np=$(t.np)")
    curr = EquilibriumModule.calc_delstar_curr(eq, Vector{Float64}(psi))
    mesh_currents = cell_centered ? zeros(Float64, t.nc) : curr
    mask = falses(t.nc)
    for (_, spec) in t.cond_dict
        eta = get(spec, "eta", -1.0)
        eta > 0 || continue
        rid = Int32(spec["reg_id"])
        @inbounds for k in 1:t.nc
            t.reg[k] == rid || continue
            mask[k] = true
            if cell_centered
                v = curr[t.lc[k, 1] + 1] + curr[t.lc[k, 2] + 1] + curr[t.lc[k, 3] + 1]
                mesh_currents[k] = v / 3.0
            end
        end
    end
    if include_Vcoils
        for (cname, cobj) in t.coil_sets
            haskey(t.vcoils, cname) || continue
            for sub in cobj["sub_coils"]
                rid = Int32(sub["reg_id"])
                @inbounds for k in 1:t.nc
                    t.reg[k] == rid || continue
                    mask[k] = true
                    if cell_centered
                        v = curr[t.lc[k, 1] + 1] + curr[t.lc[k, 2] + 1] + curr[t.lc[k, 3] + 1]
                        mesh_currents[k] = v / 3.0
                    end
                end
            end
        end
    end
    return mask, mesh_currents
end

"""
    get_conductor_source(t, dpsi_dt) -> (mask::BitVector, mesh_currents::Vector{Float64})

For a `dpsi/dt` source (e.g. from linear stability eigenvectors), return
per-cell induced toroidal current density in conducting regions. Mirrors
Python `TokaMaker.get_conductor_source`.
"""
function get_conductor_source(t::Tokamaker, dpsi_dt::AbstractVector)
    length(dpsi_dt) == t.np || error("dpsi_dt length $(length(dpsi_dt)) != np=$(t.np)")
    curr = Vector{Float64}(dpsi_dt)
    @inbounds for i in 1:t.np
        if t.r[i, 1] > 0.0
            curr[i] /= t.r[i, 1]
        end
    end
    have_nc = any(haskey(spec, "noncontinuous") for (_, spec) in t.cond_dict)
    area = zeros(Float64, t.nc)
    if have_nc
        @inbounds for k in 1:t.nc
            i1 = t.lc[k, 1] + 1; i2 = t.lc[k, 2] + 1; i3 = t.lc[k, 3] + 1
            v1r = t.r[i2, 1] - t.r[i1, 1]; v1z = t.r[i2, 2] - t.r[i1, 2]
            v2r = t.r[i3, 1] - t.r[i1, 1]; v2z = t.r[i3, 2] - t.r[i1, 2]
            area[k] = abs(v1r * v2z - v1z * v2r) / 2.0
        end
    end
    mesh_currents = zeros(Float64, t.nc)
    mask = falses(t.nc)
    for (_, spec) in t.cond_dict
        eta = get(spec, "eta", -1.0)
        eta > 0 || continue
        rid = Int32(spec["reg_id"])
        cell_idxs = Int[]
        @inbounds for k in 1:t.nc
            if t.reg[k] == rid
                v = -(curr[t.lc[k, 1] + 1] + curr[t.lc[k, 2] + 1] + curr[t.lc[k, 3] + 1]) / (3.0 * eta)
                mesh_currents[k] = v
                mask[k] = true
                push!(cell_idxs, k)
            end
        end
        if get(spec, "noncontinuous", false) && !isempty(cell_idxs)
            num = sum(mesh_currents[k] * area[k] for k in cell_idxs)
            den = sum(area[k] for k in cell_idxs)
            offset = num / den
            for k in cell_idxs
                mesh_currents[k] -= offset
            end
        end
    end
    return mask, mesh_currents
end

"""
    compute_flux_integral(t, psi_vals, field_vals) -> Float64

Integrate a flux function (defined as `(psi_vals, field_vals)` pairs)
over the plasma region. Mirrors Python
`TokaMaker_equilibrium.compute_flux_integral`, including the
`psi_convention=0` flip.
"""
function compute_flux_integral(t::Tokamaker, psi_vals::AbstractVector,
                               field_vals::AbstractVector)
    length(psi_vals) == length(field_vals) ||
        error("psi_vals and field_vals must be the same length")
    eq = t.equilibrium
    eq === nothing && error("compute_flux_integral: no equilibrium available (call solve! first)")
    pv = Vector{Float64}(psi_vals)
    fv = Vector{Float64}(field_vals)
    if t.psi_convention == 0
        pv = reverse(1.0 .- pv)
        fv = reverse(fv)
    end
    res = Ref{Float64}(0.0)
    buf = errbuf()
    c_tokamaker_flux_int(eq.eq_ptr, pv, fv, length(pv), res, buf)
    check_err(buf, "compute_flux_integral")
    return res[]
end

get_psi(t::Tokamaker; normalized::Bool=true) = EquilibriumModule.get_psi(t.equilibrium; normalized=normalized)
set_psi!(t::Tokamaker, psi::Vector{Float64}; update_bounds::Bool=false) =
    EquilibriumModule.set_psi!(t.equilibrium, psi; update_bounds=update_bounds)

# Scale forwarders matching Python's TokaMaker.ffp_scale / p_scale properties.
ffp_scale(t::Tokamaker) = EquilibriumModule.ffp_scale(t.equilibrium)
set_ffp_scale!(t::Tokamaker, v::Real) = EquilibriumModule.set_ffp_scale!(t.equilibrium, v)
p_scale(t::Tokamaker) = EquilibriumModule.p_scale(t.equilibrium)
set_p_scale!(t::Tokamaker, v::Real) = EquilibriumModule.set_p_scale!(t.equilibrium, v)
o_point(t::Tokamaker) = EquilibriumModule.o_point(t.equilibrium)
lim_point(t::Tokamaker) = EquilibriumModule.lim_point(t.equilibrium)
psi_bounds(t::Tokamaker) = EquilibriumModule.psi_bounds(t.equilibrium)
diverted(t::Tokamaker) = EquilibriumModule.diverted(t.equilibrium)

function _bounds_from_limiter(t::Tokamaker, dim::Int)
    isempty(t.lim_contour) && error("save_eqdsk: lim_contour empty; pass rbounds/zbounds explicitly")
    lo, hi = extrema(@view t.lim_contour[:, dim])
    pad = (hi - lo) * 0.05
    return [lo - pad, hi + pad]
end

function EquilibriumModule.save_eqdsk(t::Tokamaker, filename::AbstractString;
                                      rbounds::Union{Nothing,AbstractVector}=nothing,
                                      zbounds::Union{Nothing,AbstractVector}=nothing,
                                      kwargs...)
    rb = rbounds === nothing ? _bounds_from_limiter(t, 1) : rbounds
    zb = zbounds === nothing ? _bounds_from_limiter(t, 2) : zbounds
    return EquilibriumModule.save_eqdsk(t.equilibrium, filename;
                                        rbounds=rb, zbounds=zb, kwargs...)
end

# ----------------------------------------------------------------------------
# Profiles

function set_profiles!(t::Tokamaker;
                       ffp_prof::Union{Nothing,AbstractDict}=nothing,
                       foffset::Union{Nothing,Real}=nothing,
                       pp_prof::Union{Nothing,AbstractDict}=nothing,
                       ffp_NI_prof::Union{Nothing,AbstractDict}=nothing,
                       keep_files::Bool=false)
    t.equilibrium === nothing && error("Equilibrium not initialized; call setup! first")
    f_file = "none"; p_file = "none"; eta_file = "none"; f_NI_file = "none"
    tmpfiles = String[]
    if ffp_prof !== nothing
        f_file = "tokamaker_f.prof"
        write_profile_file(f_file, ffp_prof, "FF'"; psi_convention=t.psi_convention)
        push!(tmpfiles, f_file)
    end
    if pp_prof !== nothing
        p_file = "tokamaker_p.prof"
        write_profile_file(p_file, pp_prof, "P'"; psi_convention=t.psi_convention)
        push!(tmpfiles, p_file)
    end
    if ffp_NI_prof !== nothing
        f_NI_file = "tokamaker_f_NI.prof"
        write_profile_file(f_NI_file, ffp_NI_prof, "FF'_NI"; psi_convention=t.psi_convention)
        push!(tmpfiles, f_NI_file)
    end
    f_off = foffset === nothing ? t.F0 : Float64(foffset)
    buf = errbuf()
    c_tokamaker_load_profiles(t.equilibrium.eq_ptr,
                              padpath(f_file), f_off,
                              padpath(p_file), padpath(eta_file),
                              padpath(f_NI_file), buf)
    check_err(buf, "load_profiles")
    if !keep_files
        for f in tmpfiles
            isfile(f) && rm(f; force=true)
        end
    end
    return t
end

# ----------------------------------------------------------------------------
# Targets and constraints

function set_targets!(t::Tokamaker;
                      Ip::Union{Nothing,Real}=nothing,
                      Ip_ratio::Union{Nothing,Real}=nothing,
                      pax::Union{Nothing,Real}=nothing,
                      estore::Union{Nothing,Real}=nothing,
                      R0::Union{Nothing,Real}=nothing,
                      V0::Union{Nothing,Real}=nothing,
                      Z0::Union{Nothing,Real}=nothing,
                      retain_previous::Bool=false)
    DISABLED = -1e99
    # Z0 and V0 are aliases (Z0 in some docs; V0 in C ABI)
    v0 = V0 !== nothing ? Float64(V0) : (Z0 !== nothing ? Float64(Z0) : nothing)
    if !retain_previous
        empty!(t.targets)
    end
    Ip       === nothing || (t.targets["Ip"]       = Float64(Ip))
    Ip_ratio === nothing || (t.targets["Ip_ratio"] = Float64(Ip_ratio))
    pax      === nothing || (t.targets["pax"]      = Float64(pax))
    estore   === nothing || (t.targets["estore"]   = Float64(estore))
    R0       === nothing || (t.targets["R0"]       = Float64(R0))
    if v0 !== nothing
        t.targets["V0"] = v0
        t.targets["Z0"] = v0
    end
    buf = errbuf()
    c_tokamaker_set_targets(t.tmaker_ptr,
        get(t.targets, "Ip",       DISABLED),
        get(t.targets, "Ip_ratio", DISABLED),
        get(t.targets, "pax",      DISABLED),
        get(t.targets, "estore",   DISABLED),
        get(t.targets, "R0",       DISABLED),
        get(t.targets, "V0",       DISABLED),
        buf)
    check_err(buf, "set_targets")
    return t
end

"Return the most recently set global targets as a Dict (mirrors Python `get_targets`)."
get_targets(t::Tokamaker) = copy(t.targets)

function set_isoflux_constraints!(t::Tokamaker, isoflux::AbstractMatrix;
                                   weights::Union{Nothing,AbstractVector}=nothing,
                                   grad_wt_lim::Real=-1.0,
                                   ref_points::Union{Nothing,AbstractMatrix}=nothing)
    size(isoflux, 2) == 2 || error("isoflux must be [n,2]")
    iso = Matrix{Float64}(isoflux)
    w = weights === nothing ? nothing : Vector{Float64}(weights)
    if ref_points === nothing
        # First point acts as the reference for all subsequent points.
        n_in = size(iso, 1)
        n_in >= 2 || error("Need >= 2 isoflux points when ref_points unspecified")
        refs = zeros(Float64, n_in - 1, 2)
        @views refs[:, 1] .= iso[1, 1]
        @views refs[:, 2] .= iso[1, 2]
        iso = iso[2:end, :]
        w = w === nothing ? nothing : w[2:end]
    else
        refs = Matrix{Float64}(ref_points)
        size(refs, 1) == size(iso, 1) || error("ref_points first dim must match isoflux")
    end
    n = size(iso, 1)
    if w === nothing
        w = ones(Float64, n)
    end
    length(w) == n || error("weights length mismatch")
    # Fortran wants F-order shape (2, n); Julia col-major matrix (2, n) provides that layout.
    iso_f = Matrix{Float64}(iso')
    refs_f = Matrix{Float64}(refs')
    buf = errbuf()
    c_tokamaker_set_isoflux(t.tmaker_ptr, iso_f, refs_f, w, n, Float64(grad_wt_lim), buf)
    check_err(buf, "set_isoflux")
    return t
end

function set_psi_constraints!(t::Tokamaker, locations::AbstractMatrix,
                               targets::AbstractVector;
                               weights::Union{Nothing,AbstractVector}=nothing)
    n = size(locations, 1)
    size(locations, 2) == 2 || error("locations must be [n,2]")
    length(targets) == n || error("targets length mismatch")
    locs = Matrix{Float64}(locations')  # Fortran [2, n]
    w = weights === nothing ? ones(Float64, n) : Vector{Float64}(weights)
    buf = errbuf()
    c_tokamaker_set_flux(t.tmaker_ptr, locs, Vector{Float64}(targets), w, n, -1.0, buf)
    check_err(buf, "set_flux")
    return t
end

function set_saddle_constraints!(t::Tokamaker, saddles::AbstractMatrix;
                                  weights::Union{Nothing,AbstractVector}=nothing)
    n = size(saddles, 1)
    size(saddles, 2) == 2 || error("saddles must be [n,2]")
    locs = Matrix{Float64}(saddles')
    w = weights === nothing ? ones(Float64, n) : Vector{Float64}(weights)
    buf = errbuf()
    c_tokamaker_set_saddles(t.tmaker_ptr, locs, w, n, buf)
    check_err(buf, "set_saddles")
    return t
end

# ----------------------------------------------------------------------------
# Coil management

function coil_dict2vec(t::Tokamaker, coil_dict::AbstractDict;
                       keep_virtual::Bool=false, default_value::Real=0.0)
    nv = length(t.virtual_coils)
    n = t.ncoils + nv
    vec = fill(Float64(default_value), n)
    rem = zeros(Float64, n)
    for (key, val) in coil_dict
        skey = String(key)
        if haskey(t.coil_sets, skey)
            id = t.coil_sets[skey]["id"]::Int
            vec[id+1] += Float64(val)
            rem[id+1] = Float64(default_value)
        elseif haskey(t.virtual_coils, skey)
            if keep_virtual
                id = t.virtual_coils[skey]["id"]::Int
                vec[id+1] += Float64(val)
                rem[id+1] = Float64(default_value)
            else
                facs = get(t.virtual_coils[skey], "facs", Dict{String,Float64}())
                for (mk, mv) in facs
                    id = t.coil_sets[mk]["id"]::Int
                    vec[id+1] += mv * Float64(val)
                    rem[id+1] = Float64(default_value)
                end
            end
        else
            error("Unknown coil \"$skey\"")
        end
    end
    vec .-= rem
    return keep_virtual ? vec : vec[1:t.ncoils]
end

function coil_vec2dict(t::Tokamaker, coil_vec::AbstractVector;
                       always_virtual::Bool=false)
    n = length(coil_vec)
    nv = length(t.virtual_coils)
    (n == t.ncoils || n == t.ncoils + nv) || error("vector length must be $(t.ncoils) or $(t.ncoils+nv)")
    d = OrderedDict{String,Float64}()
    for (key, info) in t.coil_sets
        d[key] = coil_vec[info["id"]::Int + 1]
    end
    if n > t.ncoils
        for (key, info) in t.virtual_coils
            d[key] = coil_vec[info["id"]::Int + 1]
        end
    elseif always_virtual
        for key in keys(t.virtual_coils)
            d[key] = 0.0
        end
    end
    return d
end

function set_coil_currents!(t::Tokamaker, currents::AbstractDict)
    vec = coil_dict2vec(t, currents)
    buf = errbuf()
    c_tokamaker_set_coil_currents(t.tmaker_ptr, Vector{Float64}(vec), buf)
    check_err(buf, "set_coil_currents")
    return t
end

function get_coil_currents(t::Tokamaker)
    t.equilibrium === nothing && error("No equilibrium")
    currents = zeros(Float64, t.ncoils)
    reg_currents = zeros(Float64, max(t.nregs, 1))
    buf = errbuf()
    c_tokamaker_get_coil_currents(t.equilibrium.eq_ptr, currents, reg_currents, buf)
    check_err(buf, "get_coil_currents")
    return coil_vec2dict(t, currents), reg_currents
end

function set_coil_bounds!(t::Tokamaker, coil_bounds::Union{Nothing,AbstractDict}=nothing)
    nv = length(t.virtual_coils)
    n = t.ncoils + nv
    bounds = zeros(Float64, 2, n)
    @views bounds[1, :] .= -1e98
    @views bounds[2, :] .= 1e98
    if coil_bounds !== nothing
        for (key, bd) in coil_bounds
            skey = String(key)
            if haskey(t.coil_sets, skey)
                id = t.coil_sets[skey]["id"]::Int
                bounds[1, id+1] = Float64(bd[1]); bounds[2, id+1] = Float64(bd[2])
            elseif haskey(t.virtual_coils, skey)
                id = t.virtual_coils[skey]["id"]::Int
                bounds[1, id+1] = Float64(bd[1]); bounds[2, id+1] = Float64(bd[2])
            else
                error("Unknown coil \"$skey\"")
            end
        end
    end
    buf = errbuf()
    c_tokamaker_set_coil_bounds(t.tmaker_ptr, bounds, buf)
    check_err(buf, "set_coil_bounds")
    return t
end

function set_coil_vsc!(t::Tokamaker, coil_gains::AbstractDict)
    gains = coil_dict2vec(t, coil_gains)
    buf = errbuf()
    c_tokamaker_set_coil_vsc(t.tmaker_ptr, Vector{Float64}(gains), buf)
    check_err(buf, "set_coil_vsc")
    t.virtual_coils["#VSC"]["facs"] = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in coil_gains)
    return t
end

function set_coil_reg!(t::Tokamaker;
                       reg_mat::Union{Nothing,AbstractMatrix}=nothing,
                       reg_targets::Union{Nothing,AbstractVector}=nothing,
                       reg_weights::Union{Nothing,AbstractVector}=nothing,
                       reg_terms::Union{Nothing,AbstractVector}=nothing)
    # Fortran array layout: reg_mat(nregularize, ncoils+1) F-order. In Julia
    # (col-major), use `rmat[i, coil_idx]` where `i` is term index.
    if reg_terms !== nothing
        reg_mat === nothing || error("reg_terms and reg_mat are mutually exclusive")
        reg_targets === nothing || error("reg_terms and reg_targets are mutually exclusive")
        reg_weights === nothing || error("reg_terms and reg_weights are mutually exclusive")
        nreg = length(reg_terms)
        rmat = zeros(Float64, nreg, t.ncoils + 1)
        rtargets = zeros(Float64, nreg)
        rweights = ones(Float64, nreg)
        for (i, term) in enumerate(reg_terms)
            rtargets[i] = Float64(term.target)
            rweights[i] = Float64(term.weight)
            for (key, val) in term.coffs
                skey = String(key)
                if haskey(t.coil_sets, skey)
                    rmat[i, t.coil_sets[skey]["id"]::Int + 1] = Float64(val)
                elseif haskey(t.virtual_coils, skey)
                    rmat[i, t.virtual_coils[skey]["id"]::Int + 1] = Float64(val)
                else
                    error("Unknown coil \"$skey\"")
                end
            end
        end
    elseif reg_mat !== nothing
        size(reg_mat, 2) == t.ncoils + 1 || error("reg_mat must be [nregularize, ncoils+1]")
        nreg = size(reg_mat, 1)
        rmat = Matrix{Float64}(reg_mat)
        rtargets = reg_targets === nothing ? zeros(Float64, nreg) : Vector{Float64}(reg_targets)
        rweights = reg_weights === nothing ? ones(Float64, nreg) : Vector{Float64}(reg_weights)
    else
        error("Either reg_terms or reg_mat must be provided")
    end
    # Ensure VSC column has a regularization term if no VSC mapping is defined.
    vsc_facs = get(get(t.virtual_coils, "#VSC", Dict{String,Any}()), "facs", Dict{String,Float64}())
    if isempty(vsc_facs) && maximum(abs.(@view rmat[:, end])) < 1e-8
        nreg += 1
        new_row = zeros(Float64, 1, t.ncoils + 1); new_row[1, end] = 1.0
        rmat = vcat(rmat, new_row)
        push!(rtargets, 0.0); push!(rweights, 1.0)
    end
    buf = errbuf()
    c_tokamaker_set_coil_regmat(t.tmaker_ptr, nreg, rmat, rtargets, rweights, buf)
    check_err(buf, "set_coil_regmat")
    return t
end

# Local imports for the GC.@preserve init_psi path
import ..LibPath

end # module

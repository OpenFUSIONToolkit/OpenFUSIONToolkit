module Util

using Printf
using ..LibPath: liboftpy

export create_isoflux, create_power_flux_fun, create_spline_flux_fun,
       eval_green, read_eqdsk, read_ifile,
       read_fortran_namelist, read_mhdin, read_kfile

const _MU0 = 4π * 1e-7

"""
    create_isoflux(npts, r0, z0, a, kappa, delta; kappaL=kappa, deltaL=delta)

Generate `npts` analytic isoflux points on an elongated/triangular contour
following the Miller-form expression used by `_core.create_isoflux`.
Returns an `[npts, 2]` matrix.
"""
function create_isoflux(npts::Integer, r0::Real, z0::Real, a::Real, kappa::Real,
                        delta::Real; kappaL::Union{Nothing,Real}=nothing,
                        deltaL::Union{Nothing,Real}=nothing)
    kU = Float64(kappa); dU = Float64(delta)
    kL = kappaL === nothing ? kU : Float64(kappaL)
    dL = deltaL === nothing ? dU : Float64(deltaL)
    pts = zeros(Float64, npts, 2)
    for i in 0:(npts-1)
        θ = 2π * i / npts
        d = 0.5 * ((dU + dL) + (dU - dL) * sin(θ))
        k = 0.5 * ((kU + kL) + (kU - kL) * sin(θ))
        pts[i+1, 1] = r0 + a * cos(θ + asin(d) * sin(θ))
        pts[i+1, 2] = z0 + a * k * sin(θ)
    end
    return pts
end

"""
    create_power_flux_fun(npts, alpha, gamma)

Power-law profile `(1 - psi^alpha)^gamma`. Returns a profile dict suitable for
[`set_profiles!`](@ref).
"""
function create_power_flux_fun(npts::Integer, alpha::Real, gamma::Real)
    psi = collect(range(0.0, 1.0; length=npts))
    y = @. (1.0 - psi^alpha)^gamma
    return Dict("type" => "linterp", "x" => psi, "y" => y)
end

"""
    create_spline_flux_fun(npts, x, y; axis_bc, edge_bc, normalize=true)

Cubic-spline-interpolated profile sampled at `npts` evenly-spaced normalized
psi points. `axis_bc` and `edge_bc` are SciPy-style boundary conditions:
`(1, 0.0)` means clamped first derivative = 0 at that boundary.
"""
function create_spline_flux_fun(npts::Integer, x::AbstractVector, y::AbstractVector;
                                axis_bc=(1, 0.0), edge_bc=(1, 0.0), normalize::Bool=true)
    @assert length(x) == length(y)
    xs = collect(Float64.(x))
    ys = collect(Float64.(y))
    # Solve clamped/natural cubic spline manually: simple natural-spline tridiagonal solve
    n = length(xs)
    h = diff(xs)
    α = zeros(Float64, n)
    @inbounds for i in 2:n-1
        α[i] = 3 * ((ys[i+1] - ys[i]) / h[i] - (ys[i] - ys[i-1]) / h[i-1])
    end
    l = zeros(Float64, n); μ = zeros(Float64, n); z = zeros(Float64, n)
    if axis_bc[1] == 1
        l[1] = 2 * h[1]; μ[1] = 0.5
        z[1] = (3 * (ys[2] - ys[1]) / h[1] - 3 * Float64(axis_bc[2])) / l[1]
    else
        l[1] = 1.0; μ[1] = 0.0; z[1] = 0.0
    end
    @inbounds for i in 2:n-1
        l[i] = 2 * (xs[i+1] - xs[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end
    c = zeros(Float64, n); b = zeros(Float64, n - 1); d = zeros(Float64, n - 1)
    if edge_bc[1] == 1
        l[n] = h[n-1] * (2 - μ[n-1])
        z[n] = (3 * Float64(edge_bc[2]) - 3 * (ys[n] - ys[n-1]) / h[n-1] - h[n-1] * z[n-1]) / l[n]
        c[n] = z[n]
    else
        l[n] = 1.0; z[n] = 0.0; c[n] = 0.0
    end
    @inbounds for j in n-1:-1:1
        c[j] = z[j] - μ[j] * c[j+1]
        b[j] = (ys[j+1] - ys[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])
    end
    function eval_at(xq)
        idx = clamp(searchsortedlast(xs, xq), 1, n - 1)
        dx = xq - xs[idx]
        return ys[idx] + b[idx] * dx + c[idx] * dx^2 + d[idx] * dx^3
    end
    x_sample = collect(range(0.0, 1.0; length=npts))
    y_sample = [eval_at(xq) for xq in x_sample]
    if normalize
        y_sample ./= y_sample[1]
    end
    return Dict("type" => "linterp", "x" => x_sample, "y" => y_sample)
end

"""
    eval_green(x::Matrix{Float64}, xc::Vector{Float64}) -> Vector{Float64}

Toroidal-filament Green's function. `x` is `[n, 2]`, `xc` is `[2]`. Returns
`psi` due to a unit-current filament at `xc`, evaluated at each row of `x`.
"""
function eval_green(x::AbstractMatrix, xc::AbstractVector)
    n = size(x, 1)
    size(x, 2) == 2 || error("x must be [n, 2]")
    length(xc) == 2 || error("xc must have length 2")
    r = Vector{Float64}(x[:, 1])
    z = Vector{Float64}(x[:, 2])
    vals = zeros(Float64, n)
    ccall((:tokamaker_eval_green, liboftpy[]), Cvoid,
          (Cint, Ptr{Float64}, Ptr{Float64}, Cdouble, Cdouble, Ptr{Float64}),
          Int32(n), r, z, Float64(xc[1]), Float64(xc[2]), vals)
    return vals .* _MU0
end

# ----------------------------------------------------------------------------
# EQDSK reader (gEQDSK format) — port of util.read_eqdsk

function read_eqdsk(filename::AbstractString)
    obj = Dict{String,Any}()
    open(filename, "r") do io
        line = readline(io)
        obj["case"] = line[1:min(48, length(line))]
        rest = strip(line[min(49, length(line)+1):end])
        toks = split(rest)
        obj["nr"] = parse(Int, toks[end-1])
        obj["nz"] = parse(Int, toks[end])

        line_keys = [
            ["rdim",  "zdim",  "rcentr",  "rleft",  "zmid"],
            ["raxis", "zaxis", "psimag", "psibry", "bcentr"],
            ["ip",    "skip",  "skip",   "skip",   "skip"],
            ["skip",  "skip",  "skip",   "skip",   "skip"],
        ]
        for i in 1:4
            line = readline(io)
            for j in 1:5
                key = line_keys[i][j]
                key == "skip" && continue
                seg = line[(j-1)*16+1 : min(j*16, length(line))]
                obj[key] = parse(Float64, strip(seg))
            end
        end

        nr = obj["nr"]; nz = obj["nz"]
        function read_1d(n)
            out = zeros(Float64, n); buf = ""; pos = 1
            for i in 1:n
                if pos > length(buf)
                    buf = readline(io); pos = 1
                end
                out[i] = parse(Float64, strip(buf[pos:min(pos+15, length(buf))]))
                pos += 16
                pos > 16 * 5 && (pos = length(buf) + 1)
            end
            return out
        end
        function read_2d(n, m)
            out = zeros(Float64, n, m); buf = ""; pos = 1
            for k in 1:n, i in 1:m
                if pos > length(buf)
                    buf = readline(io); pos = 1
                end
                out[k, i] = parse(Float64, strip(buf[pos:min(pos+15, length(buf))]))
                pos += 16
                pos > 16 * 5 && (pos = length(buf) + 1)
            end
            return out
        end
        for k in ("fpol", "pres", "ffprim", "pprime")
            obj[k] = read_1d(nr)
        end
        obj["psirz"] = read_2d(nz, nr)
        obj["qpsi"] = read_1d(nr)
        line = readline(io)
        toks = split(strip(line))
        obj["nbbs"] = parse(Int, toks[1])
        obj["nlim"] = parse(Int, toks[2])
        obj["rzout"] = read_2d(obj["nbbs"], 2)
        obj["rzlim"] = read_2d(obj["nlim"], 2)
    end
    return obj
end

# ----------------------------------------------------------------------------
# Binary i-file reader

function _read_array(io::IO, var_type::Char, count::Int)
    var_size = (var_type in ('i', 'f')) ? 4 : (var_type in ('l', 'd')) ? 8 : error("Invalid type")
    array_size = var_size * count
    head = read(io, Int32)
    head == array_size || error("Frame size mismatch (got $head, expected $array_size)")
    if var_type == 'd'
        body = Vector{Float64}(undef, count); read!(io, body)
    elseif var_type == 'f'
        body = Vector{Float32}(undef, count); read!(io, body); body = Float64.(body)
    elseif var_type == 'i'
        body = Vector{Int32}(undef, count); read!(io, body)
    elseif var_type == 'l'
        body = Vector{Int64}(undef, count); read!(io, body)
    end
    tail = read(io, Int32)
    head == tail || error("Frame head/tail mismatch")
    return body
end

function read_ifile(filename::AbstractString)
    open(filename, "r") do io
        sizes = _read_array(io, 'i', 2)
        npsi = Int(sizes[1]); ntheta = Int(sizes[2])
        out = Dict{String,Any}("npsi" => npsi, "ntheta" => ntheta)
        var_type = 'd'
        try
            out["psi"] = _read_array(io, 'd', npsi)
        catch
            seek(io, position(io) - 4)
            try
                var_type = 'f'
                out["psi"] = _read_array(io, 'f', npsi)
            catch
                error("Unable to determine floating-point type")
            end
        end
        out["f"] = _read_array(io, var_type, npsi)
        out["p"] = _read_array(io, var_type, npsi)
        out["q"] = _read_array(io, var_type, npsi)
        R = _read_array(io, var_type, npsi * ntheta)
        Z = _read_array(io, var_type, npsi * ntheta)
        out["R"] = reshape(R, npsi, ntheta)
        out["Z"] = reshape(Z, npsi, ntheta)
        return out
    end
end

# ----------------------------------------------------------------------------
# Fortran namelist reader + EFIT machine/k-file parsers (port of
# OpenFUSIONToolkit.util.read_fortran_namelist and TokaMaker.util.read_mhdin /
# read_kfile). The `b_arr` (bottom-array) path is not ported.

"""
    read_fortran_namelist(file; silent=true) -> Dict{String,Any}

Parse a Fortran namelist file into a dict keyed by variable name. Numeric
values become `Float64`/`Int` (scalars) or `Vector{Float64}` (arrays); the
`N*value` repetition syntax is expanded; non-numeric values are kept as the raw
`String`. Mirrors `read_fortran_namelist(file, b_arr=False)`.
"""
function read_fortran_namelist(file::AbstractString; silent::Bool=true)
    lines = readlines(file)
    # End the file at a trailing "comment:" marker, if present.
    end_idx = length(lines) + 1
    for (i, ln) in enumerate(lines)
        if occursin("comment:", ln)
            end_idx = i
            break
        end
    end
    # Strip commas/quotes/comments; drop header (&) and break (/) lines.
    datalines = String[]
    for i in 1:(end_idx-1)
        lu = replace(lines[i], "," => " ")
        lu = replace(lu, "\"" => "")
        cpos = findfirst('!', lu)
        bpos = findfirst('/', lu)
        hpos = findfirst('&', lu)
        (hpos === nothing && bpos === nothing) || continue
        if cpos === nothing
            s = strip(lu)
            isempty(s) || push!(datalines, String(s))
        else
            s = strip(lu[1:cpos-1])
            isempty(s) || push!(datalines, String(s))
        end
    end
    # Identify "key = ..." block starts.
    block_idx = Int[]
    keys_list = String[]
    for (i, dl) in enumerate(datalines)
        eq = findfirst('=', dl)
        if eq !== nothing && eq > 1
            push!(block_idx, i)
            push!(keys_list, String(strip(dl[1:eq-1])))
        end
    end
    data = Dict{String,Any}()
    nkeys = length(keys_list)
    for i in 1:nkeys
        stop = i < nkeys ? block_idx[i+1] - 1 : length(datalines)
        block = join(datalines[block_idx[i]:stop], " ")
        eq = findfirst('=', block)
        rhs = strip(block[eq+1:end])
        data[keys_list[i]] = _parse_namelist_value(String(rhs))
    end
    if !silent
        for (k, v) in data
            println(k, " = ", v)
        end
    end
    return data
end

# Expand `N*value` repetitions then parse a namelist RHS to Float64 vector /
# scalar / raw string (matching the Python type inference).
function _parse_namelist_value(rhs::AbstractString)
    toks = split(rhs)
    if any(occursin('*', t) for t in toks)
        expanded = String[]
        for t in toks
            if occursin('*', t)
                parts = split(t, '*')
                n = tryparse(Int, parts[1])
                if n === nothing || length(parts) < 2
                    push!(expanded, t)        # leave as-is (matches Python `continue`)
                else
                    for _ in 1:n
                        push!(expanded, parts[2])
                    end
                end
            else
                push!(expanded, t)
            end
        end
        nums = tryparse.(Float64, expanded)
        return any(isnothing, nums) ? rhs : Float64.(nums)
    end
    if length(toks) > 1
        nums = tryparse.(Float64, toks)
        return any(isnothing, nums) ? rhs : Float64.(nums)
    else
        s = toks[1]
        iv = tryparse(Int, s)
        iv !== nothing && return iv
        fv = tryparse(Float64, s)
        fv !== nothing && return fv
        return rhs
    end
end

"""
    read_mhdin(path; e_coil_names=nothing, f_coil_names=nothing) -> (machine_dict, raw)

Read an EFIT `mhdin.dat` machine-description file. Returns `machine_dict`
(coil coordinates/turns keyed by name, probe/loop names, probe angles) and the
raw namelist dict. Mirrors `TokaMaker.util.read_mhdin`.
"""
function read_mhdin(path::AbstractString; e_coil_names=nothing, f_coil_names=nothing)
    raw = read_fortran_namelist(path)
    machine = Dict{String,Any}()
    for key in ("MPNAM2", "LPNAME")
        names = split(replace(String(raw[key]), "'" => " "))
        machine[key] = String.(names)
    end
    ecid = raw["ECID"]; RE = raw["RE"]; ZE = raw["ZE"]; WE = raw["WE"]; HE = raw["HE"]
    e_coil = Dict{String,Vector{Vector{Float64}}}()
    e_order = String[]
    for i in 1:length(ecid)
        idx = Int(ecid[i]) - 1
        name = e_coil_names === nothing ? @sprintf("ECOIL%03d", idx + 1) : e_coil_names[idx+1]
        if !haskey(e_coil, name)
            e_coil[name] = Vector{Float64}[]
            push!(e_order, name)
        end
        push!(e_coil[name], Float64[RE[i], ZE[i], WE[i], HE[i]])
    end
    machine["ECOIL"] = e_coil
    machine["ECOIL_order"] = e_order
    fcid = raw["FCID"]; RF = raw["RF"]; ZF = raw["ZF"]; WF = raw["WF"]; HF = raw["HF"]; TURNFC = raw["TURNFC"]
    f_coil = Dict{String,Vector{Float64}}()
    f_order = String[]
    for i in 1:length(fcid)
        name = f_coil_names === nothing ? @sprintf("FCOIL%03d", i) : f_coil_names[i]
        f_coil[name] = Float64[RF[i], ZF[i], WF[i], HF[i], TURNFC[i]]
        push!(f_order, name)
    end
    machine["FCOIL"] = f_coil
    machine["FCOIL_order"] = f_order
    amp2 = raw["AMP2"]
    probe_angles = Dict{String,Float64}()
    for (i, pname) in enumerate(machine["MPNAM2"])
        probe_angles[pname] = Float64(amp2[i])
    end
    machine["AMP2"] = probe_angles
    return machine, raw
end

"""
    read_kfile(path, machine_dict; e_coil_names=nothing, f_coil_names=nothing)

Read an EFIT k-file of measured constraints. Returns
`(probes, loops, e_coils, f_coils, raw)`, each a name -> `[value, weight]` dict
(plus the raw namelist). `machine_dict` is the result of [`read_mhdin`](@ref).
Mirrors `TokaMaker.util.read_kfile`.
"""
function read_kfile(path::AbstractString, machine_dict::AbstractDict;
                    e_coil_names=nothing, f_coil_names=nothing)
    raw = read_fortran_namelist(path)
    probe_names = machine_dict["MPNAM2"]
    probes = Dict{String,Vector{Float64}}()
    pv = raw["EXPMP2"]; pw = raw["FWTMP2"]
    for i in 1:length(probe_names)
        probes[probe_names[i]] = Float64[pv[i], pw[i]]
    end
    loop_names = machine_dict["LPNAME"]
    loops = Dict{String,Vector{Float64}}()
    lv = raw["COILS"]; lw = raw["FWTSI"]
    for i in 1:length(loop_names)
        loops[loop_names[i]] = Float64[lv[i], lw[i]]
    end
    enames = e_coil_names === nothing ? sort(collect(keys(machine_dict["ECOIL"]))) : e_coil_names
    ev = raw["ECURRT"]; ew = raw["FWTEC"]
    e_coils = Dict{String,Vector{Float64}}()
    for i in 1:length(enames)
        e_coils[enames[i]] = Float64[ev[i], ew[i]]
    end
    fnames = f_coil_names === nothing ? sort(collect(keys(machine_dict["FCOIL"]))) : f_coil_names
    fv = raw["BRSP"]; fw = raw["FWTFC"]
    f_coils = Dict{String,Vector{Float64}}()
    for i in 1:length(fnames)
        f_coils[fnames[i]] = Float64[fv[i], fw[i]]
    end
    return probes, loops, e_coils, f_coils, raw
end

end # module

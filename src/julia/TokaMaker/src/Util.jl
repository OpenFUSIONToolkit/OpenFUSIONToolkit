module Util

using ..LibPath: liboftpy

export create_isoflux, create_power_flux_fun, create_spline_flux_fun,
       eval_green, read_eqdsk, read_ifile

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

end # module

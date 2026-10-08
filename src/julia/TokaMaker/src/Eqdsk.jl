module Eqdsk

# GEQDSK (g-file) raw reader/writer + COCOS conversion — a Julia port of the
# raw-data parts of OpenFUSIONToolkit.TokaMaker.eqdsk (the flux-surface-average
# machinery is not ported). Supports reading/writing standard g-files, COCOS
# index conversion (Sauter & Medvedev sign/2pi factors), Bt/Ip flips, and
# in-memory byte serialization.

using Printf

export GEQDSKEquilibrium, read_geqdsk, cocos_params, cocosify, cocosify!,
       flip_Bt_Ip, flip_Bt_Ip!, save_geqdsk, eqdsk_to_bytes, eqdsk_from_bytes,
       eqdsk_from_raw

# ---------------------------------------------------------------------------
# COCOS parameter table

"""
    cocos_params(idx) -> NamedTuple

Return COCOS sign/exponent parameters (`sigma_Bp`, `sigma_RpZ`, `sigma_rhotp`,
`exp_Bp`) for COCOS index `idx` (1-8 or 11-18). Errors on invalid indices.
"""
function cocos_params(idx::Integer)
    (idx < 1 || idx > 18 || idx == 9 || idx == 10) &&
        error("Invalid COCOS index: $idx")
    exp_Bp = idx < 10 ? 0 : 1
    base = idx < 10 ? idx : idx - 10
    sigma_Bp = base in (1, 2, 5, 6) ? 1 : -1
    sigma_RpZ = base in (1, 2, 7, 8) ? 1 : -1
    sigma_rhotp = base in (1, 3, 5, 7) ? 1 : -1
    return (sigma_Bp=sigma_Bp, sigma_RpZ=sigma_RpZ,
            sigma_rhotp=sigma_rhotp, exp_Bp=exp_Bp)
end

# ---------------------------------------------------------------------------
# Parsing helpers (fixed 16-char fields, 5 per line)

_splitter(text::AbstractString, step::Int=16) =
    [text[step*k+1:step*(k+1)] for k in 0:(div(length(text), step) - 1)]

function _merge(lines::AbstractVector{<:AbstractString})
    isempty(lines) && return ""
    length(lines[1]) > 80 && return replace(join(lines), " " => "")
    return join(lines)
end

_parsevals(lines) = [parse(Float64, strip(c)) for c in _splitter(_merge(lines))]

# ---------------------------------------------------------------------------
# Reader

"""
    _read_geqdsk(filename) -> Dict{String,Any}

Parse a GEQDSK file into a dict with standard uppercase field names.
"""
function _read_geqdsk(filename::AbstractString)
    lines = readlines(filename)
    g = Dict{String,Any}()
    # --- Header ---
    header = lines[1]
    hb = rpad(header, 48)
    g["CASE"] = [hb[8k+1:8(k+1)] for k in 0:5]
    tail = length(header) > 48 ? header[49:end] : ""
    toks = split(tail)
    _, NW, NH = parse.(Int, toks[1:3])
    g["NW"] = NW; g["NH"] = NH
    off = 2
    # --- 20 scalars (4 lines x 5) ---
    s = _parsevals(lines[off:off+3]); off += 4
    g["RDIM"]=s[1]; g["ZDIM"]=s[2]; g["RCENTR"]=s[3]; g["RLEFT"]=s[4]; g["ZMID"]=s[5]
    g["RMAXIS"]=s[6]; g["ZMAXIS"]=s[7]; g["SIMAG"]=s[8]; g["SIBRY"]=s[9]; g["BCENTR"]=s[10]
    g["CURRENT"]=s[11]; g["SIMAG"]=s[12]; g["RMAXIS"]=s[14]; g["ZMAXIS"]=s[16]; g["SIBRY"]=s[18]
    nlNW = ceil(Int, NW / 5)
    # --- 1-D profiles ---
    for name in ("FPOL", "PRES", "FFPRIM", "PPRIME")
        g[name] = _parsevals(lines[off:off+nlNW-1]); off += nlNW
    end
    # --- PSIRZ (packed, NW*NH) stored flat in file order ---
    nlNWNH = ceil(Int, NW * NH / 5)
    g["PSIRZ"] = _parsevals(lines[off:off+nlNWNH-1])[1:NW*NH]; off += nlNWNH
    # --- Safety factor ---
    g["QPSI"] = _parsevals(lines[off:off+nlNW-1]); off += nlNW
    # --- Boundary + limiter ---
    parts = split(strip(lines[off])); off += 1
    NBBBS = parse(Int, parts[1]); LIMITR = parse(Int, parts[2])
    g["NBBBS"] = NBBBS; g["LIMITR"] = LIMITR
    nlNB = ceil(Int, NBBBS * 2 / 5)
    bnd = _parsevals(lines[off:off+nlNB-1]); off += max(nlNB, 1)
    g["RBBBS"] = bnd[1:2:2NBBBS]; g["ZBBBS"] = bnd[2:2:2NBBBS]
    nlLim = ceil(Int, LIMITR * 2 / 5)
    lim = _parsevals(lines[off:off+nlLim-1]); off += nlLim
    g["RLIM"] = lim[1:2:2LIMITR]; g["ZLIM"] = lim[2:2:2LIMITR]
    return g
end

# ---------------------------------------------------------------------------
# Writer

function _write_block(io::IO, vals)
    for (i, v) in enumerate(vals)
        @printf(io, "%16.9E", v)
        i % 5 == 0 && print(io, "\n")
    end
    length(vals) % 5 != 0 && print(io, "\n")
end

function _write_geqdsk(io::IO, g::AbstractDict)
    NW = Int(g["NW"]); NH = Int(g["NH"])
    case_arr = get(g, "CASE", fill(" "^8, 6))
    case_str = join(rpad(c, 8)[1:8] for c in case_arr)
    case_str = rpad(case_str, 48)[1:48]
    @printf(io, "%s%4d%4d%4d\n", case_str, 0, NW, NH)
    scalars = [g["RDIM"], g["ZDIM"], g["RCENTR"], g["RLEFT"], g["ZMID"],
               g["RMAXIS"], g["ZMAXIS"], g["SIMAG"], g["SIBRY"], g["BCENTR"],
               g["CURRENT"], g["SIMAG"], 0.0, g["RMAXIS"], 0.0,
               g["ZMAXIS"], 0.0, g["SIBRY"], 0.0, 0.0]
    _write_block(io, scalars)
    for name in ("FPOL", "PRES", "FFPRIM", "PPRIME")
        _write_block(io, g[name])
    end
    _write_block(io, g["PSIRZ"])
    _write_block(io, g["QPSI"])
    nbbbs = Int(get(g, "NBBBS", length(get(g, "RBBBS", Float64[]))))
    nlim = Int(get(g, "LIMITR", length(get(g, "RLIM", Float64[]))))
    @printf(io, " %5d %5d\n", nbbbs, nlim)
    if nbbbs > 0
        bnd = Vector{Float64}(undef, 2nbbbs)
        bnd[1:2:end] = g["RBBBS"][1:nbbbs]; bnd[2:2:end] = g["ZBBBS"][1:nbbbs]
        _write_block(io, bnd)
    end
    if nlim > 0
        lim = Vector{Float64}(undef, 2nlim)
        lim[1:2:end] = g["RLIM"][1:nlim]; lim[2:2:end] = g["ZLIM"][1:nlim]
        _write_block(io, lim)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Equilibrium object

mutable struct GEQDSKEquilibrium
    raw::Dict{String,Any}
    cocos::Int
end

"""
    read_geqdsk(filename; cocos=1) -> GEQDSKEquilibrium

Read a GEQDSK file. `cocos` is validated (1-8 or 11-18).
"""
function read_geqdsk(filename::AbstractString; cocos::Integer=1)
    cocos_params(cocos)  # validates
    return GEQDSKEquilibrium(_read_geqdsk(filename), Int(cocos))
end

"""
    eqdsk_from_raw(raw; cocos=1) -> GEQDSKEquilibrium

Construct directly from a raw g-file dict (deep-copies arrays).
"""
function eqdsk_from_raw(raw::AbstractDict; cocos::Integer=1)
    cocos_params(cocos)
    cp = Dict{String,Any}()
    for (k, v) in raw
        cp[k] = v isa AbstractArray ? copy(v) : v
    end
    return GEQDSKEquilibrium(cp, Int(cocos))
end

_copy_for_mutation(eq::GEQDSKEquilibrium) = eqdsk_from_raw(eq.raw; cocos=eq.cocos)

# ---------------------------------------------------------------------------
# COCOS conversion

"""
    cocosify!(eq, cocos_out) -> eq

Convert the raw g-file data in place from the current COCOS to `cocos_out`,
applying the Sauter & Medvedev (CPC 184, 2013) sign / 2π factors.
"""
function cocosify!(eq::GEQDSKEquilibrium, cocos_out::Integer)
    cc_in = cocos_params(eq.cocos)
    cc_out = cocos_params(cocos_out)
    sBp = cc_out.sigma_Bp * cc_in.sigma_Bp
    sRpZ = cc_out.sigma_RpZ * cc_in.sigma_RpZ
    srhotp = cc_out.sigma_rhotp * cc_in.sigma_rhotp
    exp_eff = cc_out.exp_Bp - cc_in.exp_Bp
    twopi_exp = (2π)^exp_eff
    psi_fac = sRpZ * sBp * twopi_exp
    dpsi_fac = sRpZ * sBp / twopi_exp
    bt_fac = sRpZ
    ip_fac = sRpZ
    q_fac = srhotp
    g = eq.raw
    g["SIMAG"] *= psi_fac
    g["SIBRY"] *= psi_fac
    g["PSIRZ"] = g["PSIRZ"] .* psi_fac
    g["PPRIME"] = g["PPRIME"] .* dpsi_fac
    g["FFPRIM"] = g["FFPRIM"] .* dpsi_fac
    g["FPOL"] = g["FPOL"] .* bt_fac
    g["BCENTR"] *= bt_fac
    g["CURRENT"] *= ip_fac
    g["QPSI"] = g["QPSI"] .* q_fac
    eq.cocos = Int(cocos_out)
    return eq
end

"""
    cocosify(eq, cocos_out) -> GEQDSKEquilibrium

Non-mutating COCOS conversion: returns a converted copy, leaving `eq` unchanged.
"""
cocosify(eq::GEQDSKEquilibrium, cocos_out::Integer) =
    cocosify!(_copy_for_mutation(eq), cocos_out)

"""
    flip_Bt_Ip!(eq) -> eq

Reverse the signs of Bt and Ip (negates BCENTR, FPOL, CURRENT, PSIRZ, SIMAG,
SIBRY, PPRIME, FFPRIM), keeping the COCOS index unchanged.
"""
function flip_Bt_Ip!(eq::GEQDSKEquilibrium)
    g = eq.raw
    for k in ("BCENTR", "CURRENT", "SIMAG", "SIBRY")
        g[k] *= -1
    end
    for k in ("FPOL", "PSIRZ", "PPRIME", "FFPRIM")
        g[k] = g[k] .* -1
    end
    return eq
end

flip_Bt_Ip(eq::GEQDSKEquilibrium) = flip_Bt_Ip!(_copy_for_mutation(eq))

# ---------------------------------------------------------------------------
# Serialization

"""
    save_geqdsk(eq, filename)

Write the (possibly modified) g-file data to `filename`.
"""
function save_geqdsk(eq::GEQDSKEquilibrium, filename::AbstractString)
    open(filename, "w") do io
        _write_geqdsk(io, eq.raw)
    end
end

"""
    eqdsk_to_bytes(eq) -> Vector{UInt8}

Serialize to in-memory bytes (round-trips with [`eqdsk_from_bytes`](@ref)).
"""
function eqdsk_to_bytes(eq::GEQDSKEquilibrium)
    io = IOBuffer()
    _write_geqdsk(io, eq.raw)
    return take!(io)
end

"""
    eqdsk_from_bytes(bytes; cocos=1) -> GEQDSKEquilibrium

Construct from raw g-file bytes (via a temporary file).
"""
function eqdsk_from_bytes(bytes::AbstractVector{UInt8}; cocos::Integer=1)
    cocos_params(cocos)
    path, io = mktemp()
    try
        write(io, bytes); close(io)
        return read_geqdsk(path; cocos=cocos)
    finally
        rm(path; force=true)
    end
end

end # module

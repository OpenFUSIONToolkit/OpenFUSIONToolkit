module Reconstruction

# Port of `reconstruction.py`. Constraint file format mirrors Python exactly:
#
#   <ncons>
#
#   <type-id>
#    <floats — multi-line per constraint type>
#
#   <type-id>
#    ...
#
# Crucially, the per-line "error" field is serialized as 1/err (inverse), not
# err itself. The Fortran reader expects this (it's the weight). See
# `reconstruction.py:Mirnov_con.write` etc.

using Printf: @printf, @sprintf
using ..LibPath: liboftpy
using ..CInterface
using ..SettingsModule: TokamakerReconSettings, tokamaker_recon_default_settings
using ..CoreModule: Tokamaker
using ..OFTEnvModule: OFTEnv, _write_oft_in

export MirnovCon, IpCon, FluxLoopCon, DFluxCon, PressCon, QCon, SaddleCon,
       ReconConstraints, run_reconstruction!, write_constraints_file

# Constraint type IDs (match Python module constants)
const _MIRNOV_ID    = 1
const _IP_ID        = 2
const _FLUX_LOOP_ID = 7
const _DFLUX_ID     = 8
const _PRESS_ID     = 9
const _Q_ID         = 10
const _SADDLE_ID    = 11

abstract type ReconConstraint end

struct MirnovCon <: ReconConstraint
    loc::NTuple{2,Float64}
    phi::Float64
    norm::NTuple{3,Float64}   # (R, phi, Z) components in R-phi-Z basis
    val::Float64
    err::Float64
end
MirnovCon(loc::AbstractVector, norm::AbstractVector, val::Real, err::Real; phi::Real=0.0) =
    MirnovCon((Float64(loc[1]), Float64(loc[2])), Float64(phi),
              (Float64(norm[1]), Float64(norm[2]), Float64(norm[3])),
              Float64(val), Float64(err))
function _write_con(io::IO, c::MirnovCon)
    println(io, _MIRNOV_ID)
    @printf io " %E %E %E\n" c.loc[1] c.loc[2] c.phi
    @printf io " %E %E %E\n" c.norm[1] c.norm[2] c.norm[3]
    @printf io " %E %E\n\n" c.val (1.0 / c.err)
end

struct IpCon <: ReconConstraint
    val::Float64
    err::Float64
end
function _write_con(io::IO, c::IpCon)
    println(io, _IP_ID)
    @printf io " %E %E\n\n" c.val (1.0 / c.err)
end

struct FluxLoopCon <: ReconConstraint
    loc::NTuple{2,Float64}
    val::Float64
    err::Float64
end
FluxLoopCon(loc::AbstractVector, val::Real, err::Real) =
    FluxLoopCon((Float64(loc[1]), Float64(loc[2])), Float64(val), Float64(err))
function _write_con(io::IO, c::FluxLoopCon)
    println(io, _FLUX_LOOP_ID)
    @printf io " %E %E\n" c.loc[1] c.loc[2]
    @printf io " %E %E\n\n" c.val (1.0 / c.err)
end

struct DFluxCon <: ReconConstraint
    val::Float64
    err::Float64
end
function _write_con(io::IO, c::DFluxCon)
    println(io, _DFLUX_ID)
    @printf io " %E %E\n\n" c.val (1.0 / c.err)
end

struct PressCon <: ReconConstraint
    loc::NTuple{2,Float64}
    val::Float64
    err::Float64
end
PressCon(loc::AbstractVector, val::Real, err::Real) =
    PressCon((Float64(loc[1]), Float64(loc[2])), Float64(val), Float64(err))
function _write_con(io::IO, c::PressCon)
    println(io, _PRESS_ID)
    @printf io " %E %E\n" c.loc[1] c.loc[2]
    @printf io " %E %E\n\n" c.val (1.0 / c.err)
end

struct QCon <: ReconConstraint
    type::Int
    loc::Float64
    val::Float64
    err::Float64
end
function _write_con(io::IO, c::QCon)
    println(io, _Q_ID)
    @printf io " %d %E\n" c.type c.loc
    @printf io " %E %E\n\n" c.val (1.0 / c.err)
end

struct SaddleCon <: ReconConstraint
    pt1::NTuple{2,Float64}
    pt2::NTuple{2,Float64}
    width::Float64
    val::Float64
    err::Float64
end
SaddleCon(pt1::AbstractVector, pt2::AbstractVector, width::Real, val::Real, err::Real) =
    SaddleCon((Float64(pt1[1]), Float64(pt1[2])),
              (Float64(pt2[1]), Float64(pt2[2])),
              Float64(width), Float64(val), Float64(err))
function _write_con(io::IO, c::SaddleCon)
    println(io, _SADDLE_ID)
    @printf io " %E %E\n" c.pt1[1] c.pt1[2]
    @printf io " %E %E\n" c.pt2[1] c.pt2[2]
    @printf io " %E %E %E\n\n" c.width c.val (1.0 / c.err)
end

# ---------------------------------------------------------------------------

mutable struct ReconConstraints
    items::Vector{ReconConstraint}
    ReconConstraints() = new(ReconConstraint[])
end

Base.push!(rc::ReconConstraints, c::ReconConstraint) = (push!(rc.items, c); rc)
Base.length(rc::ReconConstraints) = length(rc.items)

function write_constraints_file(filename::AbstractString, rc::ReconConstraints)
    open(filename, "w") do io
        println(io, length(rc.items), "\n")
        for c in rc.items
            _write_con(io, c)
        end
    end
    return filename
end

# ---------------------------------------------------------------------------

"""
    run_reconstruction!(tm, constraints; settings, vacuum=false,
                        infile="fit.in", outfile="fit.out", linearized_fit=false,
                        maxits=100, eps=1e-3, ftol=1e-3, xtol=1e-3, gtol=1e-3)

Write the constraint file, update the OFT input file's `gs_fit_options`
group with the requested solver controls, then call `tokamaker_recon_run`.
Returns the integer error flag from Fortran.
"""
function run_reconstruction!(t::Tokamaker, constraints::ReconConstraints;
                              settings::TokamakerReconSettings=tokamaker_recon_default_settings(),
                              vacuum::Bool=false,
                              infile::AbstractString="fit.in",
                              outfile::AbstractString="fit.out",
                              linearized_fit::Bool=false,
                              maxits::Integer=100,
                              eps::Real=1e-3, ftol::Real=1e-3,
                              xtol::Real=1e-3, gtol::Real=1e-3)
    write_constraints_file(infile, constraints)
    env = t.env
    env.oft_in_groups["gs_fit_options"] = Dict{String,String}(
        "linearized_fit" => linearized_fit ? "T" : "F",
        "maxfev" => string(maxits),
        "epsfcn" => @sprintf("%.5E", eps),
        "ftol"   => @sprintf("%.5E", ftol),
        "xtol"   => @sprintf("%.5E", xtol),
        "gtol"   => @sprintf("%.5E", gtol),
    )
    _write_oft_in(env.ifile, env.oft_in_groups)

    in_buf = Vector{UInt8}(codeunits(infile)); push!(in_buf, 0x00)
    out_buf = Vector{UInt8}(codeunits(outfile)); push!(out_buf, 0x00)
    err_flag = Ref{Int32}(0)
    GC.@preserve in_buf out_buf begin
        settings.infile = pointer(in_buf)
        settings.outfile = pointer(out_buf)
        ccall((:tokamaker_recon_run, liboftpy[]), Cvoid,
              (Ptr{Cvoid}, Cuchar, Ref{TokamakerReconSettings}, Ptr{Int32}),
              t.tmaker_ptr, UInt8(vacuum), settings, err_flag)
    end
    return Int(err_flag[])
end

end # module

module Reconstruction

# Port of `reconstruction.py`. Each constraint type is a small Julia struct
# with the same field semantics as the Python class. Constraint files are
# plain text; the format is:
#   <type-id> <number-of-floats> <floats...>
# matching Python's serializer.

using ..LibPath: liboftpy
using ..CInterface
using ..SettingsModule: TokamakerReconSettings, tokamaker_recon_default_settings
using ..CoreModule: Tokamaker

export MirnovCon, IpCon, FluxLoopCon, DFluxCon, PressCon, QCon, SaddleCon,
       ReconConstraints, run_reconstruction!, write_constraints_file

# Constraint types -----------------------------------------------------------

abstract type ReconConstraint end

# Type IDs match `id` field in Python classes (see reconstruction.py)
struct MirnovCon <: ReconConstraint
    loc::NTuple{2,Float64}
    phi::Float64    # toroidal angle (legacy, normally 0)
    norm::NTuple{2,Float64}  # local probe orientation
    val::Float64
    err::Float64
end
_typeid(::MirnovCon) = 1
_serialize(c::MirnovCon) = (1, [c.loc[1], c.loc[2], c.phi, c.norm[1], c.norm[2], c.val, c.err])

struct IpCon <: ReconConstraint
    val::Float64
    err::Float64
end
_typeid(::IpCon) = 2
_serialize(c::IpCon) = (2, [c.val, c.err])

struct FluxLoopCon <: ReconConstraint
    loc::NTuple{2,Float64}
    val::Float64
    err::Float64
end
_typeid(::FluxLoopCon) = 7
_serialize(c::FluxLoopCon) = (7, [c.loc[1], c.loc[2], c.val, c.err])

struct DFluxCon <: ReconConstraint
    val::Float64
    err::Float64
end
_typeid(::DFluxCon) = 8
_serialize(c::DFluxCon) = (8, [c.val, c.err])

struct PressCon <: ReconConstraint
    loc::NTuple{2,Float64}
    val::Float64
    err::Float64
end
_typeid(::PressCon) = 9
_serialize(c::PressCon) = (9, [c.loc[1], c.loc[2], c.val, c.err])

struct QCon <: ReconConstraint
    type::Int  # 0 -> q at psi value, 1 -> q at minor radius
    loc::Float64
    val::Float64
    err::Float64
end
_typeid(::QCon) = 10
_serialize(c::QCon) = (10, [Float64(c.type), c.loc, c.val, c.err])

struct SaddleCon <: ReconConstraint
    pt1::NTuple{2,Float64}
    pt2::NTuple{2,Float64}
    width::Float64
    val::Float64
    err::Float64
end
_typeid(::SaddleCon) = 11
_serialize(c::SaddleCon) = (11, [c.pt1[1], c.pt1[2], c.pt2[1], c.pt2[2], c.width, c.val, c.err])

# Container -----------------------------------------------------------------

mutable struct ReconConstraints
    items::Vector{ReconConstraint}
    ReconConstraints() = new(ReconConstraint[])
end

Base.push!(rc::ReconConstraints, c::ReconConstraint) = push!(rc.items, c)
Base.length(rc::ReconConstraints) = length(rc.items)

function write_constraints_file(filename::AbstractString, rc::ReconConstraints)
    open(filename, "w") do io
        println(io, length(rc.items))
        for c in rc.items
            tid, vals = _serialize(c)
            print(io, tid, " ", length(vals))
            for v in vals
                print(io, " ", v)
            end
            println(io)
        end
    end
    return filename
end

# Run ------------------------------------------------------------------------

"""
    run_reconstruction!(tm, constraints; settings=default, vacuum=false)

Run the reconstruction solver. Constraints are written to a temporary text
file (see `tokamaker_recon_settings_struct.infile`) then `tokamaker_recon_run`
is invoked. Returns the `error_flag` from the Fortran side.
"""
function run_reconstruction!(t::Tokamaker, constraints::ReconConstraints;
                              settings::TokamakerReconSettings=tokamaker_recon_default_settings(),
                              vacuum::Bool=false,
                              infile::AbstractString="fit.in",
                              outfile::AbstractString="fit.out")
    write_constraints_file(infile, constraints)
    in_buf = Vector{UInt8}(codeunits(infile))
    push!(in_buf, 0x00)
    out_buf = Vector{UInt8}(codeunits(outfile))
    push!(out_buf, 0x00)
    GC.@preserve in_buf out_buf begin
        settings.infile = pointer(in_buf)
        settings.outfile = pointer(out_buf)
        err_flag = Ref{Int32}(0)
        ccall((:tokamaker_recon_run, liboftpy[]), Cvoid,
              (Ptr{Cvoid}, Cuchar, Ref{TokamakerReconSettings}, Ptr{Int32}),
              t.tmaker_ptr, UInt8(vacuum), settings, err_flag)
        return Int(err_flag[])
    end
end

end # module

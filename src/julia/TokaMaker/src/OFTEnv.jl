module OFTEnvModule

using ..LibPath: liboftpy, OFT_PATH_SLEN, load_libraries!

export OFTEnv, oftpy_set_debug, oftpy_set_nthreads

mutable struct OFTEnv
    nthreads::Int32
    quiet::Bool
    ifile::String
    ifile_buf::Vector{UInt8}     # kept alive for the lifetime of OFT
    debug_level::Int
    oft_in_groups::Dict{String,Dict{String,String}}
    slens::Vector{Int32}
    initialized::Bool
end

function _write_oft_in(path::AbstractString, groups::AbstractDict)
    open(path, "w") do io
        for (name, options) in groups
            println(io, "&$name")
            for (k, v) in options
                println(io, "  $k=$v")
            end
            println(io, "/\n")
        end
    end
end

function _abort()
    flush(stdout)
    flush(stderr)
    ccall(:exit, Cvoid, (Cint,), -1)
    return nothing
end

const _ABORT_CFUN = Ref{Ptr{Cvoid}}(C_NULL)

function _path_cstr(s::AbstractString)
    n = length(s)
    n > OFT_PATH_SLEN - 1 && error("Path string longer than $(OFT_PATH_SLEN-1) bytes: $s")
    return Cstring(Base.unsafe_convert(Ptr{UInt8}, Base.cconvert(Cstring, s)))
end

const _GLOBAL_ENV = Ref{Union{OFTEnv,Nothing}}(nothing)

"""
    OFTEnv(; nthreads=-1, quiet=false, ifile="")

Initialize the OFT runtime. Equivalent to `OFT_env(...)` in Python.
Only one runtime per process is permitted; subsequent constructions return the existing env.
"""
function OFTEnv(; nthreads::Integer=-1, quiet::Bool=false,
                  ifile::AbstractString="", debug_level::Integer=0,
                  use_abort_callback::Bool=true)
    if _GLOBAL_ENV[] !== nothing
        return _GLOBAL_ENV[]
    end
    load_libraries!()

    groups = Dict{String,Dict{String,String}}(
        "runtime_options" => Dict("debug" => string(debug_level)),
        "mesh_options"    => Dict("meshname" => "'none'"),
    )
    actual_ifile = isempty(ifile) ? "oftpyin-$(getpid())" : String(ifile)
    _write_oft_in(actual_ifile, groups)

    # Path string must remain valid for the lifetime of the OFT runtime.
    ifile_buf = Vector{UInt8}(codeunits(actual_ifile))
    push!(ifile_buf, 0x00)

    slens = zeros(Int32, 4)
    if _ABORT_CFUN[] == C_NULL
        _ABORT_CFUN[] = @cfunction(_abort, Cvoid, ())
    end
    cb = use_abort_callback ? _ABORT_CFUN[] : C_NULL

    GC.@preserve ifile_buf begin
        ccall((:oftpy_init, liboftpy[]), Cvoid,
              (Cint, Cuchar, Ptr{UInt8}, Ptr{Int32}, Ptr{Cvoid}),
              Int32(nthreads), UInt8(quiet), pointer(ifile_buf), slens, cb)
    end

    env = OFTEnv(Int32(nthreads), quiet, actual_ifile, ifile_buf,
                 Int(debug_level), groups, slens, true)
    _GLOBAL_ENV[] = env
    return env
end

oftpy_set_debug(level::Integer) =
    ccall((:oftpy_set_debug, liboftpy[]), Cvoid, (Cint,), Int32(level))

oftpy_set_nthreads(nthreads::Integer) =
    ccall((:oftpy_set_nthreads, liboftpy[]), Cvoid, (Cint,), Int32(nthreads))

end # module

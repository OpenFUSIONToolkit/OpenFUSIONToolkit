module LibPath

using Libdl

export liboftpy, liboft_triangle, OFT_PATH_SLEN, OFT_ERROR_SLEN, load_libraries!

const OFT_PATH_SLEN = 200
const OFT_ERROR_SLEN = 200

function _lib_suffix()
    if Sys.islinux()
        return ".so"
    elseif Sys.isapple()
        return ".dylib"
    else
        error("Unsupported platform for liboftpy: $(Sys.KERNEL)")
    end
end

function _candidate_paths(libname::String)
    suffix = _lib_suffix()
    fname = "lib" * libname * suffix
    paths = String[]
    if haskey(ENV, "OFT_LIBRARY_PATH")
        push!(paths, joinpath(ENV["OFT_LIBRARY_PATH"], fname))
    end
    pkgsrc = @__DIR__
    push!(paths, joinpath(pkgsrc, "..", "..", "..", "python", "OpenFUSIONToolkit", fname))
    # Canonical CI build layout: <repo>/builds/install_{release,debug}/bin
    push!(paths, joinpath(pkgsrc, "..", "..", "..", "..", "builds", "install_release", "bin", fname))
    push!(paths, joinpath(pkgsrc, "..", "..", "..", "..", "builds", "install_debug", "bin", fname))
    # In-source build (older convention)
    push!(paths, joinpath(pkgsrc, "..", "..", "..", "install_release", "bin", fname))
    push!(paths, joinpath(pkgsrc, "..", "..", "..", "install_debug", "bin", fname))
    push!(paths, joinpath(pkgsrc, "..", "..", "..", "..", "bin", fname))
    return paths
end

function _locate(libname::String)
    for p in _candidate_paths(libname)
        if isfile(p)
            return abspath(p)
        end
    end
    sysfound = Libdl.find_library([libname])
    if !isempty(sysfound)
        return sysfound
    end
    error("""Unable to locate $libname shared library. Tried:
$(join(_candidate_paths(libname), "\n"))
Set ENV["OFT_LIBRARY_PATH"] or build OFT with -DOFT_BUILD_PYTHON=ON.""")
end

const liboftpy = Ref{String}("")
const liboft_triangle = Ref{String}("")
const _libhandle = Ref{Ptr{Cvoid}}(C_NULL)

"Force-load the OFT shared libraries. Called automatically by `OFTEnv`, but
exposed so user code can opt in earlier (e.g. for testing)."
function load_libraries!()
    if isempty(liboftpy[])
        liboftpy[] = _locate("oftpy")
    end
    try
        if isempty(liboft_triangle[])
            liboft_triangle[] = _locate("oft_triangle")
        end
    catch
        liboft_triangle[] = ""  # optional
    end
    if _libhandle[] == C_NULL
        _libhandle[] = Libdl.dlopen(liboftpy[], Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
    end
    return liboftpy[]
end

function __init__()
    # Best-effort: try to locate the library at module init, but don't fail
    # the import if it's missing — users may want to load TokaMaker for static
    # introspection. The first OFTEnv() call will load the library or raise.
    try
        liboftpy[] = _locate("oftpy")
        try
            liboft_triangle[] = _locate("oft_triangle")
        catch
            liboft_triangle[] = ""
        end
    catch
        liboftpy[] = ""
        liboft_triangle[] = ""
    end
end

end # module

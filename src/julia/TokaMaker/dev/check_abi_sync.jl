#!/usr/bin/env julia
# Verify the Julia CInterface.jl wrappers cover every C symbol exposed by
# the Python ctypes layer at src/python/OpenFUSIONToolkit/TokaMaker/_interface.py
# and src/python/OpenFUSIONToolkit/_interface.py. Exits 1 on drift.
#
# Run from anywhere. Usage:
#   julia --project=. dev/check_abi_sync.jl

const PKG_ROOT = abspath(joinpath(@__DIR__, ".."))
const REPO_ROOT = abspath(joinpath(PKG_ROOT, "..", "..", ".."))

py_files = [
    joinpath(REPO_ROOT, "src", "python", "OpenFUSIONToolkit", "_interface.py"),
    joinpath(REPO_ROOT, "src", "python", "OpenFUSIONToolkit", "TokaMaker", "_interface.py"),
    joinpath(REPO_ROOT, "src", "python", "OpenFUSIONToolkit", "TokaMaker", "util.py"),
    joinpath(REPO_ROOT, "src", "python", "OpenFUSIONToolkit", "TokaMaker", "reconstruction.py"),
]
jl_file = joinpath(PKG_ROOT, "src", "CInterface.jl")
jl_extras = [
    joinpath(PKG_ROOT, "src", "Util.jl"),
    joinpath(PKG_ROOT, "src", "OFTEnv.jl"),
    joinpath(PKG_ROOT, "src", "Reconstruction.jl"),
]

# Pattern: oftpy_lib.<name> — captures the Fortran C symbol name on the Python side.
const PY_RE = r"oftpy_lib\.([A-Za-z_][A-Za-z0-9_]*)"
const JL_RE = r":([A-Za-z_][A-Za-z0-9_]*),\s*lib"

function collect_symbols(file::String, regex::Regex)
    isfile(file) || error("Missing $file")
    s = read(file, String)
    syms = Set{String}()
    for m in eachmatch(regex, s)
        push!(syms, m[1])
    end
    return syms
end

py_syms = Set{String}()
for f in py_files
    union!(py_syms, collect_symbols(f, PY_RE))
end

jl_syms = Set{String}()
union!(jl_syms, collect_symbols(jl_file, JL_RE))
for f in jl_extras
    isfile(f) && union!(jl_syms, collect_symbols(f, JL_RE))
end

# Symbols Python exposes but Julia doesn't bind (yet)
missing_in_jl = setdiff(py_syms, jl_syms)
# Symbols Julia binds but Python doesn't expose (likely typo or dead code)
missing_in_py = setdiff(jl_syms, py_syms)

ok = true
if !isempty(missing_in_jl)
    println("Symbols in Python _interface.py but missing in Julia CInterface.jl:")
    for s in sort(collect(missing_in_jl))
        println("  - $s")
    end
    ok = false
end
if !isempty(missing_in_py)
    println("Symbols in Julia CInterface.jl but missing in Python _interface.py:")
    for s in sort(collect(missing_in_py))
        println("  - $s")
    end
    ok = false
end

if ok
    println("ABI sync OK: $(length(py_syms)) Python symbols, $(length(jl_syms)) Julia symbols, all aligned.")
    exit(0)
else
    println()
    println("Total Python: $(length(py_syms)), total Julia: $(length(jl_syms))")
    exit(1)
end

# Native session-file I/O round-trip (Milestone B) + cross-language parity.
#
# 1. Pure-Julia: save_tokamaker -> load (replace_eq! source_file) reproduces psi.
# 2. Julia -> Python: a session written by Julia loads in the Python binding.
# 3. Python -> Julia: a session written by Python loads in Julia.
#
# Both bindings wrap the same Fortran save/load routines, so a faithful file is
# binding-agnostic; these tests verify the Julia ccall wiring produces/consumes
# files the Python ctypes wiring agrees with. The cross-language tests are
# skipped gracefully when a usable Python + OpenFUSIONToolkit install is absent.

using Test
using TokaMaker

const SIO_MESH = joinpath(@__DIR__, "data", "spheromak_h1.h5")
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const PY_PKG = joinpath(REPO_ROOT, "src", "python")
const XLANG_PY = joinpath(@__DIR__, "xlang_session.py")

# Self-contained spheromak solve (mirrors test_spheromak.jl/run_sph_case and the
# Python xlang_session helper) returning the live solver.
function solve_spheromak_gs(mesh_file::AbstractString, order::Integer)
    pts, lc, reg, _, _ = load_gs_mesh(mesh_file)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs)
    gs.settings.free_boundary = false
    setup!(gs; order=order)
    set_p_scale!(gs, 0.0)
    set_profiles!(gs; ffp_prof=Dict("type" => "linterp", "x" => [0.0, 1.0], "y" => [1.0, 0.0]),
                  pp_prof=Dict("type" => "flat"))
    gs.settings.nl_tol = 1e-12
    gs.settings.maxits = 100
    gs.settings.urf = 0.0
    update_settings!(gs)
    init_psi!(gs)
    solve!(gs)
    return gs
end

# Build + set up the mesh and create an equilibrium (via init_psi!) ready to be
# overwritten by a loaded session file. No solve.
function setup_spheromak_gs(mesh_file::AbstractString, order::Integer)
    pts, lc, reg, _, _ = load_gs_mesh(mesh_file)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs)
    gs.settings.free_boundary = false
    setup!(gs; order=order)
    init_psi!(gs)
    return gs
end

_read_psi(path) = parse.(Float64, readlines(path))

# Build a `uv run --script` command for the cross-language helper. `uv`
# provisions the helper's PEP 723 dependencies (numpy/scipy/h5py) in an
# ephemeral env; PYTHONPATH points that env at the in-repo OpenFUSIONToolkit
# package. Returns `nothing` when `uv` is not on PATH.
function _uv_cmd(args...)
    uv = Sys.which("uv")
    uv === nothing && return nothing
    return addenv(Cmd(String[uv, "run", "--script", XLANG_PY, args...]),
                  "PYTHONPATH" => PY_PKG)
end

# Detect a usable uv + Python + OpenFUSIONToolkit binding sharing the same lib.
function _python_available()
    _uv_cmd() === nothing && return false
    try
        return success(_uv_cmd("probe"))
    catch
        return false
    end
end

@testset "Session I/O" begin
    if !isfile(SIO_MESH)
        @warn "spheromak mesh missing at $SIO_MESH"
        @test_skip true
    else
        mktempdir() do dir
            # --- 1. Pure-Julia round-trip ---------------------------------
            gs = solve_spheromak_gs(SIO_MESH, 2)
            psi0 = get_psi(gs; normalized=false)
            session_J = joinpath(dir, "session_julia.h5")
            save_tokamaker(gs.equilibrium, session_J)
            @test isfile(session_J)

            # Capability probe: liboftpy and HDF5.jl must share one libhdf5,
            # otherwise liboftpy's HDF5 writes silently no-op (empty file) and
            # all session I/O fails. See the "HDF5 setup" section in README.md.
            if !isfile(session_J) || filesize(session_J) < 1024
                @warn """Session I/O unavailable: liboftpy wrote an empty HDF5 file \
($(isfile(session_J) ? filesize(session_J) : 0) bytes). HDF5.jl is using a \
different libhdf5 than liboftpy — point HDF5.jl at liboftpy's libhdf5 via a \
LocalPreferences.toml [HDF5] block (see README "HDF5 setup"), then re-run. \
Skipping session-I/O tests."""
                @test_skip true
                return
            end

            @testset "Julia save -> Julia load (same solver)" begin
                replace_eq!(gs; source_file=session_J)
                psi_rt = get_psi(gs; normalized=false)
                @test psi_rt ≈ psi0 rtol = 1e-12
            end

            @testset "Julia save -> Julia load (fresh solver)" begin
                gs2 = setup_spheromak_gs(SIO_MESH, 2)
                replace_eq!(gs2; source_file=session_J)
                psi_fresh = get_psi(gs2; normalized=false)
                @test length(psi_fresh) == length(psi0)
                @test psi_fresh ≈ psi0 rtol = 1e-10
            end

            # --- Cross-language (gated on a usable Python install) --------
            if !_python_available()
                @info "Skipping cross-language session tests: no usable uv + OpenFUSIONToolkit on PYTHONPATH=$PY_PKG"
                @test_skip true
            else
                runpy(args...) = run(_uv_cmd(args...))

                @testset "Julia save -> Python load" begin
                    psi_PJ = joinpath(dir, "psi_py_from_julia.txt")
                    runpy("load", SIO_MESH, session_J, psi_PJ)
                    psi_py = _read_psi(psi_PJ)
                    @test length(psi_py) == length(psi0)
                    @test psi_py ≈ psi0 rtol = 1e-8
                end

                @testset "Python save -> Julia load" begin
                    session_P = joinpath(dir, "session_python.h5")
                    psi_P = joinpath(dir, "psi_python.txt")
                    runpy("save", SIO_MESH, session_P, psi_P)
                    @test isfile(session_P)
                    psi_python = _read_psi(psi_P)
                    gs3 = setup_spheromak_gs(SIO_MESH, 2)
                    replace_eq!(gs3; source_file=session_P)
                    psi_JL = get_psi(gs3; normalized=false)
                    @test length(psi_JL) == length(psi_python)
                    @test psi_JL ≈ psi_python rtol = 1e-8
                end
            end
        end
    end
end

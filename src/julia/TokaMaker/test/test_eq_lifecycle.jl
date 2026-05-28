# Equilibrium copy/replace lifecycle (Milestone B).
# Verifies copy_eq produces an independent snapshot that carries the
# targets/constraints that produced it, and that replace_eq! swaps the active
# equilibrium back in. Mirrors the intent of test_ITER_concurrent in
# src/tests/physics/test_TokaMaker.py (concurrent equilibria for one solver).

using Test
using LinearAlgebra: norm
using TokaMaker

const EQLC_DATA_DIR = joinpath(@__DIR__, "data")

# Build and solve a small fixed-boundary Solov'ev equilibrium (same setup as
# test_solovev.jl) so the lifecycle test is cheap.
function _solo_solved()
    R, a_, b_, c0 = 1.0, 1.2, -1.0, 1.1
    pts, tris, regs, _, _ = load_gs_mesh(joinpath(EQLC_DATA_DIR, "solo_h1.h5"))
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=tris, reg=regs)
    setup_regions!(gs)
    gs.settings.free_boundary = false
    gs.settings.pm = false
    setup!(gs; order=2, F0=1.0, full_domain=true)
    set_p_scale!(gs, a_)
    set_ffp_scale!(gs, b_ * R^2 * 2.0)
    set_profiles!(gs; ffp_prof=Dict("type" => "flat"), pp_prof=Dict("type" => "flat"))
    init_psi!(gs)
    ζ = @. (gs.r[:, 1]^2 - R^2) / (2R)
    psi_solo = @. (b_ + c0) * R^2 * gs.r[:, 2]^2 / 2 + c0 * R * ζ * gs.r[:, 2]^2 +
                  (a_ - c0) * R^2 * ζ^2 / 2
    set_psi!(gs, -psi_solo)
    gs.settings.nl_tol = 1e-12
    gs.settings.maxits = 40
    update_settings!(gs)
    solve!(gs)
    return gs, psi_solo
end

@testset "Equilibrium copy/replace lifecycle" begin
    gs, psi_solo = _solo_solved()

    # Targets + constraints so we can check the snapshot carries them.
    set_targets!(gs; Ip=1.0e3)
    set_psi_constraints!(gs, [1.2 0.0; 0.8 0.0], [0.0, 0.0])

    psi_orig = get_psi(gs; normalized=false)

    @testset "copy_eq snapshot" begin
        snap = copy_eq(gs)
        psi_snap = TokaMaker.EquilibriumModule.get_psi(snap; normalized=false)
        @test psi_snap ≈ psi_orig rtol = 1e-12
        @test get(snap.targets, "Ip", -1) == 1.0e3
        @test snap.psi_constraints !== nothing
        @test snap.psi_constraints.locations == [1.2 0.0; 0.8 0.0]
    end

    @testset "skip_targets / skip_constraints" begin
        snap2 = copy_eq(gs; skip_targets=true, skip_constraints=true)
        @test isempty(snap2.targets)
        @test snap2.psi_constraints === nothing
        @test snap2.isoflux_constraints === nothing
    end

    @testset "snapshot independence" begin
        snap = copy_eq(gs)
        psi_before = TokaMaker.EquilibriumModule.get_psi(snap; normalized=false)
        set_psi!(gs, -psi_solo .* 0.5)  # perturb the active equilibrium
        psi_after = TokaMaker.EquilibriumModule.get_psi(snap; normalized=false)
        @test psi_after ≈ psi_before rtol = 1e-12                      # snapshot frozen
        @test norm(get_psi(gs; normalized=false) .- psi_after) > 1e-6  # active changed
    end

    @testset "replace_eq! from source_eq" begin
        snap = copy_eq(gs)
        psi_snap = TokaMaker.EquilibriumModule.get_psi(snap; normalized=false)
        set_psi!(gs, -psi_solo .* 0.25)  # diverge active
        replace_eq!(gs; source_eq=snap)
        @test get_psi(gs; normalized=false) ≈ psi_snap rtol = 1e-10
        @test get(gs.targets, "Ip", -1) == 1.0e3   # live targets synced
        @test gs.psi_constraints !== nothing        # live constraints synced
    end

    @testset "replace_eq! argument validation" begin
        @test_throws ErrorException replace_eq!(gs)  # neither source given
    end
end

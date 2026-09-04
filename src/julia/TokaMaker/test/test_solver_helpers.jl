# Milestone C helper: set_resistivity!.
# (abspsi_to_normalized / psinorm_to_absolute and set_coil_current_dist! are
# exercised on the free-boundary ITER equilibrium in test_ITER.jl, where the
# flux bounds are finite and real coils exist.)

using Test
using TokaMaker

const HELP_DATA_DIR = joinpath(@__DIR__, "data")

function _solo_for_helpers()
    R, a_, b_, c0 = 1.0, 1.2, -1.0, 1.1
    pts, tris, regs, _, _ = load_gs_mesh(joinpath(HELP_DATA_DIR, "solo_h1.h5"))
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
    return gs
end

@testset "set_resistivity! (Milestone C)" begin
    gs = _solo_for_helpers()
    # Piecewise-linear and clearing (nothing) calls should both succeed.
    @test set_resistivity!(gs; eta_prof=Dict("type" => "linterp",
                                              "x" => [0.0, 1.0], "y" => [1.0, 10.0])) === gs
    @test !isfile("tokamaker_eta.prof")   # temp file cleaned up
    @test set_resistivity!(gs; eta_prof=nothing) === gs
end

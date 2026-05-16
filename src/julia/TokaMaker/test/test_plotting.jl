# Smoke-test the Makie extension. Loads CairoMakie (headless backend), builds
# a tiny equilibrium, and verifies plot_machine / plot_psi return without
# error and produce non-empty Figures. We deliberately do not assert pixel-
# exact output; that's brittle across Makie versions.

using Test
using TokaMaker

@testset "Makie extension (plot_machine, plot_psi)" begin
    # CairoMakie is a heavy optional dep; skip the test if unavailable.
    have_cairo = try
        @eval using CairoMakie
        true
    catch
        false
    end
    if !have_cairo
        @info "CairoMakie not installed; skipping plotting tests"
        @test_skip true
    else
        # Verify the extension was loaded
        @test length(methods(plot_machine)) > 0
        @test length(methods(plot_psi)) > 0

        # Build a simple Solov'ev equilibrium
        dom = GsDomain()
        define_region!(dom, "plasma", 0.015, "plasma")
        add_rectangle!(dom, 1.0, 0.0, 0.12, 0.15, "plasma")
        pts, lc, reg = build_mesh!(dom; require_boundary=false)
        env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
        gs = Tokamaker(env)
        setup_mesh!(gs; r=pts, lc=lc, reg=reg)
        setup_regions!(gs)
        gs.settings.free_boundary = false
        gs.settings.pm = false
        setup!(gs; order=2, F0=1.0, full_domain=true)
        set_p_scale!(gs, 1.2); set_ffp_scale!(gs, -2.0)
        set_profiles!(gs; ffp_prof=Dict("type" => "flat"),
                      pp_prof=Dict("type" => "flat"))
        init_psi!(gs; r0=1.0, z0=0.0, a=0.04, kappa=1.0, delta=0.0)
        gs.settings.maxits = 100; update_settings!(gs)
        solve!(gs)

        @testset "plot_machine returns Figure" begin
            fig, ax = plot_machine(gs)
            @test fig isa Makie.Figure
            @test ax isa Makie.Axis
        end

        @testset "plot_psi returns Figure" begin
            fig, ax = plot_psi(gs)
            @test fig isa Makie.Figure
            @test ax isa Makie.Axis
        end
    end
end

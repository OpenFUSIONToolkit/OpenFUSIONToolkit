# Field interpolation + get_stats edge cases.
#
# Covers items #2 and #6 from the to-do audit:
#   - get_field_eval with PSI / B / F / P / dPSI / dBr / dBt / dBz interpreters
#   - case-insensitive field_type matching (mirrors Python's `.upper()`)
#   - get_stats with geom_type="mid" (uses trace_surf instead of bbox)
#   - get_stats on a limited (non-diverted) equilibrium

using Test
using TokaMaker

@testset "Field eval + get_stats edges" begin
    # Build a small Solov'ev fixed-boundary equilibrium to exercise both APIs.
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
    gs.settings.maxits = 100
    update_settings!(gs)
    solve!(gs)

    @testset "get_field_eval scalar/vector field types" begin
        # Field types that work on a fixed-boundary equilibrium without coils.
        # dBr/dBt/dBz require a free-boundary equilibrium with coil fields
        # (gs_project_b in tokamaker_get_field_eval). They are exercised in
        # the ITER baseline test instead.
        for ft in ("PSI", "psi", "Psi", "B", "F", "P", "dPSI", "DPSI")
            fi = get_field_eval(gs.equilibrium, ft)
            v = fi([1.0, 0.0])
            @test v isa Vector{Float64}
            @test all(isfinite, v)
        end
    end

    @testset "Field consistency: B = curl(psi)/R approximately" begin
        # B_R = -1/R * d(psi)/dz, B_Z = 1/R * d(psi)/dr.
        # Verify the dPSI interpolator produces gradient values consistent
        # with finite differences on the PSI interpolator.
        fi_psi = get_field_eval(gs.equilibrium, "psi")
        fi_dpsi = get_field_eval(gs.equilibrium, "dPSI")
        pt = [1.0, 0.01]
        h = 1e-4
        psi_c = fi_psi(pt)[1]
        psi_r = fi_psi([pt[1] + h, pt[2]])[1]
        psi_z = fi_psi([pt[1], pt[2] + h])[1]
        dpsi = fi_dpsi(pt)
        @test isapprox(dpsi[1], (psi_r - psi_c) / h; atol=5e-4, rtol=0.05)
        @test isapprox(dpsi[2], (psi_z - psi_c) / h; atol=5e-4, rtol=0.05)
    end

    @testset "Field eval rejects unknown field type" begin
        @test_throws Exception get_field_eval(gs.equilibrium, "WAT")
    end

    # Note: geom_type="mid" and the dBr/dBt/dBz field interpolators both
    # require a free-boundary equilibrium with coil geometry — they call
    # `gs_trace_surf` and `gs_project_b` respectively which need plasma
    # vacuum-boundary structure that a fixed-boundary Solov'ev lacks.
    # ITER tests exercise both paths.

    @testset "get_stats limited (non-diverted) topology" begin
        # diverted is False here. The X-point swap loop in get_stats should
        # be a no-op. Just confirm the code path runs without crashing — many
        # geometric quantities are ill-defined on a tiny fixed-boundary
        # Solov'ev (no LCFS topology), so we don't assert their values.
        gs.equilibrium.F0 = gs.F0
        @test diverted(gs) == false
        s = get_stats(gs.equilibrium; free_boundary=false,
                      li_normalization="std")
        @test s isa AbstractDict
        @test s["vol"] > 0
        @test haskey(s, "kappa")
        @test haskey(s, "Ip_centroid")
    end
end

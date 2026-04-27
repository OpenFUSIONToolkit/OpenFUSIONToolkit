# Native triangulation tests. Verifies the Julia GsDomain pipeline produces
# meshes that drive a converged Solov'ev solve and an ITER baseline match
# without going through Python.

using Test
using LinearAlgebra: norm
using JSON3
using TokaMaker

const ITER_GEOM = abspath(joinpath(@__DIR__, "..", "..", "..", "tests",
                                   "physics", "ITER_geom.json"))

function build_native_solo_mesh(resolution::Real)
    dom = GsDomain()
    define_region!(dom, "plasma", resolution, "plasma")
    add_rectangle!(dom, 1.0, 0.0, 0.12, 0.15, "plasma")
    return build_mesh!(dom; require_boundary=false)
end

@testset "Native meshing" begin
    @testset "Solov'ev rectangle" begin
        pts, lc, reg = build_native_solo_mesh(0.015)
        @test size(pts, 1) > 100
        @test size(lc, 1) > 100
        @test size(pts, 2) == 2
        @test size(lc, 2) == 3
        @test minimum(reg) == 1
        # Solve and check convergence to Solov'ev exact solution
        R, a_, b_, c0 = 1.0, 1.2, -1.0, 1.1
        env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
        gs = Tokamaker(env)
        setup_mesh!(gs; r=pts, lc=lc, reg=reg)
        setup_regions!(gs)
        gs.settings.free_boundary = false
        gs.settings.pm = false
        setup!(gs; order=2, F0=1.0, full_domain=true)
        set_p_scale!(gs, a_); set_ffp_scale!(gs, b_ * R^2 * 2.0)
        set_profiles!(gs; ffp_prof=Dict("type" => "flat"),
                      pp_prof=Dict("type" => "flat"))
        init_psi!(gs)
        ζ = (gs.r[:, 1].^2 .- R^2) ./ (2R)
        ψ = (b_ + c0) .* R^2 .* gs.r[:, 2].^2 ./ 2 .+
            c0 .* R .* ζ .* gs.r[:, 2].^2 .+ (a_ - c0) .* R^2 .* ζ.^2 ./ 2
        set_psi!(gs, -ψ)
        gs.settings.nl_tol = 1e-14; gs.settings.maxits = 40
        update_settings!(gs)
        solve!(gs)
        psi_err = norm(get_psi(gs; normalized=false) .+ ψ)
        # Native mesh has different point count than Python; allow same order
        # of magnitude as Python's golden 3.2e-7 for h1 order 2.
        @test psi_err < 5e-7
    end

    if isfile(ITER_GEOM)
        @testset "ITER multi-region (annulus + polygon + rectangle)" begin
            geom = JSON3.read(read(ITER_GEOM, String))
            dom = GsDomain()
            define_region!(dom, "air", 0.6, "boundary")
            define_region!(dom, "plasma", 0.15, "plasma")
            define_region!(dom, "vacuum1", 0.3, "vacuum")
            define_region!(dom, "vacuum2", 0.3, "vacuum")
            define_region!(dom, "vv1", 0.3, "conductor"; eta=6.9e-7)
            define_region!(dom, "vv2", 0.3, "conductor"; eta=6.9e-7)
            for k in keys(geom["coils"])
                sk = String(k)
                startswith(sk, "VS") && continue
                define_region!(dom, sk, 0.2, "coil")
            end
            define_region!(dom, "VSU", 0.2, "coil"; coil_set="VS", nTurns=1.0)
            define_region!(dom, "VSL", 0.2, "coil"; coil_set="VS", nTurns=-1.0)
            _to_mat(a) = Matrix{Float64}(reduce(hcat, [Float64.(p) for p in a])')
            add_polygon!(dom, _to_mat(geom["limiter"]), "plasma";
                         parent_name="vacuum1")
            add_annulus!(dom, _to_mat(geom["inner_vv"][1]), "vacuum1",
                         _to_mat(geom["inner_vv"][2]), "vv1";
                         parent_name="vacuum2")
            add_annulus!(dom, _to_mat(geom["outer_vv"][1]), "vacuum2",
                         _to_mat(geom["outer_vv"][2]), "vv2";
                         parent_name="air")
            for (k, c) in geom["coils"]
                sk = String(k)
                parent = startswith(sk, "VS") ? "vacuum1" : "air"
                add_rectangle!(dom, c["rc"], c["zc"], c["w"], c["h"], sk;
                               parent_name=parent)
            end
            pts, lc, reg = build_mesh!(dom)
            # Within 5% of Python's count of (4757, 9400)
            @test 4500 <= size(pts, 1) <= 5000
            @test 9000 <= size(lc, 1) <= 9800
            @test maximum(reg) == 20
        end
    else
        @info "ITER_geom.json not found at $ITER_GEOM; skipping ITER mesh test"
    end
end

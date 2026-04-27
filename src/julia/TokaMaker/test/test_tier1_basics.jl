# Tier 1 API additions: get_targets, vac_solve!, compute_area_integral,
# compute_flux_integral. Exercised on the cached ITER mesh.

using Test
using TokaMaker

const ITER_MESH = joinpath(@__DIR__, "data", "ITER_mesh.h5")

function _build_iter_baseline()
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(ITER_MESH)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
    setup!(gs; order=2, F0=5.3 * 6.2)

    set_coil_vsc!(gs, Dict("VS" => 1.0))
    set_coil_bounds!(gs, Dict(name => [-50e6, 50e6] for name in keys(gs.coil_sets)))
    set_targets!(gs; Ip=15.6e6, pax=6.2e5)

    isoflux = Float64[
        8.20  0.41; 8.06  1.46; 7.51  2.62; 6.14  3.78
        4.51  3.02; 4.26  1.33; 4.28  0.08; 4.49 -1.34
        7.28 -1.89; 8.00 -0.68; 5.125 -3.4
    ]
    set_isoflux_constraints!(gs, isoflux)
    set_saddle_constraints!(gs, Float64[5.125 -3.4][:, :])

    reg_terms = Any[]
    for n in keys(gs.coil_sets)
        if startswith(n, "CS1")
            push!(reg_terms, (coffs=Dict(n => 1.0), target=0.0, weight=2e-2))
        elseif startswith(n, "CS") || startswith(n, "PF") || startswith(n, "VS")
            push!(reg_terms, (coffs=Dict(n => 1.0), target=0.0, weight=1e-2))
        end
    end
    push!(reg_terms, (coffs=Dict("#VSC" => 1.0), target=0.0, weight=1e2))
    set_coil_reg!(gs; reg_terms=reg_terms)

    set_profiles!(gs;
        ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
        pp_prof=create_power_flux_fun(40, 4.0, 1.0))
    init_psi!(gs; r0=6.3, z0=0.5, a=2.0, kappa=1.4, delta=0.0)
    solve!(gs)
    return gs
end

@testset "Tier 1 API" begin
    if !isfile(ITER_MESH)
        @info "Skipping Tier 1 tests: $(ITER_MESH) not found"
        return
    end
    gs = _build_iter_baseline()

    @testset "get_targets" begin
        tg = get_targets(gs)
        @test tg["Ip"] ≈ 15.6e6
        @test tg["pax"] ≈ 6.2e5
        @test !haskey(tg, "estore")

        # retain_previous=false clears
        set_targets!(gs; Ip=12e6)
        tg2 = get_targets(gs)
        @test tg2["Ip"] ≈ 12e6
        @test !haskey(tg2, "pax")

        # retain_previous=true preserves prior values
        set_targets!(gs; pax=5.5e5, retain_previous=true)
        tg3 = get_targets(gs)
        @test tg3["Ip"] ≈ 12e6
        @test tg3["pax"] ≈ 5.5e5

        # V0 / Z0 alias
        set_targets!(gs; Ip=15.6e6, pax=6.2e5, V0=0.5)
        tg4 = get_targets(gs)
        @test tg4["V0"] ≈ 0.5
        @test tg4["Z0"] ≈ 0.5
    end

    @testset "compute_area_integral" begin
        ones_field = ones(Float64, gs.np)
        A_plasma = compute_area_integral(gs, ones_field; reg_mask=1)
        A_total  = compute_area_integral(gs, ones_field; reg_mask=-1)
        @test A_plasma > 0
        @test A_total > A_plasma
        @test 10 < A_plasma < 100
        @test 100 < A_total < 1000
        @test_throws ErrorException compute_area_integral(gs, ones(gs.np + 1))
    end

    @testset "compute_flux_integral" begin
        psi_grid = collect(range(0.0, 1.0; length=33))
        press = @. 6.2e5 * (1 - psi_grid^2)^2
        P_int = compute_flux_integral(gs, psi_grid, press)
        @test P_int > 0
        @test_throws ErrorException compute_flux_integral(gs, [0.0, 1.0], [1.0])
    end

    @testset "vac_solve!" begin
        cur = Dict{String,Float64}(name => 0.0 for name in keys(gs.coil_sets))
        first_pf = first(filter(n -> startswith(n, "PF"), collect(keys(cur))))
        cur[first_pf] = 1e6
        set_coil_currents!(gs, cur)
        vac_solve!(gs)
        psi = get_psi(gs)
        @test maximum(abs, psi) > 0
        @test length(psi) == gs.np
    end
end

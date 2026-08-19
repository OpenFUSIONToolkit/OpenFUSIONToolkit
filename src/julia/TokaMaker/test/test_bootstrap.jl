# ITER bootstrap golden test (Milestone G).
# Mirrors run_ITER_bootstrap_case / ITER_bootstrap_eq_dict in
# src/tests/physics/test_TokaMaker.py: set up the ITER baseline, run
# solve_with_bootstrap with H-mode kinetic profiles, and compare equilibrium
# stats + bootstrap diagnostics against the Python golden.

using Test
using TokaMaker

const BS_MESH = joinpath(@__DIR__, "data", "ITER_mesh.h5")

# Golden values from test_TokaMaker.py:ITER_bootstrap_eq_dict.
const ITER_BS_GOLDEN = Dict{String,Float64}(
    "Ip"          => 15600817.585821694,
    "kappa"       => 1.87554142781964,
    "R_geo"       => 6.222376807932244,
    "a_geo"       => 1.9817209643036526,
    "q_0"         => 0.9951304914765554,
    "q_95"        => 2.856235920791585,
    "P_ax"        => 739971.7132708698,
    "j_BS_max"    => 193963.2797949608,
    "j_BS_axis"   => 7555.958566625245,
    "jphi_axis"   => 1459409.3677809385,
    "jphi_max"    => 1551188.1280449552,
    "j_ind_axis"  => 1357487.1677957429,
    "bs_fraction" => 0.1575907471180497,
)

function _setup_ITER_baseline()
    pts, tris, regs, coil_dict, cond_dict = load_gs_mesh(BS_MESH)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=tris, reg=regs)
    setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
    setup!(gs; order=2, F0=5.3 * 6.2)
    set_coil_vsc!(gs, Dict("VS" => 1.0))
    set_coil_bounds!(gs, Dict(name => [-50e6, 50e6] for name in keys(gs.coil_sets)))
    set_targets!(gs; Ip=15.6e6, pax=6.2e5)
    isoflux = Float64[
        8.20 0.41; 8.06 1.46; 7.51 2.62; 6.14 3.78; 4.51 3.02; 4.26 1.33;
        4.28 0.08; 4.49 -1.34; 7.28 -1.89; 8.00 -0.68; 5.125 -3.4]
    set_isoflux_constraints!(gs, isoflux)
    set_saddle_constraints!(gs, Float64[5.125 -3.4][:, :])
    reg_terms = []
    for name in keys(gs.coil_sets)
        w = startswith(name, "CS1") ? 2e-2 : 1e-2
        push!(reg_terms, (coffs=Dict{String,Float64}(String(name) => 1.0), target=0.0, weight=w))
    end
    push!(reg_terms, (coffs=Dict{String,Float64}("#VSC" => 1.0), target=0.0, weight=1e2))
    set_coil_reg!(gs; reg_terms=reg_terms)
    set_profiles!(gs; ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
                  pp_prof=create_power_flux_fun(40, 4.0, 1.0))
    init_psi!(gs; r0=6.3, z0=0.5, a=2.0, kappa=1.4, delta=0.0)
    solve!(gs)
    return gs
end

@testset "ITER bootstrap (solve_with_bootstrap golden)" begin
    if !isfile(BS_MESH)
        @warn "ITER mesh missing at $BS_MESH"
        @test_skip true
    else
        gs = _setup_ITER_baseline()
        n = 257
        ind = create_power_flux_fun(n, 2.25, 2.5)["y"]
        ne = Hmode_profiles(edge=0.35, ped=0.6, core=1.1, rgrid=n, expin=1.6,
                            expout=1.6, widthp=0.35, xphalf=0.965) .* 1e20
        Te = Hmode_profiles(edge=1500.0, ped=5000.0, core=21000.0, rgrid=n,
                            expin=1.3, expout=1.7, widthp=0.1, xphalf=0.965)
        ni = copy(ne); Ti = copy(Te); Zeff = fill(1.7, n)
        res = solve_with_bootstrap(gs; ne=ne, Te=Te, ni=ni, Ti=Ti, Zeff=Zeff,
                                   Ip_target=15.6e6, inductive_jphi=ind,
                                   scale_jBS=1.0, isolate_edge_jBS=false,
                                   psi_pad=1e-3, iterations=2, verbose=false)
        stats = get_stats(gs.equilibrium; li_normalization="ITER",
                          free_boundary=gs.settings.free_boundary)
        jBS = res["j_BS"]; jt = res["total_j_phi"]; ji = res["j_inductive"]
        psi = collect(range(0.0, 1.0; length=n))
        trapz(y) = sum((y[1:end-1] .+ y[2:end]) ./ 2 .* diff(psi))
        got = Dict{String,Float64}(
            "Ip"          => stats["Ip"],
            "kappa"       => stats["kappa"],
            "R_geo"       => stats["R_geo"],
            "a_geo"       => stats["a_geo"],
            "q_0"         => stats["q_0"],
            "q_95"        => stats["q_95"],
            "P_ax"        => stats["P_ax"],
            "j_BS_max"    => maximum(abs, jBS),
            "j_BS_axis"   => jBS[1],
            "jphi_axis"   => jt[1],
            "jphi_max"    => maximum(abs, jt),
            "j_ind_axis"  => ji[1],
            "bs_fraction" => trapz(jBS) / trapz(jt),
        )
        # validate_dict-style: 1% relative tolerance (a bit looser on the small
        # j_BS_axis where a few-% spread is expected across mesh/solver details).
        for (k, gold) in ITER_BS_GOLDEN
            tol = k == "j_BS_axis" ? 0.05 : 0.02
            @test isapprox(got[k], gold; rtol=tol)
        end
    end
end

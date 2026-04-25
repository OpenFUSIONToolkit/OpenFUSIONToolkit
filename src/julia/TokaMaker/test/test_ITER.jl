# ITER baseline equilibrium test.
# Mirrors run_ITER_case (eig=false, stab=false, recon=false) in
# src/tests/physics/test_TokaMaker.py:462-678.

using Test
using OrderedCollections: OrderedDict
using TokaMaker

const ITER_DATA_DIR = joinpath(@__DIR__, "data")
const ITER_MESH = joinpath(ITER_DATA_DIR, "ITER_mesh.h5")

# Golden values from test_TokaMaker.py:701-726.
const ITER_EQ_DICT = Dict{String,Any}(
    "Ip" => 15599996.692463942,
    "Ip_centroid" => [6.20273409, 0.52959503],
    "kappa" => 1.8728151512244395,
    "kappaU" => 1.7634853298971116,
    "kappaL" => 1.9821449725517677,
    "delta" => 0.4721203463868737,
    "deltaL" => 0.5390229746861446,
    "R_geo" => 6.222328618622752,
    "a_geo" => 1.9835670211775267,
    "vol" => 820.212921617247,
    "q_0" => 0.8232444101221106,
    "q_95" => 2.7602989886308738,
    "P_ax" => 619225.017325726,
    "W_MHD" => 242986393.97329777,
    "beta_pol" => 42.443587820489455,
    "dflux" => 1.540293464599462,
    "tflux" => 121.86081608235014,
    "l_i" => 0.9054096856166233,
    "beta_tor" => 1.7798144109869558,
    "beta_n" => 1.1951205307278518,
    "LCS1" => 2.4858609418809336e-06,
    "MCS1_plasma" => 8.931779419000401e-07,
    "Lplasma" => 1.1900576990802187e-05,
)

# Tolerance per validate_dict (test_TokaMaker.py:47-78).
const ITER_TOL = Dict{String,Float64}("delta" => 5e-2, "deltaU" => 5e-2, "deltaL" => 5e-2)

function _check_within_tol(results::AbstractDict, expected::AbstractDict)
    bad = String[]
    for (k, exp) in expected
        haskey(results, k) || (push!(bad, "$k missing"); continue)
        actual = results[k]
        rel = get(ITER_TOL, k, 1e-2)
        if exp isa AbstractArray
            for i in eachindex(exp)
                if abs((actual[i] - exp[i]) / exp[i]) > rel
                    push!(bad, "$k[$i]: expected $(exp[i]), got $(actual[i])")
                end
            end
        else
            if abs((actual - exp) / exp) > rel
                push!(bad, "$k: expected $exp, got $actual (rel=$(round(abs((actual-exp)/exp); sigdigits=3)))")
            end
        end
    end
    return bad
end

function build_coil_reg_term(; coffs::AbstractDict, target::Real=0.0, weight::Real=1.0)
    return (coffs=Dict{String,Float64}(String(k) => Float64(v) for (k, v) in coffs),
            target=Float64(target), weight=Float64(weight))
end

function run_ITER_case(fe_order::Integer)
    pts, tris, regs, coil_dict, cond_dict = load_gs_mesh(ITER_MESH)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=tris, reg=regs)
    setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
    setup!(gs; order=fe_order, F0=5.3 * 6.2)

    set_coil_vsc!(gs, Dict("VS" => 1.0))
    set_coil_bounds!(gs,
        Dict(name => [-50e6, 50e6] for name in keys(gs.coil_sets)))

    Ip_target = 15.6e6
    P0_target = 6.2e5
    set_targets!(gs; Ip=Ip_target, pax=P0_target)

    isoflux = Float64[
        8.20  0.41
        8.06  1.46
        7.51  2.62
        6.14  3.78
        4.51  3.02
        4.26  1.33
        4.28  0.08
        4.49 -1.34
        7.28 -1.89
        8.00 -0.68
        5.125 -3.4   # X-point appended (mirrors Python test:539)
    ]
    x_point = Float64[5.125  -3.4][:, :]
    set_isoflux_constraints!(gs, isoflux)
    set_saddle_constraints!(gs, x_point)

    # Regularization terms - mirror Python lines 542-553
    reg_terms = []
    for name in keys(gs.coil_sets)
        if startswith(name, "CS1")
            push!(reg_terms, build_coil_reg_term(; coffs=Dict(name => 1.0), weight=2e-2))
        elseif startswith(name, "CS")
            push!(reg_terms, build_coil_reg_term(; coffs=Dict(name => 1.0), weight=1e-2))
        elseif startswith(name, "PF") || startswith(name, "VS")
            push!(reg_terms, build_coil_reg_term(; coffs=Dict(name => 1.0), weight=1e-2))
        end
    end
    push!(reg_terms, build_coil_reg_term(; coffs=Dict("#VSC" => 1.0), weight=1e2))
    set_coil_reg!(gs; reg_terms=reg_terms)

    set_profiles!(gs;
        ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
        pp_prof=create_power_flux_fun(40, 4.0, 1.0))

    init_psi!(gs; r0=6.3, z0=0.5, a=2.0, kappa=1.4, delta=0.0)
    solve!(gs)

    eq = gs.equilibrium
    eq.F0 = gs.F0
    stats = get_stats(eq; li_normalization="ITER", free_boundary=gs.settings.free_boundary)
    Lp, M_pc = calc_inductance(eq, gs.ncoils)
    cs1u_id = gs.coil_sets["CS1U"]["id"]::Int
    stats["LCS1"] = gs.Lcoils[cs1u_id + 1, cs1u_id + 1]
    stats["MCS1_plasma"] = M_pc[cs1u_id + 1]
    stats["Lplasma"] = Lp
    return stats
end

@testset "ITER baseline" begin
    if !isfile(ITER_MESH)
        @warn "ITER mesh missing at $ITER_MESH; generate via Python pipeline"
        @test_skip true
    else
        for order in (2,)
            @testset "order=$order" begin
                stats = run_ITER_case(order)
                bad = _check_within_tol(stats, ITER_EQ_DICT)
                isempty(bad) || foreach(b -> @info("FAIL: $b"), bad)
                @test isempty(bad)
            end
        end
    end
end

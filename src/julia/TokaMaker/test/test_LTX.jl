# LTX equilibrium / eig / stability tests (Milestone G).
# Mirrors run_LTX_case (test_type='' / 'eig' / 'stab') in
# src/tests/physics/test_TokaMaker.py, using the shared LTX_mesh.h5 fixture
# (generated from Python's gs_Domain via the same construction as create_mesh).

using Test
using TokaMaker

const LTX_MESH = joinpath(@__DIR__, "data", "LTX_mesh.h5")

# Golden from test_TokaMaker.py:LTX_eq_dict (delta* entries are commented out
# there; Ip_centroid[2] is None and not compared).
const LTX_EQ_GOLDEN = Dict{String,Any}(
    "Ip"        => 90002.51679781199,
    "kappa"     => 1.525595596236063,
    "kappaU"    => 1.5256161060199729,
    "kappaL"    => 1.5255750864521527,
    "R_geo"     => 0.39198831687969443,
    "a_geo"     => 0.2378387910900877,
    "vol"       => 0.6507554668762836,
    "q_0"       => 1.3280540982807334,
    "q_95"      => 5.901997513881755,
    "P_ax"      => 1720.958666106632,
    "W_MHD"     => 563.2852958452944,
    "beta_pol"  => 41.4052788464515,
    "dflux"     => 0.0009601908294886685,
    "tflux"     => 0.08551496495989133,
    "l_i"       => 1.0271521431711803,
    "beta_tor"  => 1.9276444168027145,
    "beta_n"    => 1.3972402635015146,
)
const LTX_IP_CENTROID_R = 0.40545876749

# Golden eigenvalues from test_TokaMaker.py:test_LTX_eig / test_LTX_stability.
const LTX_TAU_W = [5.152566e-3, 3.953030e-3, 2.536384e-3, 2.172948e-3, 1.853882e-3]
const LTX_GAMMA = [234.1051, -214.4196, -282.0877, -388.7592, -388.7592]

build_reg(coffs; target=0.0, weight=1.0) =
    (coffs=Dict{String,Float64}(String(k) => Float64(v) for (k, v) in coffs),
     target=Float64(target), weight=Float64(weight))

# Mesh setup + FEM (no equilibrium solve): enough for wall-mode analysis.
function build_LTX(fe_order::Integer)
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(LTX_MESH)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
    setup!(gs; order=fe_order, F0=0.10752)
    return gs
end

# Everything from VSC definition through the first equilibrium solve (mirrors
# test_TokaMaker.py:846-876). The stability spectrum is taken here.
function solve_LTX!(gs::Tokamaker)
    set_coil_vsc!(gs, Dict("INTERNALU" => 1.0, "INTERNALL" => -1.0))
    set_targets!(gs; Ip=8.0e4, Ip_ratio=2.0)
    set_isoflux_constraints!(gs, create_isoflux(20, 0.40, 0.0, 0.22, 1.5, 0.1))
    # Regularization terms (mirror test_TokaMaker.py:853-866).
    disable_list = ("YELLOW",)
    reg_terms = []
    for name in keys(gs.coil_sets)
        if name[1:end-1] in disable_list
            push!(reg_terms, build_reg(Dict(name => 1.0); weight=1e4)); continue
        end
        if name == "OH"
            push!(reg_terms, build_reg(Dict(name => 1.0); weight=1e-1)); continue
        elseif name[end] == 'L'
            continue
        end
        push!(reg_terms, build_reg(Dict(name => 1.0); weight=1e-1))
        push!(reg_terms, build_reg(Dict(name => 1.0, name[1:end-1] * "L" => -1.0); weight=1e2))
    end
    push!(reg_terms, build_reg(Dict("#VSC" => 1.0); weight=1e-4))
    set_coil_reg!(gs; reg_terms=reg_terms)
    set_profiles!(gs; ffp_prof=create_power_flux_fun(50, 1.5, 2.0),
                  pp_prof=create_power_flux_fun(50, 4.0, 1.0))
    init_psi!(gs; r0=0.42, z0=0.0, a=0.15, kappa=1.5, delta=0.6)
    solve!(gs)
    return gs
end

function run_LTX_eq(fe_order::Integer)
    gs = build_LTX(fe_order)
    solve_LTX!(gs)
    # Second solve with eddy-current time step and raised Ip target (the golden
    # is taken after this step; mirrors test_TokaMaker.py:883-888).
    psi_last = get_psi(gs; normalized=false)
    set_psi_dt!(gs, psi_last, 5.0e-3)
    set_targets!(gs; Ip=9.0e4, Ip_ratio=2.0)
    solve!(gs)
    eq = gs.equilibrium; eq.F0 = gs.F0
    return get_stats(eq; li_normalization="std", free_boundary=gs.settings.free_boundary)
end

@testset "LTX equilibrium" begin
    if !isfile(LTX_MESH)
        @warn "LTX mesh missing at $LTX_MESH"
        @test_skip true
    else
        for order in (2, 3)
            @testset "order=$order" begin
                stats = run_LTX_eq(order)
                for (k, gold) in LTX_EQ_GOLDEN
                    @test isapprox(stats[k], gold; rtol=1e-2)
                end
                @test isapprox(stats["Ip_centroid"][1], LTX_IP_CENTROID_R; rtol=1e-2)
            end
        end
    end
end

# Wall L/R time constants (Python compute_wall_modes -> 1/eig_vals[:,0]).
@testset "LTX wall eigenvalues" begin
    if !isfile(LTX_MESH)
        @test_skip true
    else
        for order in (2, 3)
            @testset "order=$order" begin
                gs = build_LTX(order)
                res = eig_wall(gs; neigs=10)
                Tau_w = 1.0 ./ res.eigs[1:5, 1]
                for (k, gold) in enumerate(LTX_TAU_W)
                    @test isapprox(Tau_w[k], gold; rtol=1e-2)
                end
            end
        end
    end
end

# Linearized stability spectrum (Python compute_linear_stability(1e3,10,False),
# taken after the first equilibrium solve). gamma = eig_vals[:,0] (no sign flip
# since omega is supplied positive, matching compute_linear_stability).
@testset "LTX linear stability" begin
    if !isfile(LTX_MESH)
        @test_skip true
    else
        for order in (2, 3)
            @testset "order=$order" begin
                gs = build_LTX(order)
                solve_LTX!(gs)
                res = eig_td(gs; omega=1.0e3, neigs=10, include_bounds=false)
                gamma = res.eigs[1:5, 1]
                for (k, gold) in enumerate(LTX_GAMMA)
                    @test isapprox(gamma[k], gold; rtol=1e-2)
                end
            end
        end
    end
end

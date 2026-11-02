# Spheromak analytic verification (Milestone G).
# Mirrors run_sph_case / test_spheromak_h1 in src/tests/physics/test_TokaMaker.py.
# Uses the shared mesh fixture (spheromak_h1.h5, generated from Python's
# gs_Domain) so the solve is identical to Python's and matches its goldens.

using Test
using LinearAlgebra: norm, dot
using SpecialFunctions: besselj
using TokaMaker

const SPH_MESH = joinpath(@__DIR__, "data", "spheromak_h1.h5")

# First zeros of J1 and J0 (jn_zeros(1,1)[0], jn_zeros(0,1)[0]).
const _J1_ZERO = 3.8317059702075125
const _J0_ZERO = 2.404825557695773

# Cerfon analytic spheromak eigenfunction (test_TokaMaker.py:261-265).
function spheromak_psi(r::AbstractVector, z::AbstractVector, a::Real, h::Real)
    g11 = _J1_ZERO .* r ./ a
    nrm = _J0_ZERO * besselj(1, _J0_ZERO)
    return g11 .* besselj.(1, g11) .* sin.(π .* z ./ h) ./ nrm
end

function run_sph_case(mesh_file::AbstractString, fe_order::Integer)
    pts, lc, reg, _, _ = load_gs_mesh(mesh_file)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs)
    gs.settings.free_boundary = false
    setup!(gs; order=fe_order)
    set_p_scale!(gs, 0.0)
    set_profiles!(gs; ffp_prof=Dict("type" => "linterp", "x" => [0.0, 1.0], "y" => [1.0, 0.0]),
                  pp_prof=Dict("type" => "flat"))
    gs.settings.nl_tol = 1e-12
    gs.settings.maxits = 100
    gs.settings.urf = 0.0
    update_settings!(gs)
    init_psi!(gs)
    solve!(gs)
    psi_TM = get_psi(gs; normalized=false)
    psi_eig = spheromak_psi(gs.r[:, 1], gs.r[:, 2], 1.0, 1.0)
    psi_TM .*= dot(psi_eig, psi_eig) / dot(psi_TM, psi_eig)
    return norm(psi_TM .- psi_eig) / norm(psi_eig)
end

# Golden psi errors from test_TokaMaker.py:test_spheromak_h1 (orders 2,3,4).
const SPH_H1_GOLDEN = Dict(
    2 => 2.039674417912789e-05,
    3 => 5.103597862537552e-07,
    4 => 8.088772274705608e-09,
)

@testset "Spheromak analytic" begin
    if !isfile(SPH_MESH)
        @warn "spheromak mesh missing at $SPH_MESH"
        @test_skip true
    else
        for order in (2, 3, 4)
            @testset "order=$order" begin
                err = run_sph_case(SPH_MESH, order)
                # validate_sph allows up to 1.1x the expected error.
                @test err <= SPH_H1_GOLDEN[order] * 1.1
                @test err ≈ SPH_H1_GOLDEN[order] rtol = 1e-3
            end
        end
    end
end

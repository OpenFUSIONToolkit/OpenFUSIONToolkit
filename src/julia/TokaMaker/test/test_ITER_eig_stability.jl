# ITER wall eigenvalue test — mirrors test_ITER_eig from
# src/tests/physics/test_TokaMaker.py (lines 684-689).
#
# Runs eig_wall on a freshly-set-up ITER mesh (no solve required) and
# compares the first 5 wall L/R eigenvalues against the Python golden
# values from test_TokaMaker.py.

using Test
using TokaMaker

const ITER_MESH = joinpath(@__DIR__, "data", "ITER_mesh.h5")
const ITER_TAU_W_GOLDEN = [1.51083009, 2.87431718, 3.91493237, 5.23482507, 5.61049374]

@testset "ITER eig_wall" begin
    if !isfile(ITER_MESH)
        @info "Skipping eig_wall test: $(ITER_MESH) not found"
        return
    end
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(ITER_MESH)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
    setup!(gs; order=2, F0=5.3 * 6.2)
    result = eig_wall(gs; neigs=10)
    tau_w = result.eigs[1:5, 1]
    @test size(result.eigs) == (10, 2)
    @test length(tau_w) == 5
    for k in 1:5
        @test tau_w[k] ≈ ITER_TAU_W_GOLDEN[k] rtol = 1e-3
    end
end

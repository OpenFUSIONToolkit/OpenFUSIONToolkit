# Solov'ev analytical equilibrium verification.
# Mirrors run_solo_case in src/tests/physics/test_TokaMaker.py:150-200.

using Test
using LinearAlgebra: norm
using TokaMaker

const SOLO_DATA_DIR = joinpath(@__DIR__, "data")

# Cerfon-Freidberg analytic Solov'ev solution (test_TokaMaker.py:151-158).
function solovev_psi(r_grid::AbstractVector, z_grid::AbstractVector,
                     R::Real, a::Real, b::Real, c0::Real)
    ζ = @. (r_grid^2 - R^2) / (2 * R)
    ψx = (a - c0) * (b + c0)^2 * R^4 / (8 * c0^2)
    ζx = -(b + c0) * R / (2 * c0)
    Zx = sqrt((b + c0) * (a - c0) / (2 * c0^2)) * R
    ψ_grid = @. (b + c0) * R^2 * z_grid^2 / 2 + c0 * R * ζ * z_grid^2 +
                (a - c0) * R^2 * ζ^2 / 2
    return ψ_grid, ψx, [sqrt(ζx * 2 * R + R^2), Zx]
end

function run_solo_case(mesh_file::AbstractString, fe_order::Integer)
    R, a_, b_, c0 = 1.0, 1.2, -1.0, 1.1
    pts, tris, regs, _, _ = load_gs_mesh(mesh_file)

    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=tris, reg=regs)
    setup_regions!(gs)
    gs.settings.free_boundary = false
    gs.settings.pm = false
    setup!(gs; order=fe_order, F0=1.0, full_domain=true)

    # Solov'ev scalings
    set_p_scale!(gs, a_)
    set_ffp_scale!(gs, b_ * R^2 * 2.0)
    set_profiles!(gs;
        ffp_prof=Dict("type" => "flat"),
        pp_prof=Dict("type" => "flat"))

    init_psi!(gs)  # bare init

    # Seed exact Solov'ev solution as initial guess (matching Python)
    psi_solo, _, rz_x = solovev_psi(gs.r[:, 1], gs.r[:, 2], R, a_, b_, c0)
    set_psi!(gs, -psi_solo)

    gs.settings.nl_tol = 1e-14
    gs.settings.maxits = 40
    update_settings!(gs)
    try
        solve!(gs)
    catch e
        @info "solve failed: $e"
        return nothing
    end

    psi_TM = get_psi(gs; normalized=false)
    psi_err = norm(psi_TM .+ psi_solo)

    x_pts, _ = get_xpoints(gs.equilibrium)
    X_err = 0.0
    # Mirror Python: only the first two candidate X-points (above & below).
    nrows = min(size(x_pts, 1), 2)
    for i in 1:nrows
        diff = [x_pts[i, 1] - rz_x[1], x_pts[i, 2] - rz_x[2]]
        if x_pts[i, 2] < 0.0
            diff[2] = x_pts[i, 2] + rz_x[2]
        end
        X_err += norm(diff)
    end
    return (psi_err=psi_err, X_err=X_err)
end

# Golden values from test_TokaMaker.py:222-250 (h1 = mesh_resolution 0.015).
const SOLO_H1_GOLDEN = Dict(
    2 => (psi_err=3.2048631614233643e-07, X_err=0.00014929412629149645),
    3 => (psi_err=8.919954733021135e-10,  X_err=4.659825491095631e-07),
    4 => (psi_err=5.084454462338564e-15,  X_err=4.224329766330554e-12),
)

@testset "Solov'ev h1" begin
    mesh_file = joinpath(SOLO_DATA_DIR, "solo_h1.h5")
    if !isfile(mesh_file)
        @warn "Solov'ev mesh missing at $mesh_file; skipping. Generate via Python."
        @test_skip true
    else
        for order in (2, 3, 4)
            @testset "order=$order" begin
                res = run_solo_case(mesh_file, order)
                @test res !== nothing
                if res !== nothing
                    g = SOLO_H1_GOLDEN[order]
                    # 10% margin per validate_solo (test_TokaMaker.py:203-218)
                    @test abs(res.psi_err) <= abs(g.psi_err) * 1.1
                    @test abs(res.X_err)   <= abs(g.X_err)   * 1.1
                end
            end
        end
    end
end

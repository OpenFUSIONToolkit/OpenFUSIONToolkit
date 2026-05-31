# Coil vacuum-solve analytic parity (Milestone G).
# Mirrors run_coil_case / test_coil_h1 in src/tests/physics/test_TokaMaker.py,
# using the shared coil_h1.h5 fixture (built from Python's gs_Domain). Compares
# the FEM vacuum boundary psi and the coil-coil mutual inductance against an
# analytic reference: a quadrature of the axisymmetric coil Green's function
# (eval_green), replacing Python's scipy.dblquad with a tensor Gauss-Legendre
# rule (the integrands are smooth, so a fixed rule reaches machine precision).

using Test
using LinearAlgebra: norm, dot, SymTridiagonal, eigen
using TokaMaker

const COIL_MESH = joinpath(@__DIR__, "data", "coil_h1.h5")

# Coil geometry (test_TokaMaker.py:338-340).
const CX1, CY1, CDX, CDY = 0.8, 0.8, 0.1, 0.1
const CX2, CY2 = 0.8, 0.4

# Expected FEM-vs-analytic errors for the h1 mesh, no current distribution
# (test_TokaMaker.py:test_coil_h1, orders 2/3/4).
const COIL_PSI_ERR = [9.62847e-03, 7.45252e-04, 4.29408e-05]
const COIL_MUT_ERR = [1.01747e-03, 3.92000e-04, 5.27038e-05]

# Gauss-Legendre nodes/weights on [-1,1] via Golub-Welsch (eigen of the Jacobi
# matrix). Returns (nodes, weights).
function gauss_legendre(n::Integer)
    k = collect(1.0:(n - 1))
    beta = k ./ sqrt.(4 .* k .^ 2 .- 1)
    J = SymTridiagonal(zeros(n), beta)
    E = eigen(J)
    return E.values, 2.0 .* (E.vectors[1, :]) .^ 2
end

# Tensor-product GL nodes/weights over the rectangle centered at (cx,cy) with
# full widths (dx,dy). Returns (pts[N,2], wts[N]) where N = n^2.
function rect_quadrature(cx, cy, dx, dy, n)
    nodes, wts = gauss_legendre(n)
    rs = cx .+ (dx / 2) .* nodes
    zs = cy .+ (dy / 2) .* nodes
    wr = (dx / 2) .* wts
    wz = (dy / 2) .* wts
    pts = Matrix{Float64}(undef, n * n, 2)
    w = Vector{Float64}(undef, n * n)
    idx = 1
    for i in 1:n, j in 1:n
        pts[idx, 1] = rs[i]; pts[idx, 2] = zs[j]
        w[idx] = wr[i] * wz[j]
        idx += 1
    end
    return pts, w
end

# Flux at field points `pts[N,2]` from a uniform unit current density over
# coil1, i.e. green[i] = ∫∫_coil1 G(filament@(rc,zc); field@pt_i). Uses the
# symmetry G(filament; field) = G(field; filament) so each field point needs a
# single eval_green call over the coil1 quadrature nodes.
function coil1_green(pts::AbstractMatrix, c1pts::AbstractMatrix, c1w::AbstractVector)
    n = size(pts, 1)
    out = Vector{Float64}(undef, n)
    for i in 1:n
        g = eval_green(c1pts, [pts[i, 1], pts[i, 2]])  # G(c1 node; filament@pt_i)
        out[i] = dot(c1w, g)
    end
    return out
end

# Analytic coil1-coil2 mutual (test_TokaMaker.py:analytic_mutual): nested
# integral of the Green's function over both coil cross-sections, normalized.
function analytic_mutual(n)
    c1pts, c1w = rect_quadrature(CX1, CY1, CDX, CDY, n)
    c2pts, c2w = rect_quadrature(CX2, CY2, CDX, CDY, n)
    inner = coil1_green(c2pts, c1pts, c1w)          # ∫∫_coil1 G at each coil2 node
    mutual = dot(c2w, inner)
    return mutual * 2.0 * π / (CDX * CDY) / (CDX * CDY)
end

# Boundary error contribution for one edge: select boundary nodes by `mask`,
# sort by coordinate column `sortcol`, optionally drop the first (shared corner
# at R=0), then return (green, psi) at those points.
function edge_green_psi(r, psi, mask, sortcol, drop_first, c1pts, c1w)
    idx = findall(mask)
    pts = r[idx, :]
    psib = psi[idx]
    order = sortperm(pts[:, sortcol])
    pts = pts[order, :]; psib = psib[order]
    if drop_first
        pts = pts[2:end, :]; psib = psib[2:end]
    end
    green = coil1_green(pts, c1pts, c1w)
    return green, psib
end

function run_coil_case(fe_order::Integer; nquad::Integer=24)
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(COIL_MESH)
    env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
    gs = Tokamaker(env)
    setup_mesh!(gs; r=pts, lc=lc, reg=reg)
    setup_regions!(gs; cond_dict=cond_dict, coil_dict=coil_dict)
    setup!(gs; order=fe_order)
    set_coil_currents!(gs, Dict("COIL1" => CDX * CDY))
    vac_eq = vac_solve!(gs)
    psi0 = TokaMaker.EquilibriumModule.get_psi(vac_eq; normalized=false)

    # Coil mutual: Lcoils[1,2] vs the analytic Green's-function integral.
    Mcc = analytic_mutual(nquad)
    mutual_err = abs((gs.Lcoils[1, 2] + Mcc) / Mcc)

    # Boundary psi error vs analytic green on the three solid edges
    # (test_TokaMaker.py:395-401). Julia mesh columns: r[:,1]=R, r[:,2]=Z.
    c1pts, c1w = rect_quadrature(CX1, CY1, CDX, CDY, nquad)
    r = gs.r
    g1, p1 = edge_green_psi(r, psi0, r[:, 2] .== 1.0, 1, true, c1pts, c1w)   # top (z=1), sort R
    g2, p2 = edge_green_psi(r, psi0, r[:, 1] .== 1.0, 2, false, c1pts, c1w)  # right (r=1), sort Z
    g3, p3 = edge_green_psi(r, psi0, r[:, 2] .== 0.0, 1, true, c1pts, c1w)   # bottom (z=0), sort R
    green_full = vcat(g1, g2, g3)
    psi_full = vcat(p1, p2, p3)
    psi_err = norm(green_full .+ psi_full) / norm(green_full)
    return psi_err, mutual_err
end

@testset "Coil analytic parity" begin
    if !isfile(COIL_MESH)
        @warn "coil mesh missing at $COIL_MESH"
        @test_skip true
    else
        for order in (2, 3, 4)
            @testset "order=$order" begin
                psi_err, mutual_err = run_coil_case(order)
                # validate_coil allows up to 1.1x the expected error.
                @test psi_err <= 1.1 * COIL_PSI_ERR[order-1]
                @test mutual_err <= 1.1 * COIL_MUT_ERR[order-1]
            end
        end
    end
end

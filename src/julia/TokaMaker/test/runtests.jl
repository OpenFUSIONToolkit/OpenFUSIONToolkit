using Test
using TokaMaker

@testset "TokaMaker.jl" begin
    @testset "Settings struct ABI" begin
        s = tokamaker_default_settings()
        @test sizeof(s) == 56
        @test s.pm == true
        @test s.free_boundary == true
        @test s.limited_only == false
        @test s.maxits == 40
        @test s.mode == 1
        @test s.urf == 0.2
        @test s.nl_tol ≈ 1e-6
        @test s.lim_zmax ≈ 1e99

        s.pm = false
        @test s.pm == false
        s.urf = 0.5
        @test s.urf == 0.5
    end

    @testset "Recon settings ABI" begin
        rs = tokamaker_recon_default_settings()
        @test rs.fitI == true
        @test rs.fitP == true
        @test rs.fit_Pscale == false
        @test rs.pm == true
        rs.fitR0 = true
        @test rs.fitR0 == true
    end

    @testset "Bootstrap helpers" begin
        # Hmode profile sanity
        prof = Hmode_profiles(edge=0.08, ped=0.4, core=2.5, rgrid=11)
        @test length(prof) == 11
        @test prof[1] ≈ 2.5 rtol = 1e-2
        @test prof[end] ≈ 0.08 rtol = 1e-2

        # Coulomb log
        ln_e, ln_i = calculate_ln_lambda([1000.0, 2000.0], [800.0, 1500.0],
                                          [1e19, 5e19], [1e19, 5e19])
        @test all(ln_e .> 10.0)
        @test all(ln_i .> 10.0)
        @test length(ln_e) == 2

        # Redl bootstrap on a tiny synthetic profile (all arrays length 3)
        n = 3
        psi_N = [0.1, 0.5, 0.9]
        Te = fill(2000.0, n); Ti = fill(1500.0, n)
        ne = fill(1e19, n);  ni = fill(1e19, n)
        pe = ne .* Te .* 1.602e-19
        pi_arr = ni .* Ti .* 1.602e-19
        Zeff = fill(1.5, n); Rmaj = fill(6.2, n); q = fill(2.5, n)
        eps = fill(0.3, n); fT = fill(0.5, n); I_psi = fill(20.0, n)
        dTe = fill(-100.0, n); dTi = fill(-80.0, n); dne = fill(-1e17, n)
        j_BS, coeffs = redl_bootstrap(; psi_N=psi_N, Te=Te, Ti=Ti, ne=ne, ni=ni,
            pe=pe, pi=pi_arr, Zeff=Zeff, R=Rmaj, q=q, eps=eps, fT=fT,
            I_psi=I_psi, dT_e_dpsi=dTe, dT_i_dpsi=dTi, dn_e_dpsi=dne)
        @test length(j_BS) == n
        @test all(isfinite, j_BS)
        @test haskey(coeffs, "L31")
        @test haskey(coeffs, "L32")
        @test haskey(coeffs, "alpha")
    end

    @testset "Mesh I/O round-trip" begin
        mktempdir() do dir
            fn = joinpath(dir, "mesh.h5")
            pts = Float64[1.0 0.0; 2.0 0.0; 1.5 1.0]
            tris = Int32[1 2 3]
            regs = Int32[1]
            cd = Dict("plasma" => Dict("reg_id" => 1))
            cdc = Dict{String,Any}()
            save_gs_mesh(pts, tris, regs, cd, cdc, fn)
            p2, t2, r2, cd2, cdc2 = load_gs_mesh(fn)
            @test p2 == pts
            @test t2 == tris
            @test r2 == regs
            @test haskey(cd2, "plasma")
        end
    end

    @testset "Util helpers" begin
        pts = create_isoflux(8, 6.2, 0.5, 2.0, 1.7, 0.33)
        @test size(pts) == (8, 2)
        @test pts[1, 1] ≈ 6.2 + 2.0 atol = 1e-10
        @test pts[1, 2] ≈ 0.5 atol = 1e-10

        prof = create_power_flux_fun(40, 1.5, 2.0)
        @test prof["type"] == "linterp"
        @test length(prof["x"]) == 40
        @test prof["y"][1] == 1.0
        @test prof["y"][end] == 0.0

        spl = create_spline_flux_fun(20, [0.0, 0.5, 1.0], [1.0, 0.5, 0.0])
        @test spl["type"] == "linterp"
        @test length(spl["x"]) == 20
        @test spl["y"][1] ≈ 1.0
    end

    # Tests below require liboftpy to be built; they short-circuit otherwise.
    if !isempty(TokaMaker.LibPath.liboftpy[])
        include("test_solovev.jl")
        include("test_ITER.jl")
        include("test_meshing.jl")
        include("test_fields_and_stats.jl")
        include("test_ITER_eig_stability.jl")
        include("test_plotting.jl")
        include("test_reconstruction.jl")
        include("test_bootstrap_composers.jl")
        include("test_tier1_basics.jl")
    else
        @info "Skipping runtime tests: liboftpy not located. Build OFT with -DOFT_BUILD_PYTHON=ON."
    end
end

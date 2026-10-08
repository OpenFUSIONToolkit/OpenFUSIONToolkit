# GEQDSK COCOS / serialization tests (Milestone G).
# Mirrors the COCOS round-trip / sign / invalid / flip / save-load / bytes tests
# in src/tests/physics/test_TokaMaker.py (test_eqdsk_cocos_*), using the shared
# ITER_test.eqdsk fixture (COCOS 1).

using Test
using TokaMaker

const EQDSK_FILE = joinpath(@__DIR__, "data", "ITER_test.eqdsk")
const COCOS_VALID = vcat(1:8, 11:18)

# Independent reference of the COCOS sign/exp table (test_TokaMaker.py:1540).
function ref_cocos(idx)
    (idx < 1 || idx > 18 || idx == 9 || idx == 10) && error("bad")
    exp_Bp = idx < 10 ? 0 : 1
    base = idx < 10 ? idx : idx - 10
    (sigma_Bp = base in (1, 2, 5, 6) ? 1 : -1,
     sigma_RpZ = base in (1, 2, 7, 8) ? 1 : -1,
     sigma_rhotp = base in (1, 3, 5, 7) ? 1 : -1,
     exp_Bp = exp_Bp)
end

# Multiplicative factors cocosify() should apply, keyed by raw field name
# (test_TokaMaker.py:_expected_cocos_factors).
function expected_factors(cin, cout)
    a = ref_cocos(cin); b = ref_cocos(cout)
    sBp = b.sigma_Bp * a.sigma_Bp
    sRpZ = b.sigma_RpZ * a.sigma_RpZ
    srhotp = b.sigma_rhotp * a.sigma_rhotp
    twopi = (2π)^(b.exp_Bp - a.exp_Bp)
    psi = sRpZ * sBp * twopi
    dpsi = sRpZ * sBp / twopi
    return Dict("SIMAG" => psi, "SIBRY" => psi, "PSIRZ" => psi,
                "PPRIME" => dpsi, "FFPRIM" => dpsi,
                "FPOL" => sRpZ, "BCENTR" => sRpZ, "CURRENT" => sRpZ,
                "QPSI" => srhotp)
end

@testset "GEQDSK COCOS" begin
    if !isfile(EQDSK_FILE)
        @warn "ITER_test.eqdsk missing at $EQDSK_FILE"
        @test_skip true
    else
        @testset "cocos_params table" begin
            for n in COCOS_VALID
                p = cocos_params(n)
                r = ref_cocos(n)
                @test (p.sigma_Bp, p.sigma_RpZ, p.sigma_rhotp, p.exp_Bp) ==
                      (r.sigma_Bp, r.sigma_RpZ, r.sigma_rhotp, r.exp_Bp)
            end
            for bad in (0, 9, 10, 19, -1)
                @test_throws ErrorException cocos_params(bad)
            end
        end

        @testset "cocos roundtrip 1->N->1" begin
            eq_ref = read_geqdsk(EQDSK_FILE)
            for cout in COCOS_VALID
                eq_rt = cocosify(cocosify(eq_ref, cout), eq_ref.cocos)
                @test eq_rt.cocos == eq_ref.cocos
                for key in ("SIMAG", "SIBRY", "BCENTR", "CURRENT",
                            "FPOL", "PRES", "PPRIME", "FFPRIM", "QPSI", "PSIRZ")
                    @test all(isapprox.(eq_rt.raw[key], eq_ref.raw[key];
                                        rtol=1e-12, atol=1e-12))
                end
            end
        end

        @testset "cocos signs 1->N" begin
            eq_ref = read_geqdsk(EQDSK_FILE)
            for cout in COCOS_VALID
                eq_n = cocosify(eq_ref, cout)
                @test eq_n.cocos == cout
                for (key, fac) in expected_factors(eq_ref.cocos, cout)
                    @test all(isapprox.(eq_n.raw[key], eq_ref.raw[key] .* fac;
                                        rtol=1e-12, atol=1e-12))
                end
            end
        end

        @testset "invalid cocos raises" begin
            for bad in (0, 9, 10, 19, -1)
                @test_throws ErrorException read_geqdsk(EQDSK_FILE; cocos=bad)
                eq = read_geqdsk(EQDSK_FILE)
                @test_throws ErrorException cocosify(eq, bad)
            end
        end

        @testset "flip_Bt_Ip twice is identity" begin
            eq_ref = read_geqdsk(EQDSK_FILE)
            eq_rt = flip_Bt_Ip(flip_Bt_Ip(eq_ref))
            for key in ("BCENTR", "FPOL", "CURRENT", "SIMAG", "SIBRY",
                        "PSIRZ", "PPRIME", "FFPRIM")
                @test all(isapprox.(eq_rt.raw[key], eq_ref.raw[key];
                                    rtol=1e-12, atol=1e-12))
            end
        end

        @testset "save / load roundtrip" begin
            eq = read_geqdsk(EQDSK_FILE)
            mktempdir() do dir
                out = joinpath(dir, "rt.eqdsk")
                save_geqdsk(eq, out)
                eq_rt = read_geqdsk(out)
                for key in ("NW", "NH", "RDIM", "ZDIM", "RCENTR", "RLEFT", "ZMID",
                            "RMAXIS", "ZMAXIS", "SIMAG", "SIBRY", "BCENTR", "CURRENT")
                    @test isapprox(eq.raw[key], eq_rt.raw[key]; rtol=1e-6, atol=0)
                end
                for key in ("FPOL", "PRES", "PPRIME", "FFPRIM", "QPSI", "PSIRZ",
                            "RBBBS", "ZBBBS", "RLIM", "ZLIM")
                    @test all(isapprox.(eq.raw[key], eq_rt.raw[key];
                                        rtol=1e-6, atol=1e-10))
                end
            end
        end

        @testset "bytes roundtrip" begin
            eq = read_geqdsk(EQDSK_FILE)
            eq_rt = eqdsk_from_bytes(eqdsk_to_bytes(eq); cocos=eq.cocos)
            for key in ("FPOL", "PRES", "PSIRZ", "QPSI")
                @test all(isapprox.(eq.raw[key], eq_rt.raw[key];
                                    rtol=1e-6, atol=1e-10))
            end
        end
    end
end

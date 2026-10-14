# EFIT namelist / mhdin / k-file parser tests (Milestone E).
# Mirrors OpenFUSIONToolkit.util.read_fortran_namelist and
# TokaMaker.util.read_mhdin / read_kfile. No fixtures ship for these in the
# Python suite, so we use small self-contained synthetic namelists written to a
# temp dir (the parser is the thing under test).

using Test
using TokaMaker

const _MHDIN = """
&in3
 MPNAM2 = 'MP1' 'MP2'
 LPNAME = 'FL1' 'FL2'
 AMP2 = 0.0 90.0
 ECID = 1 1 2
 RE = 1.0 1.0 2.0
 ZE = 0.1 0.2 0.3
 WE = 0.05 0.05 0.05
 HE = 0.05 0.05 0.05
 FCID = 1 2
 RF = 1.5 1.6
 ZF = 0.5 -0.5
 WF = 0.1 0.1
 HF = 0.1 0.1
 TURNFC = 100.0 200.0
/
"""

const _KFILE = """
&in1
 EXPMP2 = 0.1 0.2
 FWTMP2 = 1.0 1.0
 COILS = 0.5 0.6
 FWTSI = 1.0 0.0
 ECURRT = 1000.0 2000.0
 FWTEC = 1.0 1.0
 BRSP = 500.0 600.0
 FWTFC = 1.0 1.0
/
"""

@testset "EFIT parsers" begin
    @testset "read_fortran_namelist" begin
        mktempdir() do dir
            f = joinpath(dir, "nl.dat")
            # Exercise scalar, array, N*repetition, and string values.
            write(f, """
            &test
             SCALAR_I = 7
             SCALAR_F = 3.5
             ARR = 1.0 2.0 3.0
             REPEAT = 4*0.0
             NAME = 'foo'
            /
            """)
            nl = read_fortran_namelist(f)
            @test nl["SCALAR_I"] == 7
            @test nl["SCALAR_I"] isa Integer
            @test nl["SCALAR_F"] ≈ 3.5
            @test nl["ARR"] == [1.0, 2.0, 3.0]
            @test nl["REPEAT"] == zeros(4)
            @test occursin("foo", String(nl["NAME"]))   # non-numeric kept as string
        end
    end

    @testset "read_mhdin + read_kfile" begin
        mktempdir() do dir
            mfile = joinpath(dir, "mhdin.dat"); write(mfile, _MHDIN)
            kfile = joinpath(dir, "k.dat");     write(kfile, _KFILE)

            machine, raw = read_mhdin(mfile)
            @test machine["MPNAM2"] == ["MP1", "MP2"]
            @test machine["LPNAME"] == ["FL1", "FL2"]
            @test machine["AMP2"]["MP2"] ≈ 90.0
            # ECID = [1,1,2] -> ECOIL001 has 2 sub-entries, ECOIL002 has 1.
            @test haskey(machine["ECOIL"], "ECOIL001")
            @test length(machine["ECOIL"]["ECOIL001"]) == 2
            @test length(machine["ECOIL"]["ECOIL002"]) == 1
            @test machine["ECOIL"]["ECOIL002"][1] ≈ [2.0, 0.3, 0.05, 0.05]
            # FCID length 2 -> FCOIL001, FCOIL002 with [R,Z,W,H,turns].
            @test machine["FCOIL"]["FCOIL002"] ≈ [1.6, -0.5, 0.1, 0.1, 200.0]

            probes, loops, ecoils, fcoils, kraw = read_kfile(kfile, machine)
            @test probes["MP1"] ≈ [0.1, 1.0]
            @test loops["FL2"] ≈ [0.6, 0.0]
            @test ecoils["ECOIL001"] ≈ [1000.0, 1.0]
            @test fcoils["FCOIL002"] ≈ [600.0, 1.0]
        end
    end
end

# ---------------------------------------------------------------------------
# Infinite-n (ideal ballooning) stability boundary driver for the
# TokaMaker -> GPEC negative-triangularity pedestal workflow.
#
# For each gEQDSK given on the command line this script loads the equilibrium
# with GPEC and computes, on every flux surface of the GPEC radial grid:
#
#   alpha           experimental normalized pressure gradient
#   alpha_critical  first (low-alpha) infinite-n ballooning stability boundary
#
# Results are written next to each gEQDSK as
# `ballooning_boundary_<case>.json`, which the notebook reads back.
#
# Usage:
#   julia -t auto --project=$GPEC_PROJECT compute_ballooning_boundary.jl <geqdsk> [<geqdsk> ...]
#
# `ballooning_alpha_boundary` is the fast first-boundary-only driver: the
# alpha scan at each surface stops at the first Delta' zero crossing, so the
# cost scales with the location of the boundary rather than with the full
# alpha range. That is all the pedestal-height prediction needs -- for
# negative triangularity the pedestal has no access to the second stable
# region, so the first boundary is the operative limit.
# ---------------------------------------------------------------------------

using GeneralizedPerturbedEquilibrium
using GeneralizedPerturbedEquilibrium: Equilibrium, LocalStability

# [Equilibrium] settings mirroring GPEC's own `DIIID-like_ideal_example`. The
# `log_asymptotic` auto grid (mpsi=0) packs surfaces toward the edge, which is
# where the pedestal ballooning boundary lives, and uses far fewer surfaces
# than a uniform grid of the same edge resolution.
function make_eq_dict(geqdsk_path)
    return Dict{String,Any}(
        "eq_filename" => geqdsk_path,
        "eq_type" => "efit",
        "jac_type" => "hamada",
        "grid_type" => "log_asymptotic",
        "psilow" => 1e-4,
        "psihigh" => 0.9995,
        "mpsi" => 0,
        "psi_accuracy" => 0.001,
        "mtheta" => 256,
        "newq0" => 0,
        "etol" => 1e-10,
    )
end

# Standard JSON has no NaN/Inf, so non-finite entries (surfaces with no
# boundary in range, or a failed surface) are written as `null` and read back
# as NaN on the python side.
fmt(v) = isfinite(v) ? repr(v) : "null"

function write_json(path, bnd)
    open(path, "w") do f
        write(f, "{\n")
        write(f, "  \"psi\": [", join(fmt.(bnd.psi), ", "), "],\n")
        write(f, "  \"alpha\": [", join(fmt.(bnd.alpha), ", "), "],\n")
        write(f, "  \"alpha_critical\": [", join(fmt.(bnd.alpha_critical), ", "), "]\n")
        write(f, "}\n")
    end
end

isempty(ARGS) && error("Usage: julia -t auto --project=<GPEC repo> compute_ballooning_boundary.jl <geqdsk> [...]")

for geqdsk_file in ARGS
    geqdsk_path = abspath(geqdsk_file)
    isfile(geqdsk_path) || error("File not found: $geqdsk_path")
    case_label = splitext(basename(geqdsk_path))[1]
    println("Ballooning boundary: $case_label")

    eq_config = Equilibrium.EquilibriumConfig(make_eq_dict(geqdsk_path), dirname(geqdsk_path))
    equil = Equilibrium.setup_equilibrium(eq_config)

    # n_scan=64: the pedestal surfaces have closely spaced Delta' poles, where a
    # coarse scan can step over the first crossing. Early stopping keeps the
    # cost of the finer scan modest.
    bnd = LocalStability.ballooning_alpha_boundary(equil; n_scan=64)

    json_path = joinpath(dirname(geqdsk_path), "ballooning_boundary_$(case_label).json")
    write_json(json_path, bnd)
    println("  saved: $json_path  ($(length(bnd.psi)) surfaces)")
end

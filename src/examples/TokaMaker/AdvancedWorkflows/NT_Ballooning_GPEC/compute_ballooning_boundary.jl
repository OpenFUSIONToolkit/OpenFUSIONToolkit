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
#   julia -t auto --project=$GPEC_PROJECT compute_ballooning_boundary.jl [--both] <geqdsk> [...]
#
# With `--both`, the second (upper) ballooning boundary `alpha_critical2` is
# computed as well and added to the JSON. This is the slower.
#
# `ballooning_alpha_boundary` is the fast first-boundary-only driver: the
# alpha scan at each surface stops at the first Delta' zero crossing, so the
# cost scales with the location of the boundary rather than with the full
# alpha range. 
# ---------------------------------------------------------------------------

using GeneralizedPerturbedEquilibrium
using GeneralizedPerturbedEquilibrium: Equilibrium, LocalStability

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

fmt(v) = isfinite(v) ? repr(v) : "null"

function write_json(path, psi, alpha, alpha_critical, alpha_critical2)
    open(path, "w") do f
        write(f, "{\n")
        write(f, "  \"psi\": [", join(fmt.(psi), ", "), "],\n")
        write(f, "  \"alpha\": [", join(fmt.(alpha), ", "), "],\n")
        write(f, "  \"alpha_critical\": [", join(fmt.(alpha_critical), ", "), "]")
        if alpha_critical2 !== nothing
            write(f, ",\n  \"alpha_critical2\": [", join(fmt.(alpha_critical2), ", "), "]\n")
        else
            write(f, "\n")
        end
        write(f, "}\n")
    end
end

both_boundaries = "--both" in ARGS
geqdsk_files = filter(arg -> arg != "--both", ARGS)
isempty(geqdsk_files) && error("Usage: julia -t auto --project=<GPEC repo> compute_ballooning_boundary.jl [--both] <geqdsk> [...]")

for geqdsk_file in geqdsk_files
    geqdsk_path = abspath(geqdsk_file)
    isfile(geqdsk_path) || error("File not found: $geqdsk_path")
    case_label = splitext(basename(geqdsk_path))[1]
    println("Ballooning boundary: $case_label")

    eq_config = Equilibrium.EquilibriumConfig(make_eq_dict(geqdsk_path), dirname(geqdsk_path))
    equil = Equilibrium.setup_equilibrium(eq_config)

    if both_boundaries
        bnd = LocalStability.ballooning_alpha_boundaries(equil; n_scan=64)
        psi, alpha = bnd.psi, bnd.alpha
        alpha_critical, alpha_critical2 = bnd.alpha_critical1, bnd.alpha_critical2
    else
        bnd = LocalStability.ballooning_alpha_boundary(equil; n_scan=64)
        psi, alpha = bnd.psi, bnd.alpha
        alpha_critical, alpha_critical2 = bnd.alpha_critical, nothing
    end

    json_path = joinpath(dirname(geqdsk_path), "ballooning_boundary_$(case_label).json")
    write_json(json_path, psi, alpha, alpha_critical, alpha_critical2)
    println("  saved: $json_path  ($(length(psi)) surfaces)")
end

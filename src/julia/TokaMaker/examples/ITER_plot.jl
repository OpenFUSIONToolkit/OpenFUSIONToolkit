# ITER baseline equilibrium + plots, runnable as a script.
#
# Usage (from src/julia/TokaMaker/):
#   julia --project=. examples/ITER_plot.jl
#
# Requires:
#   - liboftpy built (cmake -DOFT_BUILD_PYTHON=ON in <repo>/builds)
#   - CairoMakie installed in the project env:
#       julia --project=. -e 'using Pkg; Pkg.add("CairoMakie")'
#   - ITER_geom.json from src/tests/physics/

using TokaMaker
using JSON3
using CairoMakie

const ITER_GEOM = abspath(joinpath(@__DIR__, "..", "..", "..", "..",
                                   "src", "tests", "physics", "ITER_geom.json"))

println("Building ITER mesh natively...")
geom = JSON3.read(read(ITER_GEOM, String))
dom = GsDomain()
define_region!(dom, "air",     0.6, "boundary")
define_region!(dom, "plasma",  0.15, "plasma")
define_region!(dom, "vacuum1", 0.3, "vacuum")
define_region!(dom, "vacuum2", 0.3, "vacuum")
define_region!(dom, "vv1",     0.3, "conductor"; eta=6.9e-7)
define_region!(dom, "vv2",     0.3, "conductor"; eta=6.9e-7)
for k in keys(geom["coils"])
    sk = String(k)
    startswith(sk, "VS") && continue
    define_region!(dom, sk, 0.2, "coil")
end
define_region!(dom, "VSU", 0.2, "coil"; coil_set="VS", nTurns=1.0)
define_region!(dom, "VSL", 0.2, "coil"; coil_set="VS", nTurns=-1.0)

_to_mat(a) = Matrix{Float64}(reduce(hcat, [Float64.(p) for p in a])')
add_polygon!(dom, _to_mat(geom["limiter"]), "plasma"; parent_name="vacuum1")
add_annulus!(dom, _to_mat(geom["inner_vv"][1]), "vacuum1",
             _to_mat(geom["inner_vv"][2]), "vv1"; parent_name="vacuum2")
add_annulus!(dom, _to_mat(geom["outer_vv"][1]), "vacuum2",
             _to_mat(geom["outer_vv"][2]), "vv2"; parent_name="air")
for (k, c) in geom["coils"]
    sk = String(k)
    parent = startswith(sk, "VS") ? "vacuum1" : "air"
    add_rectangle!(dom, c["rc"], c["zc"], c["w"], c["h"], sk; parent_name=parent)
end
pts, lc, reg = build_mesh!(dom)
println("  $(size(pts,1)) points, $(size(lc,1)) cells, $(maximum(reg)) regions")

println("\nInitializing solver and solving baseline equilibrium...")
env = OFTEnv(nthreads=2, quiet=true, use_abort_callback=false)
gs = Tokamaker(env)
setup_mesh!(gs; r=pts, lc=lc, reg=reg)
setup_regions!(gs; cond_dict=get_conductors(dom), coil_dict=get_coils(dom))
setup!(gs; order=2, F0=5.3 * 6.2)

set_coil_vsc!(gs, Dict("VS" => 1.0))
set_coil_bounds!(gs, Dict(name => [-50e6, 50e6] for name in keys(gs.coil_sets)))
set_targets!(gs; Ip=15.6e6, pax=6.2e5)

isoflux = Float64[
    8.20  0.41
    8.06  1.46
    7.51  2.62
    6.14  3.78
    4.51  3.02
    4.26  1.33
    4.28  0.08
    4.49 -1.34
    7.28 -1.89
    8.00 -0.68
    5.125 -3.4
]
set_isoflux_constraints!(gs, isoflux)
set_saddle_constraints!(gs, Float64[5.125 -3.4][:, :])

reg_terms = Any[]
for n in keys(gs.coil_sets)
    if startswith(n, "CS1")
        push!(reg_terms, (coffs=Dict(n => 1.0), target=0.0, weight=2e-2))
    elseif startswith(n, "CS") || startswith(n, "PF") || startswith(n, "VS")
        push!(reg_terms, (coffs=Dict(n => 1.0), target=0.0, weight=1e-2))
    end
end
push!(reg_terms, (coffs=Dict("#VSC" => 1.0), target=0.0, weight=1e2))
set_coil_reg!(gs; reg_terms=reg_terms)

set_profiles!(gs; ffp_prof=create_power_flux_fun(40, 1.5, 2.0),
              pp_prof=create_power_flux_fun(40, 4.0, 1.0))

init_psi!(gs; r0=6.3, z0=0.5, a=2.0, kappa=1.4, delta=0.0)
solve!(gs)

gs.equilibrium.F0 = gs.F0
stats = get_stats(gs.equilibrium; li_normalization="ITER")
println("\nEquilibrium stats:")
for k in ("Ip", "kappa", "delta", "R_geo", "a_geo", "vol", "q_0", "q_95",
          "P_ax", "W_MHD", "beta_n", "l_i")
    v = stats[k]
    if v isa Real
        println("  $(rpad(k, 10)) = $(round(v; sigdigits=5))")
    end
end

println("\nPlotting...")
fig = Figure(size=(1200, 600))

ax1 = Axis(fig[1, 1]; title="ITER machine geometry",
           xlabel="R (m)", ylabel="Z (m)")
plot_machine(fig, ax1, gs)

ax2 = Axis(fig[1, 2]; title="Poloidal flux contours",
           xlabel="R (m)", ylabel="Z (m)")
plot_machine(fig, ax2, gs; vacuum_color=nothing, cond_color=:lightgray,
             coil_color=:gray, limiter_color=:black)
plot_psi(fig, ax2, gs)
plot_constraints(fig, ax2, gs; isoflux=isoflux,
                 saddles=Float64[5.125 -3.4][:, :])

mkpath("plots")
out = joinpath("plots", "iter_baseline.png")
save(out, fig)
save_eqdsk(gs, joinpath("plots", "iter_baseline.eqdsk"); nr=129, nz=129)
println("\nSaved: $(abspath(out))")
println("Saved: $(abspath(joinpath("plots", "iter_baseline.eqdsk")))")

module TokaMaker

include("LibPath.jl")
using .LibPath: liboftpy, liboft_triangle, OFT_PATH_SLEN, OFT_ERROR_SLEN

include("Settings.jl")
using .SettingsModule: TokamakerSettings, tokamaker_default_settings,
                       TokamakerReconSettings, tokamaker_recon_default_settings

include("CInterface.jl")
using .CInterface

include("OFTEnv.jl")
using .OFTEnvModule: OFTEnv, oftpy_set_debug, oftpy_set_nthreads

include("Equilibrium.jl")
using .EquilibriumModule: TokaMakerEquilibrium

include("Profiles.jl")
using .Profiles: write_profile_file

include("Core.jl")
using .CoreModule

include("Util.jl")
using .Util: create_isoflux, create_power_flux_fun, create_spline_flux_fun,
             eval_green, read_eqdsk, read_ifile

include("FieldEval.jl")
using .FieldEvalModule: TokamakerFieldInterpolator, get_field_eval

include("TimeDep.jl")
using .TimeDep: setup_td!, step_td!, set_psi_dt!, eig_td, eig_wall

include("Reconstruction.jl")
using .Reconstruction: MirnovCon, IpCon, FluxLoopCon, DFluxCon, PressCon,
                        QCon, SaddleCon, ReconConstraints, run_reconstruction!,
                        write_constraints_file

include("Bootstrap.jl")
using .Bootstrap: Hmode_profiles, parameterize_edge_jBS, calculate_ln_lambda,
                   redl_bootstrap, get_jphi_from_GS, find_peaks, peak_widths,
                   analyze_bootstrap_edge_spike, find_optimal_scale, solve_jphi,
                   solve_with_bootstrap

# Wire Bootstrap's runtime callbacks now that Core/Equilibrium are loaded.
function __init__()
    Bootstrap.set_targets!_anon[]   = (gs, Ip, pax) -> CoreModule.set_targets!(gs; Ip=Ip, pax=pax)
    Bootstrap.set_profiles!_anon[]  = (gs, ffp, pp) -> CoreModule.set_profiles!(gs; ffp_prof=ffp, pp_prof=pp)
    Bootstrap.solve!_anon[]         = (gs)         -> CoreModule.solve!(gs)
    Bootstrap.get_profiles_anon[]   = (eq; npsi, psi_pad) -> EquilibriumModule.get_profiles(eq; npsi=npsi, psi_pad=psi_pad)
    Bootstrap.get_q_anon[]          = (eq; npsi, psi_pad) -> EquilibriumModule.get_q(eq; npsi=npsi, psi_pad=psi_pad, compute_geo=true)
    Bootstrap.get_stats_anon[]      = (eq; lcfs_pad) -> EquilibriumModule.get_stats(eq; lcfs_pad=lcfs_pad)
    Bootstrap.sauter_fc_anon[]      = (eq; npsi, psi_pad) -> EquilibriumModule.calc_sauter_fc(eq; npsi=npsi, psi_pad=psi_pad)
    Bootstrap.psi_bounds_anon[]     = (gs) -> CoreModule.psi_bounds(gs)
    Bootstrap.compute_flux_integral_anon[] = (gs, psi, field) -> CoreModule.compute_flux_integral(gs, psi, field)
    return nothing
end

include("Meshing.jl")
using .Meshing: save_gs_mesh, load_gs_mesh, GsDomain, define_region!,
                 add_polygon!, add_rectangle!, add_annulus!, add_enclosed!,
                 build_mesh!, get_coils, get_conductors, save_json

include("Eqdsk.jl")
using .Eqdsk: GEQDSKEquilibrium, read_geqdsk, cocos_params, cocosify, cocosify!,
              flip_Bt_Ip, flip_Bt_Ip!, save_geqdsk, eqdsk_to_bytes, eqdsk_from_bytes,
              eqdsk_from_raw

"""
    compute_forces_components(t, psi; cell_centered=false)
        -> (J_cond, B_cond, mask, R)

Compute terms needed to evaluate forces in passively conducting regions for the
given `psi`: conductor toroidal current density `J_cond`, the in-conductor
B-field `B_cond`, a per-cell `mask` of conducting cells, and coordinates `R`
(cell centers if `cell_centered`, else node points). Mirrors Python
`util.compute_forces_components`. Defined here (top module) since it composes
`get_conductor_currents` (Core) with `get_field_eval` (FieldEval).
"""
function compute_forces_components(t::Tokamaker, psi::AbstractVector;
                                   cell_centered::Bool=false)
    mask, J_cond = get_conductor_currents(t, psi; cell_centered=cell_centered)
    np = size(t.r, 1)
    # Vertices belonging to conducting cells (t.lc is 0-based -> +1).
    pt_mask = falses(np)
    @inbounds for c in 1:t.nc
        if mask[c]
            pt_mask[t.lc[c, 1] + 1] = true
            pt_mask[t.lc[c, 2] + 1] = true
            pt_mask[t.lc[c, 3] + 1] = true
        end
    end
    # Evaluate B for the supplied psi in conducting regions, then restore psi.
    psi_save = get_psi(t; normalized=false)
    set_psi!(t, psi)
    fi = get_field_eval(t.equilibrium, "B")
    B_cond = zeros(Float64, np, 3)
    @inbounds for i in 1:np
        pt_mask[i] || continue
        B_cond[i, :] = fi(t.r[i, 1:2])
    end
    set_psi!(t, psi_save)
    if cell_centered
        Bv = zeros(Float64, t.nc, 3)
        rcc = zeros(Float64, t.nc, 3)
        @inbounds for c in 1:t.nc
            i1 = t.lc[c, 1] + 1; i2 = t.lc[c, 2] + 1; i3 = t.lc[c, 3] + 1
            Bv[c, :] = (B_cond[i1, :] .+ B_cond[i2, :] .+ B_cond[i3, :]) ./ 3.0
            rcc[c, :] = (t.r[i1, :] .+ t.r[i2, :] .+ t.r[i3, :]) ./ 3.0
        end
        return (J_cond=J_cond, B_cond=Bv, mask=mask, R=rcc)
    end
    return (J_cond=J_cond, B_cond=B_cond, mask=mask, R=t.r)
end
export compute_forces_components

export OFTEnv, oftpy_set_debug, oftpy_set_nthreads
export TokamakerSettings, tokamaker_default_settings
export TokamakerReconSettings, tokamaker_recon_default_settings
export Tokamaker, TokaMakerEquilibrium
export setup_mesh!, setup_regions!, setup!, init_psi!, solve!, vac_solve!, reset!,
       copy_eq, replace_eq!,
       update_settings!, set_profiles!, set_targets!, get_targets,
       set_isoflux_constraints!, set_psi_constraints!, set_saddle_constraints!,
       set_coil_currents!, set_coil_bounds!, set_coil_vsc!, set_coil_reg!,
       set_resistivity!, set_coil_current_dist!,
       abspsi_to_normalized, psinorm_to_absolute,
       get_coil_currents, get_psi, set_psi!,
       compute_area_integral, compute_flux_integral,
       get_conductor_currents, get_conductor_source,
       coil_dict2vec, coil_vec2dict
export create_isoflux, create_power_flux_fun, create_spline_flux_fun,
       eval_green, read_eqdsk, read_ifile
export GEQDSKEquilibrium, read_geqdsk, cocos_params, cocosify, cocosify!,
       flip_Bt_Ip, flip_Bt_Ip!, save_geqdsk, eqdsk_to_bytes, eqdsk_from_bytes,
       eqdsk_from_raw

# Re-export the constraint helpers from EquilibriumModule
import .EquilibriumModule: get_globals, get_refs, get_xpoints, get_profiles,
                            get_q, calc_sauter_fc, calc_loopvoltage,
                            calc_delstar_curr, calc_jtor_plasma,
                            save_eqdsk, save_ifile, save_mug, save_tokamaker, load_refs!,
                            get_stats, calc_inductance, trace_surf
import .CoreModule: ffp_scale, set_ffp_scale!, p_scale, set_p_scale!,
                    o_point, lim_point, psi_bounds, diverted

export get_globals, get_refs, get_xpoints, get_profiles, get_q, calc_sauter_fc,
       calc_loopvoltage, calc_delstar_curr, calc_jtor_plasma,
       save_eqdsk, save_ifile, save_mug, save_tokamaker, get_stats, calc_inductance, trace_surf
export ffp_scale, set_ffp_scale!, p_scale, set_p_scale!,
       o_point, lim_point, psi_bounds, diverted
export TokamakerFieldInterpolator, get_field_eval
export setup_td!, step_td!, set_psi_dt!, eig_td, eig_wall
export MirnovCon, IpCon, FluxLoopCon, DFluxCon, PressCon, QCon, SaddleCon,
       ReconConstraints, run_reconstruction!, write_constraints_file
export Hmode_profiles, parameterize_edge_jBS, calculate_ln_lambda,
       redl_bootstrap, get_jphi_from_GS, find_peaks, peak_widths,
       analyze_bootstrap_edge_spike, find_optimal_scale, solve_jphi,
       solve_with_bootstrap
export save_gs_mesh, load_gs_mesh, GsDomain, define_region!,
       add_polygon!, add_rectangle!, add_annulus!, add_enclosed!,
       build_mesh!, get_coils, get_conductors, save_json

# Plotting functions implemented in the Makie package extension
# (`ext/TokaMakerMakieExt.jl`). Stubs here so the names exist for
# `import` and `methodswith` queries before Makie is loaded.
function plot_machine end
function plot_mesh end
function plot_topology end
function plot_psi end
function plot_constraints end
function plot_eddy end
export plot_machine, plot_mesh, plot_topology, plot_psi, plot_constraints, plot_eddy

end # module

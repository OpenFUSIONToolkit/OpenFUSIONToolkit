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
                   analyze_bootstrap_edge_spike, find_optimal_scale, solve_jphi

# Wire Bootstrap's runtime callbacks now that Core/Equilibrium are loaded.
function __init__()
    Bootstrap.set_targets!_anon[]   = (gs, Ip, pax) -> CoreModule.set_targets!(gs; Ip=Ip, pax=pax)
    Bootstrap.set_profiles!_anon[]  = (gs, ffp, pp) -> CoreModule.set_profiles!(gs; ffp_prof=ffp, pp_prof=pp)
    Bootstrap.solve!_anon[]         = (gs)         -> CoreModule.solve!(gs)
    Bootstrap.get_profiles_anon[]   = (eq; npsi, psi_pad) -> EquilibriumModule.get_profiles(eq; npsi=npsi, psi_pad=psi_pad)
    Bootstrap.get_q_anon[]          = (eq; npsi, psi_pad) -> EquilibriumModule.get_q(eq; npsi=npsi, psi_pad=psi_pad, compute_geo=true)
    Bootstrap.get_stats_anon[]      = (eq; lcfs_pad) -> EquilibriumModule.get_stats(eq; lcfs_pad=lcfs_pad)
    return nothing
end

include("Meshing.jl")
using .Meshing: save_gs_mesh, load_gs_mesh, GsDomain, define_region!,
                 add_polygon!, add_rectangle!, add_annulus!, add_enclosed!,
                 build_mesh!, get_coils, get_conductors

export OFTEnv, oftpy_set_debug, oftpy_set_nthreads
export TokamakerSettings, tokamaker_default_settings
export TokamakerReconSettings, tokamaker_recon_default_settings
export Tokamaker, TokaMakerEquilibrium
export setup_mesh!, setup_regions!, setup!, init_psi!, solve!, reset!,
       update_settings!, set_profiles!, set_targets!,
       set_isoflux_constraints!, set_psi_constraints!, set_saddle_constraints!,
       set_coil_currents!, set_coil_bounds!, set_coil_vsc!, set_coil_reg!,
       get_coil_currents, get_psi, set_psi!,
       coil_dict2vec, coil_vec2dict
export create_isoflux, create_power_flux_fun, create_spline_flux_fun,
       eval_green, read_eqdsk, read_ifile

# Re-export the constraint helpers from EquilibriumModule
import .EquilibriumModule: get_globals, get_refs, get_xpoints, get_profiles,
                            get_q, calc_sauter_fc, calc_loopvoltage,
                            calc_delstar_curr, calc_jtor_plasma,
                            save_eqdsk, save_ifile, save_mug, load_refs!,
                            get_stats, calc_inductance, trace_surf
import .CoreModule: ffp_scale, set_ffp_scale!, p_scale, set_p_scale!,
                    o_point, lim_point, psi_bounds, diverted

export get_globals, get_refs, get_xpoints, get_profiles, get_q, calc_sauter_fc,
       calc_loopvoltage, calc_delstar_curr, calc_jtor_plasma,
       save_eqdsk, save_ifile, save_mug, get_stats, calc_inductance, trace_surf
export ffp_scale, set_ffp_scale!, p_scale, set_p_scale!,
       o_point, lim_point, psi_bounds, diverted
export TokamakerFieldInterpolator, get_field_eval
export setup_td!, step_td!, set_psi_dt!, eig_td, eig_wall
export MirnovCon, IpCon, FluxLoopCon, DFluxCon, PressCon, QCon, SaddleCon,
       ReconConstraints, run_reconstruction!, write_constraints_file
export Hmode_profiles, parameterize_edge_jBS, calculate_ln_lambda,
       redl_bootstrap, get_jphi_from_GS, find_peaks, peak_widths,
       analyze_bootstrap_edge_spike, find_optimal_scale, solve_jphi
export save_gs_mesh, load_gs_mesh, GsDomain, define_region!,
       add_polygon!, add_rectangle!, add_annulus!, add_enclosed!,
       build_mesh!, get_coils, get_conductors

# Plotting functions implemented in the Makie package extension
# (`ext/TokaMakerMakieExt.jl`). Stubs here so the names exist for
# `import` and `methodswith` queries before Makie is loaded.
function plot_machine end
function plot_psi end
function plot_constraints end
function plot_eddy end
export plot_machine, plot_psi, plot_constraints, plot_eddy

end # module

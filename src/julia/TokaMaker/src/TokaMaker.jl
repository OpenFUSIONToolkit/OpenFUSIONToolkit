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
                   redl_bootstrap, get_jphi_from_GS

include("Meshing.jl")
using .Meshing: save_gs_mesh, load_gs_mesh, GsDomain

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
       ReconConstraints, run_reconstruction!
export Hmode_profiles, parameterize_edge_jBS, calculate_ln_lambda,
       redl_bootstrap, get_jphi_from_GS
export save_gs_mesh, load_gs_mesh, GsDomain

end # module

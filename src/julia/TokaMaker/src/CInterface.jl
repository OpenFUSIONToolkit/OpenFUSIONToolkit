module CInterface

# Raw 1:1 ccall wrappers for the BIND(C) entry points in
# src/python/wrappers/tokamaker_f.F90 plus the shared OFT entry points used
# by TokaMaker. Naming and argument order match _interface.py exactly so
# drift is easy to detect.

using ..LibPath: liboftpy, OFT_ERROR_SLEN, OFT_PATH_SLEN
using ..SettingsModule: TokamakerSettings, TokamakerReconSettings

export errbuf, check_err, padpath,
       c_oft_setup_smesh, c_oft_setup_vmesh, c_oftpy_load_xml,
       c_tokamaker_alloc, c_tokamaker_destroy,
       c_tokamaker_equil_copy, c_tokamaker_equil_set, c_tokamaker_equil_destroy,
       c_tokamaker_setup_regions, c_tokamaker_setup,
       c_tokamaker_load_profiles, c_tokamaker_init_psi,
       c_tokamaker_solve, c_tokamaker_vac_solve,
       c_tokamaker_setup_td, c_tokamaker_eig_td, c_tokamaker_eig_wall,
       c_tokamaker_step_td,
       c_tokamaker_get_mesh, c_tokamaker_get_limiter,
       c_tokamaker_get_psi, c_tokamaker_get_dels_curr, c_tokamaker_get_jtor,
       c_tokamaker_area_int, c_tokamaker_flux_int,
       c_tokamaker_get_coil_currents, c_tokamaker_get_plasma_Lmat,
       c_tokamaker_get_refs, c_tokamaker_trace_surf,
       c_tokamaker_get_q, c_tokamaker_sauter_fc,
       c_tokamaker_get_globals, c_tokamaker_gs_calc_vloop,
       c_tokamaker_get_profs, c_tokamaker_get_vfixed,
       c_tokamaker_get_field_eval, c_tokamaker_apply_field_eval,
       c_tokamaker_set_psi, c_tokamaker_set_psi_dt,
       c_tokamaker_set_settings, c_tokamaker_set_dipole_a,
       c_tokamaker_set_mirror_slosh, c_tokamaker_set_targets,
       c_tokamaker_set_isoflux, c_tokamaker_set_flux, c_tokamaker_set_saddles,
       c_tokamaker_set_coil_currents, c_tokamaker_set_coil_regmat,
       c_tokamaker_set_coil_bounds, c_tokamaker_set_coil_vsc, c_tokamaker_set_vcoil,
       c_tokamaker_save_eqdsk, c_tokamaker_save_ifile, c_tokamaker_save_mug,
       c_tokamaker_set_coil_current_dist

# ---------------------------------------------------------------------------
# Error handling helpers

errbuf() = zeros(UInt8, OFT_ERROR_SLEN)

function check_err(buf::Vector{UInt8}, where::AbstractString)
    if buf[1] != 0
        msg = unsafe_string(pointer(buf))
        error("$where: $msg")
    end
    return nothing
end

"Pad a path string into a NUL-terminated UInt8 vector of length OFT_PATH_SLEN."
function padpath(s::AbstractString, len::Int=OFT_PATH_SLEN)
    bytes = codeunits(s)
    length(bytes) >= len && error("Path string too long ($(length(bytes)) >= $len): $s")
    out = zeros(UInt8, len)
    @inbounds for i in eachindex(bytes)
        out[i] = bytes[i]
    end
    return out
end

# ---------------------------------------------------------------------------
# Shared OFT runtime entry points used by TokaMaker

function c_oft_setup_smesh(ndim, np, r_loc::Matrix{Float64}, npc, nc,
                           lc_loc::Matrix{Int32}, reg_loc::Vector{Int32},
                           order_out::Ref{Int32}, mesh_ptr::Ref{Ptr{Cvoid}})
    ccall((:oft_setup_smesh, liboftpy[]), Cvoid,
          (Cint, Cint, Ptr{Float64}, Cint, Cint, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Ptr{Cvoid}}),
          Int32(ndim), Int32(np), r_loc, Int32(npc), Int32(nc),
          lc_loc, reg_loc, order_out, mesh_ptr)
end

function c_oft_setup_vmesh(ndim, np, r_loc::Matrix{Float64}, npc, nc,
                           lc_loc::Matrix{Int32}, reg_loc::Vector{Int32},
                           order_out::Ref{Int32}, mesh_ptr::Ref{Ptr{Cvoid}})
    ccall((:oft_setup_vmesh, liboftpy[]), Cvoid,
          (Cint, Cint, Ptr{Float64}, Cint, Cint, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Ptr{Cvoid}}),
          Int32(ndim), Int32(np), r_loc, Int32(npc), Int32(nc),
          lc_loc, reg_loc, order_out, mesh_ptr)
end

function c_oftpy_load_xml(xml_file::Vector{UInt8}, oft_node_ptr::Ref{Ptr{Cvoid}})
    ccall((:oftpy_load_xml, liboftpy[]), Cvoid,
          (Ptr{UInt8}, Ptr{Ptr{Cvoid}}), xml_file, oft_node_ptr)
end

c_dump_cov() = ccall((:dump_cov, liboftpy[]), Cvoid, ())

# ---------------------------------------------------------------------------
# TokaMaker setup / lifecycle

function c_tokamaker_alloc(tMaker_ptr::Ref{Ptr{Cvoid}}, mesh_ptr::Ptr{Cvoid},
                           error_str::Vector{UInt8})
    ccall((:tokamaker_alloc, liboftpy[]), Cvoid,
          (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{UInt8}),
          tMaker_ptr, mesh_ptr, error_str)
end

c_tokamaker_destroy(tMaker_ptr::Ptr{Cvoid}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_destroy, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}), tMaker_ptr, error_str)

function c_tokamaker_equil_copy(tMaker_ptr::Ptr{Cvoid}, old_equil::Ptr{Cvoid},
                                new_equil::Ref{Ptr{Cvoid}}, error_str::Vector{UInt8})
    ccall((:tokamaker_equil_copy, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{UInt8}),
          tMaker_ptr, old_equil, new_equil, error_str)
end

c_tokamaker_equil_set(tMaker_ptr::Ptr{Cvoid}, new_equil::Ptr{Cvoid},
                      error_str::Vector{UInt8}) =
    ccall((:tokamaker_equil_set, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}),
          tMaker_ptr, new_equil, error_str)

c_tokamaker_equil_destroy(eq_ptr::Ptr{Cvoid}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_equil_destroy, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}), eq_ptr, error_str)

function c_tokamaker_setup_regions(tMaker_ptr::Ptr{Cvoid}, coil_file::Vector{UInt8},
                                   reg_eta::Vector{Float64}, contig_flag::Vector{Int32},
                                   xpoint_mask::Vector{Int32}, coil_nturns::Matrix{Float64},
                                   ncoils::Integer, error_str::Vector{UInt8})
    ccall((:tokamaker_setup_regions, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
           Ptr{Float64}, Cint, Ptr{UInt8}),
          tMaker_ptr, coil_file, reg_eta, contig_flag, xpoint_mask,
          coil_nturns, Int32(ncoils), error_str)
end

function c_tokamaker_setup(tMaker_ptr::Ptr{Cvoid}, order::Integer, full_domain::Bool,
                           ncoils::Ref{Int32}, coil_Lmat::Ref{Ptr{Float64}},
                           error_str::Vector{UInt8})
    ccall((:tokamaker_setup, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Cuchar, Ptr{Int32}, Ptr{Ptr{Float64}}, Ptr{UInt8}),
          tMaker_ptr, Int32(order), UInt8(full_domain), ncoils, coil_Lmat, error_str)
end

# ---------------------------------------------------------------------------
# Init / profiles / solve

c_tokamaker_load_profiles(eq_ptr::Ptr{Cvoid}, f_file::Vector{UInt8}, f_offset::Float64,
                          p_file::Vector{UInt8}, eta_file::Vector{UInt8},
                          f_NI_file::Vector{UInt8}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_load_profiles, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}, Cdouble, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}),
          eq_ptr, f_file, f_offset, p_file, eta_file, f_NI_file, error_str)

c_tokamaker_init_psi(tMaker_ptr::Ptr{Cvoid}, r0, z0, a, kappa, delta,
                     rhs_source::Ref{Float64}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_init_psi, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
           Ptr{Cdouble}, Ptr{UInt8}),
          tMaker_ptr, Float64(r0), Float64(z0), Float64(a), Float64(kappa),
          Float64(delta), rhs_source, error_str)

c_tokamaker_solve(tMaker_ptr::Ptr{Cvoid}, vacuum::Bool, error_str::Vector{UInt8}) =
    ccall((:tokamaker_solve, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cuchar, Ptr{UInt8}),
          tMaker_ptr, UInt8(vacuum), error_str)

c_tokamaker_vac_solve(tMaker_ptr::Ptr{Cvoid}, psi_in::Vector{Float64},
                      rhs_source::Ptr{Float64}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_vac_solve, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, psi_in, rhs_source, error_str)

# ---------------------------------------------------------------------------
# Time-dependent

c_tokamaker_setup_td(tMaker_ptr::Ptr{Cvoid}, dt, lin_tol, nl_tol, pre_plasma::Bool,
                    error_str::Vector{UInt8}) =
    ccall((:tokamaker_setup_td, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cuchar, Ptr{UInt8}),
          tMaker_ptr, Float64(dt), Float64(lin_tol), Float64(nl_tol),
          UInt8(pre_plasma), error_str)

c_tokamaker_eig_td(tMaker_ptr::Ptr{Cvoid}, omega, neigs::Integer,
                   eigs::Matrix{Float64}, eig_vecs::Matrix{Float64},
                   include_bounds::Bool, eta_plasma, pm::Bool,
                   error_str::Vector{UInt8}) =
    ccall((:tokamaker_eig_td, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Cint, Ptr{Float64}, Ptr{Float64}, Cuchar,
           Cdouble, Cuchar, Ptr{UInt8}),
          tMaker_ptr, Float64(omega), Int32(neigs), eigs, eig_vecs,
          UInt8(include_bounds), Float64(eta_plasma), UInt8(pm), error_str)

c_tokamaker_eig_wall(tMaker_ptr::Ptr{Cvoid}, neigs::Integer,
                    eigs::Matrix{Float64}, eig_vecs::Matrix{Float64},
                    pm::Bool, error_str::Vector{UInt8}) =
    ccall((:tokamaker_eig_wall, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Float64}, Cuchar, Ptr{UInt8}),
          tMaker_ptr, Int32(neigs), eigs, eig_vecs, UInt8(pm), error_str)

c_tokamaker_step_td(tMaker_ptr::Ptr{Cvoid}, curr::Vector{Float64},
                   volt::Vector{Float64}, time::Ref{Float64}, dt::Ref{Float64},
                   nl_its::Ref{Int32}, lin_its::Ref{Int32}, nretry::Ref{Int32},
                   error_str::Vector{UInt8}) =
    ccall((:tokamaker_step_td, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{UInt8}),
          tMaker_ptr, curr, volt, time, dt, nl_its, lin_its, nretry, error_str)

# ---------------------------------------------------------------------------
# Mesh / get accessors

c_tokamaker_get_mesh(tMaker_ptr::Ptr{Cvoid}, np::Ref{Int32},
                    r_loc::Ref{Ptr{Float64}}, nc::Ref{Int32},
                    lc_loc::Ref{Ptr{Int32}}, reg_loc::Ref{Ptr{Int32}},
                    error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_mesh, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Int32}, Ptr{Ptr{Float64}}, Ptr{Int32},
           Ptr{Ptr{Int32}}, Ptr{Ptr{Int32}}, Ptr{UInt8}),
          tMaker_ptr, np, r_loc, nc, lc_loc, reg_loc, error_str)

c_tokamaker_get_limiter(tMaker_ptr::Ptr{Cvoid}, np::Ref{Int32},
                       r_loc::Ref{Ptr{Float64}}, nloops::Ref{Int32},
                       loop_ptr::Ref{Ptr{Int32}}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_limiter, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Int32}, Ptr{Ptr{Float64}}, Ptr{Int32},
           Ptr{Ptr{Int32}}, Ptr{UInt8}),
          tMaker_ptr, np, r_loc, nloops, loop_ptr, error_str)

c_tokamaker_get_psi(eq_ptr::Ptr{Cvoid}, psi_vals::Vector{Float64},
                   psi_lim::Ref{Float64}, psi_max::Ref{Float64},
                   error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_psi, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, psi_vals, psi_lim, psi_max, error_str)

c_tokamaker_get_dels_curr(eq_ptr::Ptr{Cvoid}, psi_vals::Vector{Float64},
                         error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_dels_curr, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, psi_vals, error_str)

c_tokamaker_get_jtor(eq_ptr::Ptr{Cvoid}, jtor::Vector{Float64},
                    error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_jtor, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, jtor, error_str)

c_tokamaker_area_int(tMaker_ptr::Ptr{Cvoid}, vec_vals::Vector{Float64},
                    reg_ind::Integer, result::Ref{Float64},
                    error_str::Vector{UInt8}) =
    ccall((:tokamaker_area_int, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, vec_vals, Int32(reg_ind), result, error_str)

c_tokamaker_flux_int(eq_ptr::Ptr{Cvoid}, psi_vals::Vector{Float64},
                    field_vals::Vector{Float64}, nvals::Integer,
                    result::Ref{Float64}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_flux_int, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, psi_vals, field_vals, Int32(nvals), result, error_str)

c_tokamaker_get_coil_currents(eq_ptr::Ptr{Cvoid}, currents::Vector{Float64},
                             reg_currents::Vector{Float64},
                             error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_coil_currents, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, currents, reg_currents, error_str)

c_tokamaker_get_plasma_Lmat(eq_ptr::Ptr{Cvoid}, Lmat::Vector{Float64},
                           error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_plasma_Lmat, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, Lmat, error_str)

c_tokamaker_get_refs(eq_ptr::Ptr{Cvoid}, o_point::Ref{Ptr{Float64}},
                    lim_point::Ref{Ptr{Float64}}, x_points::Ref{Ptr{Float64}},
                    diverted::Ref{Ptr{UInt8}}, plasma_bounds::Ref{Ptr{Float64}},
                    ffp_scale::Ref{Ptr{Float64}}, p_scale::Ref{Ptr{Float64}},
                    has_plasma::Ref{Ptr{UInt8}}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_refs, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Ptr{Float64}}, Ptr{Ptr{Float64}}, Ptr{Ptr{Float64}},
           Ptr{Ptr{UInt8}}, Ptr{Ptr{Float64}}, Ptr{Ptr{Float64}}, Ptr{Ptr{Float64}},
           Ptr{Ptr{UInt8}}, Ptr{UInt8}),
          eq_ptr, o_point, lim_point, x_points, diverted, plasma_bounds,
          ffp_scale, p_scale, has_plasma, error_str)

c_tokamaker_trace_surf(tMaker_ptr::Ptr{Cvoid}, psi_surf::Float64,
                      points::Ref{Ptr{Float64}}, npoints::Ref{Int32},
                      error_str::Vector{UInt8}) =
    ccall((:tokamaker_trace_surf, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Ptr{Ptr{Float64}}, Ptr{Int32}, Ptr{UInt8}),
          tMaker_ptr, psi_surf, points, npoints, error_str)

c_tokamaker_get_q(eq_ptr::Ptr{Cvoid}, npsi::Integer, psi_q::Vector{Float64},
                 qvals::Vector{Float64}, ravgs::Matrix{Float64},
                 dl::Ref{Float64}, rbounds::Matrix{Float64}, zbounds::Matrix{Float64},
                 error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_q, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, Int32(npsi), psi_q, qvals, ravgs, dl, rbounds, zbounds, error_str)

c_tokamaker_sauter_fc(eq_ptr::Ptr{Cvoid}, npsi::Integer, psi_saut::Vector{Float64},
                    fc::Vector{Float64}, r_avgs::Matrix{Float64},
                    modb_avgs::Matrix{Float64}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_sauter_fc, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, Int32(npsi), psi_saut, fc, r_avgs, modb_avgs, error_str)

c_tokamaker_get_globals(eq_ptr::Ptr{Cvoid}, Itor::Ref{Float64},
                       centroid::Vector{Float64}, vol::Ref{Float64},
                       pvol::Ref{Float64}, dflux::Ref{Float64},
                       tflux::Ref{Float64}, bp_vol::Ref{Float64},
                       error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_globals, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, Itor, centroid, vol, pvol, dflux, tflux, bp_vol, error_str)

c_tokamaker_gs_calc_vloop(eq_ptr::Ptr{Cvoid}, vloop::Ref{Float64},
                         error_str::Vector{UInt8}) =
    ccall((:tokamaker_gs_calc_vloop, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, vloop, error_str)

c_tokamaker_get_profs(eq_ptr::Ptr{Cvoid}, npsi::Integer, psi_in::Vector{Float64},
                     f::Vector{Float64}, fp::Vector{Float64},
                     p::Vector{Float64}, pp::Vector{Float64},
                     error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_profs, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          eq_ptr, Int32(npsi), psi_in, f, fp, p, pp, error_str)

c_tokamaker_get_vfixed(eq_ptr::Ptr{Cvoid}, npts::Ref{Int32},
                      pts::Ref{Ptr{Float64}}, fluxes::Ref{Ptr{Float64}},
                      error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_vfixed, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Int32}, Ptr{Ptr{Float64}}, Ptr{Ptr{Float64}}, Ptr{UInt8}),
          eq_ptr, npts, pts, fluxes, error_str)

c_tokamaker_get_field_eval(eq_ptr::Ptr{Cvoid}, imode::Integer,
                          int_obj::Ref{Ptr{Cvoid}}, error_str::Vector{UInt8}) =
    ccall((:tokamaker_get_field_eval, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Ptr{Cvoid}}, Ptr{UInt8}),
          eq_ptr, Int32(imode), int_obj, error_str)

c_tokamaker_apply_field_eval(eq_ptr::Ptr{Cvoid}, int_obj::Ptr{Cvoid},
                            int_type::Integer, pt::Vector{Float64},
                            fbary_tol::Float64, cell::Ref{Int32},
                            dim::Integer, field::Vector{Float64}) =
    ccall((:tokamaker_apply_field_eval, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Float64}, Cdouble,
           Ptr{Int32}, Cint, Ptr{Float64}),
          eq_ptr, int_obj, Int32(int_type), pt, fbary_tol, cell, Int32(dim), field)

# ---------------------------------------------------------------------------
# Setters

c_tokamaker_set_psi(eq_ptr::Ptr{Cvoid}, psi_vals::Vector{Float64},
                   update_bounds::Bool, error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_psi, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Cuchar, Ptr{UInt8}),
          eq_ptr, psi_vals, UInt8(update_bounds), error_str)

c_tokamaker_set_psi_dt(tMaker_ptr::Ptr{Cvoid}, psi_vals::Vector{Float64},
                      icoils::Vector{Float64}, vcoils::Vector{Float64},
                      dt::Float64, error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_psi_dt, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cdouble, Ptr{UInt8}),
          tMaker_ptr, psi_vals, icoils, vcoils, dt, error_str)

c_tokamaker_set_settings(tMaker_ptr::Ptr{Cvoid}, settings::TokamakerSettings,
                        error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_settings, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ref{TokamakerSettings}, Ptr{UInt8}),
          tMaker_ptr, settings, error_str)

c_tokamaker_set_dipole_a(tMaker_ptr::Ptr{Cvoid}, dipole_a::Float64,
                        error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_dipole_a, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Ptr{UInt8}),
          tMaker_ptr, dipole_a, error_str)

c_tokamaker_set_mirror_slosh(tMaker_ptr::Ptr{Cvoid}, mirror_n::Float64,
                            mirror_bturn::Float64, mirror_zthroat::Float64,
                            error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_mirror_slosh, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Ptr{UInt8}),
          tMaker_ptr, mirror_n, mirror_bturn, mirror_zthroat, error_str)

c_tokamaker_set_targets(tMaker_ptr::Ptr{Cvoid}, ip::Float64, ip_ratio::Float64,
                       pax::Float64, estore::Float64, R0::Float64, V0::Float64,
                       error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_targets, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{UInt8}),
          tMaker_ptr, ip, ip_ratio, pax, estore, R0, V0, error_str)

c_tokamaker_set_isoflux(tMaker_ptr::Ptr{Cvoid}, targets::Matrix{Float64},
                       ref_points::Matrix{Float64}, weights::Vector{Float64},
                       ntargets::Integer, grad_wt_lim::Float64,
                       error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_isoflux, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cint, Cdouble, Ptr{UInt8}),
          tMaker_ptr, targets, ref_points, weights, Int32(ntargets), grad_wt_lim, error_str)

c_tokamaker_set_flux(tMaker_ptr::Ptr{Cvoid}, locations::Matrix{Float64},
                    targets::Vector{Float64}, weights::Vector{Float64},
                    ntargets::Integer, grad_wt_lim::Float64,
                    error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_flux, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cint, Cdouble, Ptr{UInt8}),
          tMaker_ptr, locations, targets, weights, Int32(ntargets), grad_wt_lim, error_str)

c_tokamaker_set_saddles(tMaker_ptr::Ptr{Cvoid}, targets::Matrix{Float64},
                       weights::Vector{Float64}, ntargets::Integer,
                       error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_saddles, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{UInt8}),
          tMaker_ptr, targets, weights, Int32(ntargets), error_str)

c_tokamaker_set_coil_currents(tMaker_ptr::Ptr{Cvoid}, currents::Vector{Float64},
                             error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_coil_currents, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, currents, error_str)

c_tokamaker_set_coil_regmat(tMaker_ptr::Ptr{Cvoid}, nregularize::Integer,
                           coil_reg_mat::Matrix{Float64},
                           coil_reg_targets::Vector{Float64},
                           coil_reg_weights::Vector{Float64},
                           error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_coil_regmat, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, Int32(nregularize), coil_reg_mat, coil_reg_targets,
          coil_reg_weights, error_str)

c_tokamaker_set_coil_bounds(tMaker_ptr::Ptr{Cvoid}, coil_bounds::Matrix{Float64},
                           error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_coil_bounds, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, coil_bounds, error_str)

c_tokamaker_set_coil_vsc(tMaker_ptr::Ptr{Cvoid}, coil_gains::Vector{Float64},
                        error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_coil_vsc, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, coil_gains, error_str)

c_tokamaker_set_vcoil(tMaker_ptr::Ptr{Cvoid}, rcoils::Vector{Float64},
                     error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_vcoil, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{UInt8}),
          tMaker_ptr, rcoils, error_str)

c_tokamaker_set_coil_current_dist(tMaker_ptr::Ptr{Cvoid}, iCoil::Integer,
                                 curr_dist::Vector{Float64},
                                 dist_pointer::Ref{Ptr{Float64}},
                                 normalize::Bool, error_str::Vector{UInt8}) =
    ccall((:tokamaker_set_coil_current_dist, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Ptr{Float64}}, Cuchar, Ptr{UInt8}),
          tMaker_ptr, Int32(iCoil), curr_dist, dist_pointer, UInt8(normalize), error_str)

# ---------------------------------------------------------------------------
# I/O

c_tokamaker_save_eqdsk(eq_ptr::Ptr{Cvoid}, filename::Vector{UInt8}, nr::Integer,
                      nz::Integer, rbounds::Vector{Float64}, zbounds::Vector{Float64},
                      run_info::Vector{UInt8}, psi_pad::Float64, rcentr::Float64,
                      trunc_eq::Bool, lim_filename::Vector{UInt8},
                      lcfs_press::Float64, cocos::Integer,
                      error_str::Vector{UInt8}) =
    ccall((:tokamaker_save_eqdsk, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Ptr{Float64},
           Ptr{UInt8}, Cdouble, Cdouble, Cuchar, Ptr{UInt8}, Cdouble, Cint, Ptr{UInt8}),
          eq_ptr, filename, Int32(nr), Int32(nz), rbounds, zbounds, run_info,
          psi_pad, rcentr, UInt8(trunc_eq), lim_filename, lcfs_press,
          Int32(cocos), error_str)

c_tokamaker_save_ifile(eq_ptr::Ptr{Cvoid}, filename::Vector{UInt8}, npsi::Integer,
                      ntheta::Integer, psi_pad::Float64, lcfs_press::Float64,
                      pack_lcfs::Bool, single_prec::Bool,
                      error_str::Vector{UInt8}) =
    ccall((:tokamaker_save_ifile, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Cdouble, Cdouble, Cuchar, Cuchar, Ptr{UInt8}),
          eq_ptr, filename, Int32(npsi), Int32(ntheta), psi_pad, lcfs_press,
          UInt8(pack_lcfs), UInt8(single_prec), error_str)

c_tokamaker_save_mug(eq_ptr::Ptr{Cvoid}, filename::Vector{UInt8},
                    error_str::Vector{UInt8}) =
    ccall((:tokamaker_save_mug, liboftpy[]), Cvoid,
          (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}),
          eq_ptr, filename, error_str)

end # module

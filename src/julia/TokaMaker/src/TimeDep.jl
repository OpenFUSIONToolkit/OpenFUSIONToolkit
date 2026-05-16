module TimeDep

using ..CInterface
using ..CoreModule: Tokamaker, coil_dict2vec, get_coil_currents

export setup_td!, step_td!, set_psi_dt!, eig_td, eig_wall

function setup_td!(t::Tokamaker, dt::Real, lin_tol::Real, nl_tol::Real;
                   pre_plasma::Bool=false)
    buf = errbuf()
    c_tokamaker_setup_td(t.tmaker_ptr, dt, lin_tol, nl_tol, pre_plasma, buf)
    check_err(buf, "setup_td")
    return t
end

function step_td!(t::Tokamaker, time::Real, dt::Real;
                  coil_currents::Union{Nothing,AbstractDict}=nothing,
                  coil_voltages::Union{Nothing,AbstractDict}=nothing)
    if coil_currents === nothing
        coil_currents, _ = get_coil_currents(t)
    end
    curr_arr = Vector{Float64}(coil_dict2vec(t, coil_currents; keep_virtual=true))
    volt_arr = if coil_voltages !== nothing
        Vector{Float64}(coil_dict2vec(t, coil_voltages; keep_virtual=true))
    else
        Vector{Float64}(coil_dict2vec(t, Dict{String,Float64}(); keep_virtual=true))
    end
    time_ref = Ref{Float64}(Float64(time))
    dt_ref = Ref{Float64}(Float64(dt))
    nl_its = Ref{Int32}(0); lin_its = Ref{Int32}(0); nretry = Ref{Int32}(0)
    buf = errbuf()
    c_tokamaker_step_td(t.tmaker_ptr, curr_arr, volt_arr, time_ref, dt_ref,
                        nl_its, lin_its, nretry, buf)
    check_err(buf, "step_td")
    return (time=time_ref[], dt=dt_ref[],
            nl_its=Int(nl_its[]), lin_its=Int(lin_its[]),
            nretry=Int(nretry[]))
end

function set_psi_dt!(t::Tokamaker, psi0::AbstractVector, dt::Real;
                     coil_currents::Union{Nothing,AbstractDict}=nothing,
                     coil_voltages::Union{Nothing,AbstractDict}=nothing)
    length(psi0) == t.np || error("psi0 length must equal np=$(t.np)")
    psi = Vector{Float64}(psi0)
    if coil_currents === nothing
        coil_currents, _ = get_coil_currents(t)
    end
    curr_arr = Vector{Float64}(coil_dict2vec(t, coil_currents; keep_virtual=true))
    volt_arr = if coil_voltages !== nothing
        Vector{Float64}(coil_dict2vec(t, coil_voltages; keep_virtual=true))
    else
        Vector{Float64}(coil_dict2vec(t, Dict{String,Float64}(); keep_virtual=true))
    end
    buf = errbuf()
    c_tokamaker_set_psi_dt(t.tmaker_ptr, psi, curr_arr, volt_arr, Float64(dt), buf)
    check_err(buf, "set_psi_dt")
    return t
end

function eig_td(t::Tokamaker; omega::Real=-1e4, neigs::Integer=4,
                include_bounds::Bool=true, eta_plasma::Real=-1.0,
                pm::Bool=false)
    n_psi = t.np
    eigs = zeros(Float64, 2, neigs)
    eig_vecs = zeros(Float64, n_psi, neigs)
    buf = errbuf()
    c_tokamaker_eig_td(t.tmaker_ptr, Float64(omega), neigs, eigs, eig_vecs,
                       include_bounds, Float64(eta_plasma), pm, buf)
    check_err(buf, "eig_td")
    return (eigs=permutedims(eigs), eig_vecs=eig_vecs)
end

function eig_wall(t::Tokamaker; neigs::Integer=4, pm::Bool=false)
    n_psi = t.np
    eigs = zeros(Float64, 2, neigs)
    eig_vecs = zeros(Float64, n_psi, neigs)
    buf = errbuf()
    c_tokamaker_eig_wall(t.tmaker_ptr, neigs, eigs, eig_vecs, pm, buf)
    check_err(buf, "eig_wall")
    return (eigs=permutedims(eigs), eig_vecs=eig_vecs)
end

end # module

module FieldEvalModule

using ..CInterface
using ..EquilibriumModule: TokaMakerEquilibrium

export TokamakerFieldInterpolator, get_field_eval

# Field type → (imode, dim_return)
# Mirrors the imode mapping in _core.py:get_field_eval.
const _FIELD_TYPES = Dict{String,Tuple{Int,Int}}(
    "B"   => (1, 2),
    "psi" => (2, 1),
    "F"   => (3, 1),
    "P"   => (4, 1),
    "dPSI"=> (5, 2),
    "dBr" => (6, 2),
    "dBt" => (7, 2),
    "dBz" => (8, 2),
)

mutable struct TokamakerFieldInterpolator
    eq::TokaMakerEquilibrium
    int_obj::Ptr{Cvoid}
    int_type::Int
    dim_return::Int
    dim_eval::Int
    cell::Base.RefValue{Int32}
    fbary_tol::Float64
    val::Vector{Float64}
    finalized::Bool

    function TokamakerFieldInterpolator(eq::TokaMakerEquilibrium, int_obj::Ptr{Cvoid},
                                         int_type::Integer, dim::Integer;
                                         fbary_tol::Real=1e-8)
        dim_eval = dim == 2 ? 3 : Int(dim)
        fi = new(eq, int_obj, Int(int_type), Int(dim), dim_eval,
                 Ref{Int32}(-1), Float64(fbary_tol), zeros(Float64, dim_eval), false)
        finalizer(_destroy_field_eval!, fi)
        return fi
    end
end

function _destroy_field_eval!(fi::TokamakerFieldInterpolator)
    fi.finalized && return
    if fi.int_obj != C_NULL
        try
            pt = zeros(Float64, 3)
            c_tokamaker_apply_field_eval(fi.eq.eq_ptr, fi.int_obj, -fi.int_type, pt,
                                          fi.fbary_tol, fi.cell, fi.dim_eval, fi.val)
        catch
        end
        fi.int_obj = C_NULL
    end
    fi.finalized = true
    return nothing
end

function get_field_eval(eq::TokaMakerEquilibrium, field_type::AbstractString)
    # Match Python's case-insensitive lookup
    canonical = Dict("B" => "B", "PSI" => "psi", "F" => "F", "P" => "P",
                     "DPSI" => "dPSI", "DBR" => "dBr", "DBT" => "dBt", "DBZ" => "dBz")
    key = get(canonical, uppercase(field_type), nothing)
    key === nothing && error("Unknown field_type \"$field_type\" (expected B, psi, F, P, dPSI, dBr, dBt, dBz)")
    imode, dim = _FIELD_TYPES[key]
    int_obj = Ref{Ptr{Cvoid}}(C_NULL)
    buf = errbuf()
    c_tokamaker_get_field_eval(eq.eq_ptr, imode, int_obj, buf)
    check_err(buf, "get_field_eval")
    return TokamakerFieldInterpolator(eq, int_obj[], imode, dim)
end

"""
    (fi::TokamakerFieldInterpolator)(pt) -> Vector{Float64}

Evaluate the interpolated field at point `pt = [r, z]`. Returns a vector of
length `dim_return`. Repeated calls reuse a cached starting cell guess for
fast lookup near previously-evaluated points.
"""
function (fi::TokamakerFieldInterpolator)(pt::AbstractVector)
    pt_eval = zeros(Float64, 3)
    pt_eval[1] = Float64(pt[1])
    pt_eval[2] = Float64(pt[2])
    c_tokamaker_apply_field_eval(fi.eq.eq_ptr, fi.int_obj, fi.int_type, pt_eval,
                                  fi.fbary_tol, fi.cell, fi.dim_eval, fi.val)
    return copy(fi.val[1:fi.dim_return])
end

end # module

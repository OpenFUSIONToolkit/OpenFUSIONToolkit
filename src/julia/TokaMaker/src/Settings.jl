module SettingsModule

export TokamakerSettings, tokamaker_default_settings, TokamakerReconSettings,
       tokamaker_recon_default_settings

# Layout matches C struct in tokamaker_f.F90:
#   bool pm, free_boundary, limited_only, dipole_mode, mirror_mode  (5 bytes + 3 pad)
#   int32 maxits, mode                                               (8 bytes)
#   double urf, nl_tol, rmin, lim_zmax                               (32 bytes)
#   char* limiter_file                                                (8 bytes)
# Total: 56 bytes on 64-bit.
#
# Julia mutable structs do not auto-pad mid-struct, so we pack the 5 bools as
# UInt8 and emit explicit 3-byte padding to match the C ABI. We expose property
# accessors that present the bools as `Bool` for the user.

mutable struct TokamakerSettings
    _pm::UInt8
    _free_boundary::UInt8
    _limited_only::UInt8
    _dipole_mode::UInt8
    _mirror_mode::UInt8
    _pad1::UInt8
    _pad2::UInt8
    _pad3::UInt8
    maxits::Int32
    mode::Int32
    urf::Float64
    nl_tol::Float64
    rmin::Float64
    lim_zmax::Float64
    limiter_file::Ptr{UInt8}
end

const _BOOL_FIELDS = (:pm, :free_boundary, :limited_only, :dipole_mode, :mirror_mode)
const _BOOL_BACKING = Dict(:pm => :_pm, :free_boundary => :_free_boundary,
                          :limited_only => :_limited_only, :dipole_mode => :_dipole_mode,
                          :mirror_mode => :_mirror_mode)

function Base.getproperty(s::TokamakerSettings, name::Symbol)
    if name in _BOOL_FIELDS
        return getfield(s, _BOOL_BACKING[name]) != 0
    end
    return getfield(s, name)
end

function Base.setproperty!(s::TokamakerSettings, name::Symbol, value)
    if name in _BOOL_FIELDS
        setfield!(s, _BOOL_BACKING[name], UInt8(value ? 1 : 0))
        return value
    end
    setfield!(s, name, convert(fieldtype(TokamakerSettings, name), value))
end

function Base.propertynames(::TokamakerSettings, private::Bool=false)
    public = (_BOOL_FIELDS..., :maxits, :mode, :urf, :nl_tol, :rmin, :lim_zmax, :limiter_file)
    private && return (fieldnames(TokamakerSettings)..., public...)
    return public
end

"""
Module-level pinned NUL-terminated `"none"` byte buffer. Settings structs need
`limiter_file` to point at valid memory; the default Python value is `path2c('none')`
which is a c_char_p with that string. We pin a single buffer here so the
returned `TokamakerSettings` doesn't dangle when its constructor's locals
fall out of scope.
"""
const _NONE_PATH = let v = Vector{UInt8}(undef, 5); v .= UInt8.(('n','o','n','e','\0')); v end

function tokamaker_default_settings()
    return TokamakerSettings(
        UInt8(1),    # pm
        UInt8(1),    # free_boundary
        UInt8(0),    # limited_only
        UInt8(0),    # dipole_mode
        UInt8(0),    # mirror_mode
        UInt8(0), UInt8(0), UInt8(0),  # padding
        Int32(40),   # maxits
        Int32(1),    # mode
        0.2,         # urf
        1.0e-6,      # nl_tol
        0.0,         # rmin
        1.0e99,      # lim_zmax
        pointer(_NONE_PATH),  # limiter_file -> "none"
    )
end

# Reconstruction settings struct (12 fields per reconstruction.py)
mutable struct TokamakerReconSettings
    _fitI::UInt8
    _fitP::UInt8
    _fit_Pscale::UInt8
    _fit_FFPscale::UInt8
    _fitR0::UInt8
    _fitV0::UInt8
    _fitCoils::UInt8
    _fitF0::UInt8
    _fixedCentering::UInt8
    _pm::UInt8
    _pad1::UInt8
    _pad2::UInt8
    _pad3::UInt8
    _pad4::UInt8
    _pad5::UInt8
    _pad6::UInt8
    infile::Ptr{UInt8}
    outfile::Ptr{UInt8}
end

const _RECON_BOOL_FIELDS = (:fitI, :fitP, :fit_Pscale, :fit_FFPscale, :fitR0,
                            :fitV0, :fitCoils, :fitF0, :fixedCentering, :pm)
const _RECON_BOOL_BACKING = Dict(s => Symbol("_" * String(s)) for s in _RECON_BOOL_FIELDS)

function Base.getproperty(s::TokamakerReconSettings, name::Symbol)
    if name in _RECON_BOOL_FIELDS
        return getfield(s, _RECON_BOOL_BACKING[name]) != 0
    end
    return getfield(s, name)
end

function Base.setproperty!(s::TokamakerReconSettings, name::Symbol, value)
    if name in _RECON_BOOL_FIELDS
        setfield!(s, _RECON_BOOL_BACKING[name], UInt8(value ? 1 : 0))
        return value
    end
    setfield!(s, name, convert(fieldtype(TokamakerReconSettings, name), value))
end

function tokamaker_recon_default_settings()
    return TokamakerReconSettings(
        UInt8(1),  # fitI
        UInt8(1),  # fitP
        UInt8(0),  # fit_Pscale
        UInt8(0),  # fit_FFPscale
        UInt8(0),  # fitR0
        UInt8(0),  # fitV0
        UInt8(0),  # fitCoils
        UInt8(0),  # fitF0
        UInt8(0),  # fixedCentering
        UInt8(1),  # pm
        UInt8(0), UInt8(0), UInt8(0), UInt8(0), UInt8(0), UInt8(0),  # padding
        C_NULL,    # infile
        C_NULL,    # outfile
    )
end

end # module

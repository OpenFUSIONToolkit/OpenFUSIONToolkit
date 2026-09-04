module Profiles

# Profile dict handling, mirrors `create_prof_file` in _core.py.
# Profile dicts have:
#   "type"  ::String  -> "flat" | "linterp" | "jphi-linterp"
#   "x"     ::Vector  -> normalized psi knots (only for linterp variants)
#   "y"     ::Vector  -> values at knots

export write_profile_file

function write_profile_file(filename::AbstractString, prof::AbstractDict, name::AbstractString;
                           psi_convention::Int=0)
    ptype = String(prof["type"])
    lines = String[ptype]
    if ptype == "flat"
        # nothing else
    elseif ptype == "linterp" || ptype == "jphi-linterp"
        haskey(prof, "x") || error("Profile $name missing key 'x'")
        haskey(prof, "y") || error("Profile $name missing key 'y'")
        x = collect(Float64.(prof["x"]))
        y = collect(Float64.(prof["y"]))
        any(diff(x) .< 0.0) && error("psi values in $name profile must be monotonically increasing")
        (x[1] < 0.0 || x[end] > 1.0) && error("Invalid psi values in $name profile ($(x[1]), $(x[end]))")
        if psi_convention == 0
            x = 1.0 .- x
            order = sortperm(x)
            x = x[order]
            y = y[order]
        elseif psi_convention != 1
            error("Unknown psi_convention $psi_convention; must be 0 (tokamak) or 1 (spheromak)")
        end
        push!(lines, "$(length(x)-1) $(y[1])")
        push!(lines, join(string.(x[2:end]), " "))
        push!(lines, join(string.(y[2:end]), " "))
    else
        error("Invalid profile type \"$ptype\" (expected flat, linterp, jphi-linterp)")
    end
    open(filename, "w") do io
        write(io, join(lines, "\n"))
    end
    return filename
end

end # module

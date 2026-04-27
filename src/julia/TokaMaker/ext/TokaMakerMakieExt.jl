module TokaMakerMakieExt

# Makie plotting extension for TokaMaker.jl.
# Mirrors the Python `plot_machine`, `plot_psi`, and `plot_constraints` methods
# on `OpenFUSIONToolkit.TokaMaker._core.TokaMaker`.

using TokaMaker
using TokaMaker: Tokamaker, TokaMakerEquilibrium
import TokaMaker.EquilibriumModule
using Makie
using Makie: Point2f, GeometryBasics

# ---------------------------------------------------------------------------
# Common helpers

function _rzplot(gs::Tokamaker)
    if gs.settings.mirror_mode
        return Float64.(gs.r[:, 2]), Float64.(gs.r[:, 1])
    end
    return Float64.(gs.r[:, 1]), Float64.(gs.r[:, 2])
end

_faces1(gs::Tokamaker) = Int.(gs.lc) .+ 1   # 0-based mesh storage -> 1-based

# Shade a subset of triangles in a single solid color (Python tricontourf w/ const).
function _shade_solid!(ax, r::AbstractVector, z::AbstractVector,
                       faces::AbstractMatrix{<:Integer},
                       mask::AbstractVector{Bool}, color)
    color === nothing && return nothing
    any(mask) || return nothing
    sub = faces[mask, :]
    pts = [Point2f(r[i], z[i]) for i in eachindex(r)]
    tris = [GeometryBasics.TriangleFace{Int}(sub[k, 1], sub[k, 2], sub[k, 3])
            for k in 1:size(sub, 1)]
    Makie.mesh!(ax, pts, tris; color=color, shading=Makie.NoShading)
    return nothing
end

# Per-face colormap shading (analog of Python tripcolor) by expanding each
# triangle to its own three vertices so vertex-color = face-color.
function _tripcolor!(ax, r::AbstractVector, z::AbstractVector,
                     faces::AbstractMatrix{<:Integer},
                     mask::AbstractVector{Bool}, vals::AbstractVector{<:Real};
                     colormap, colorrange=nothing)
    any(mask) || return nothing
    sub = faces[mask, :]
    sub_vals = vals[mask]
    n = size(sub, 1)
    pts = Vector{Point2f}(undef, 3n)
    cols = Vector{Float64}(undef, 3n)
    for k in 1:n
        for j in 1:3
            i = sub[k, j]
            pts[3(k - 1) + j] = Point2f(r[i], z[i])
            cols[3(k - 1) + j] = sub_vals[k]
        end
    end
    tris = [GeometryBasics.TriangleFace{Int}(3(k - 1) + 1, 3(k - 1) + 2, 3(k - 1) + 3)
            for k in 1:n]
    if colorrange === nothing
        return Makie.mesh!(ax, pts, tris; color=cols, colormap=colormap,
                           shading=Makie.NoShading)
    end
    return Makie.mesh!(ax, pts, tris; color=cols, colormap=colormap,
                       colorrange=colorrange, shading=Makie.NoShading)
end

# Marching triangles: emit pairs of Point2f per triangle that crosses `level`.
function _contour_segments(r::AbstractVector, z::AbstractVector,
                           faces::AbstractMatrix{<:Integer},
                           psi::AbstractVector, level::Real;
                           cell_mask::Union{Nothing,AbstractVector{Bool}}=nothing)
    segs = Point2f[]
    @inbounds for k in 1:size(faces, 1)
        cell_mask === nothing || cell_mask[k] || continue
        i1, i2, i3 = faces[k, 1], faces[k, 2], faces[k, 3]
        v = (psi[i1] - level, psi[i2] - level, psi[i3] - level)
        crossings = NTuple{2,Float64}[]
        for (a, b) in ((1, 2), (2, 3), (3, 1))
            va, vb = v[a], v[b]
            if (va <= 0.0 && vb > 0.0) || (va > 0.0 && vb <= 0.0)
                ia = (i1, i2, i3)[a]
                ib = (i1, i2, i3)[b]
                t = va / (va - vb)
                push!(crossings,
                      (r[ia] + t * (r[ib] - r[ia]),
                       z[ia] + t * (z[ib] - z[ia])))
            end
        end
        if length(crossings) >= 2
            push!(segs, Point2f(crossings[1][1], crossings[1][2]))
            push!(segs, Point2f(crossings[2][1], crossings[2][2]))
        end
    end
    return segs
end

function _level_color(color, colormap, k::Int, n::Int)
    color !== nothing && return color
    colormap === nothing && return nothing
    grad = Makie.cgrad(colormap)
    t = n <= 1 ? 0.5 : (k - 1) / (n - 1)
    return grad[t]
end

function _draw_levels!(ax, r, z, faces, psi, levels, color, colormap,
                       linestyle, linewidth;
                       cell_mask::Union{Nothing,AbstractVector{Bool}}=nothing)
    isempty(levels) && return nothing
    n = length(levels)
    for (k, c) in enumerate(levels)
        col = _level_color(color, colormap, k, n)
        col === nothing && continue
        segs = _contour_segments(r, z, faces, psi, c; cell_mask=cell_mask)
        isempty(segs) && continue
        kw = (color=col, linewidth=linewidth)
        if linestyle !== nothing
            kw = merge(kw, (linestyle=linestyle,))
        end
        Makie.linesegments!(ax, segs; kw...)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# plot_machine

const JULIA_BLUE   = Makie.RGBf(0.251, 0.388, 0.847)   # #4063D8
const JULIA_RED    = Makie.RGBf(0.796, 0.235, 0.20)    # #CB3C33
const JULIA_GREEN  = Makie.RGBf(0.220, 0.596, 0.149)   # #389826
const JULIA_PURPLE = Makie.RGBf(0.584, 0.345, 0.698)   # #9558B2

function TokaMaker.plot_machine(fig, ax::Makie.AbstractAxis, gs::Tokamaker;
                                equilibrium::Union{Nothing,TokaMakerEquilibrium}=nothing,
                                vacuum_color=Makie.RGBf(0.96, 0.96, 0.97),
                                cond_color=Makie.RGBf(0.45, 0.50, 0.58),
                                limiter_color=:black,
                                coil_color=JULIA_PURPLE,
                                coil_colormap=nothing,
                                coil_symmap::Bool=false,
                                coil_scale::Real=1.0,
                                coil_clabel::Union{Nothing,AbstractString}=raw"$I_C$ [A]",
                                colorbar=nothing)
    r, z = _rzplot(gs)
    faces = _faces1(gs)

    # Vacuum shading (regions 2..nvac+1)
    if vacuum_color !== nothing
        mask = Vector((gs.reg .> 1) .& (gs.reg .<= gs.nvac + 1))
        _shade_solid!(ax, r, z, faces, mask, vacuum_color)
    end

    # Coil shading
    if coil_colormap !== nothing
        eq = equilibrium === nothing ? gs.equilibrium : equilibrium
        eq === nothing && error("coil_colormap requires a solved equilibrium")
        _, region_currents = TokaMaker.get_coil_currents(eq)
        mesh_currents = zeros(Float64, gs.nc)
        if gs.ncoils > 0
            for k in 1:gs.nc
                rid = Int(gs.reg[k])
                if rid >= 1 && rid <= length(region_currents)
                    mesh_currents[k] = region_currents[rid]
                end
            end
        end
        mask = Vector(abs.(mesh_currents) .> 0.0)
        if any(mask)
            mesh_currents .*= coil_scale
            crange = nothing
            if coil_symmap
                m = maximum(abs.(mesh_currents))
                crange = (-m, m)
            end
            plt = _tripcolor!(ax, r, z, faces, mask, mesh_currents;
                              colormap=coil_colormap, colorrange=crange)
            if coil_clabel !== nothing && fig !== nothing && plt !== nothing
                Makie.Colorbar(fig[1, end + 1], plt; label=coil_clabel)
            end
        end
    elseif coil_color !== nothing
        for (_, info) in gs.coil_dict
            mask = Vector(gs.reg .== Int32(info["reg_id"]))
            _shade_solid!(ax, r, z, faces, mask, coil_color)
        end
    end

    # Conductor shading (excluding entries flagged as vacuum stand-ins)
    if cond_color !== nothing
        for (_, info) in gs.cond_dict
            haskey(info, "vac_id") && continue
            mask = Vector(gs.reg .== Int32(info["reg_id"]))
            _shade_solid!(ax, r, z, faces, mask, cond_color)
        end
    end

    # Limiter contour(s)
    if limiter_color !== nothing && !isempty(gs.lim_contours)
        for contour in gs.lim_contours
            if gs.settings.mirror_mode
                Makie.lines!(ax, contour[:, 2], contour[:, 1]; color=limiter_color)
            else
                Makie.lines!(ax, contour[:, 1], contour[:, 2]; color=limiter_color)
            end
        end
    end

    ax.aspect = Makie.DataAspect()
    return colorbar
end

TokaMaker.plot_machine(ax::Makie.AbstractAxis, gs::Tokamaker; kwargs...) =
    TokaMaker.plot_machine(nothing, ax, gs; kwargs...)

function TokaMaker.plot_machine(gs::Tokamaker; kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; xlabel="R (m)", ylabel="Z (m)")
    TokaMaker.plot_machine(fig, ax, gs; kwargs...)
    return fig, ax
end

# ---------------------------------------------------------------------------
# plot_psi

function TokaMaker.plot_psi(fig, ax::Makie.AbstractAxis, gs::Tokamaker;
                            equilibrium::Union{Nothing,TokaMakerEquilibrium}=nothing,
                            psi::Union{Nothing,AbstractVector}=nothing,
                            normalized::Bool=true,
                            plasma_color=nothing,
                            plasma_colormap=:plasma,
                            plasma_nlevels::Int=8,
                            plasma_levels::Union{Nothing,AbstractVector}=nothing,
                            plasma_linestyle=nothing,
                            plasma_linewidth::Real=1.5,
                            vacuum_color=Makie.RGBf(0.40, 0.40, 0.45),
                            vacuum_colormap=nothing,
                            vacuum_nlevels::Int=8,
                            vacuum_levels::Union{Nothing,AbstractVector}=nothing,
                            vacuum_linestyle=nothing,
                            vacuum_linewidth::Real=0.8,
                            xpoint_color=JULIA_RED,
                            xpoint_marker=:xcross,
                            xpoint_inactive_alpha::Real=0.5,
                            opoint_color=JULIA_GREEN,
                            opoint_marker=:star5)
    r, z = _rzplot(gs)
    faces = _faces1(gs)
    eq = equilibrium === nothing ? gs.equilibrium : equilibrium
    eq === nothing && error("No equilibrium available for plot_psi")

    psi_vals = if psi === nothing
        Vector{Float64}(EquilibriumModule.get_psi(eq; normalized=normalized))
    else
        Vector{Float64}(psi)
    end
    if normalized && gs.psi_convention == 0
        psi_vals = 1.0 .- psi_vals
    end

    plvls = if plasma_levels === nothing
        normalized ?
            collect(range(0.0, 1.0; length=plasma_nlevels)) :
            collect(range(minimum(psi_vals), maximum(psi_vals); length=plasma_nlevels))
    else
        l = collect(Float64.(plasma_levels))
        if normalized && gs.psi_convention == 0
            sort(1.0 .- l)
        else
            sort(l)
        end
    end

    vlvls = if vacuum_levels === nothing
        v = Float64[]
        if normalized
            if minimum(psi_vals) < -0.1
                append!(v, collect(range(minimum(psi_vals), 0.0;
                                          length=vacuum_nlevels + 1))[1:end - 1])
            end
            if maximum(psi_vals) > 1.1
                append!(v, collect(range(1.0, maximum(psi_vals);
                                          length=vacuum_nlevels + 1))[1:end - 1])
            end
        end
        v
    else
        l = collect(Float64.(vacuum_levels))
        if normalized && gs.psi_convention == 0
            sort(1.0 .- l)
        else
            sort(l)
        end
    end

    plasma_mask = Vector(gs.reg .== Int32(1))
    vacuum_mask = .!plasma_mask
    _draw_levels!(ax, r, z, faces, psi_vals, vlvls,
                  vacuum_color, vacuum_colormap, vacuum_linestyle, vacuum_linewidth;
                  cell_mask=vacuum_mask)
    _draw_levels!(ax, r, z, faces, psi_vals, plvls,
                  plasma_color, plasma_colormap, plasma_linestyle, plasma_linewidth;
                  cell_mask=plasma_mask)

    # X-points
    if xpoint_color !== nothing
        xp, diverted = EquilibriumModule.get_xpoints(eq)
        if size(xp, 1) > 0
            if diverted
                Makie.scatter!(ax, [xp[end, 1]], [xp[end, 2]];
                               color=xpoint_color, marker=xpoint_marker, markersize=15)
                if size(xp, 1) > 1
                    Makie.scatter!(ax, xp[1:end - 1, 1], xp[1:end - 1, 2];
                                   color=(xpoint_color, xpoint_inactive_alpha),
                                   marker=xpoint_marker, markersize=15)
                end
            else
                Makie.scatter!(ax, xp[:, 1], xp[:, 2];
                               color=(xpoint_color, xpoint_inactive_alpha),
                               marker=xpoint_marker, markersize=15)
            end
        end
    end

    # O-point
    if opoint_color !== nothing
        op = EquilibriumModule.o_point(eq)
        if op[1] > 0.0
            Makie.scatter!(ax, [op[1]], [op[2]];
                           color=opoint_color, marker=opoint_marker, markersize=15)
        end
    end

    ax.aspect = Makie.DataAspect()
    return nothing
end

TokaMaker.plot_psi(ax::Makie.AbstractAxis, gs::Tokamaker; kwargs...) =
    TokaMaker.plot_psi(nothing, ax, gs; kwargs...)

function TokaMaker.plot_psi(gs::Tokamaker; kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; xlabel="R (m)", ylabel="Z (m)")
    TokaMaker.plot_psi(fig, ax, gs; kwargs...)
    return fig, ax
end

# ---------------------------------------------------------------------------
# plot_constraints
#
# The Julia port does not currently store constraints on the equilibrium
# object, so callers pass them in directly. This mirrors what
# `set_isoflux_constraints!` / `set_saddle_constraints!` accept.

function TokaMaker.plot_constraints(fig, ax::Makie.AbstractAxis, gs::Tokamaker;
                                    isoflux::Union{Nothing,AbstractMatrix}=nothing,
                                    saddles::Union{Nothing,AbstractMatrix}=nothing,
                                    isoflux_color=JULIA_BLUE,
                                    isoflux_marker=:cross,
                                    saddle_color=JULIA_GREEN,
                                    saddle_marker=:xcross)
    if isoflux !== nothing && isoflux_color !== nothing && size(isoflux, 1) > 0
        Makie.scatter!(ax, Float64.(isoflux[:, 1]), Float64.(isoflux[:, 2]);
                       color=isoflux_color, marker=isoflux_marker, markersize=12)
    end
    if saddles !== nothing && saddle_color !== nothing && size(saddles, 1) > 0
        Makie.scatter!(ax, Float64.(saddles[:, 1]), Float64.(saddles[:, 2]);
                       color=saddle_color, marker=saddle_marker, markersize=12)
    end
    return nothing
end

TokaMaker.plot_constraints(ax::Makie.AbstractAxis, gs::Tokamaker; kwargs...) =
    TokaMaker.plot_constraints(nothing, ax, gs; kwargs...)

end # module

module Meshing

# Native Julia port of TokaMaker's `meshing.py`. Three layers:
#
# 1. `GsDomain` — top-level user API mirroring Python's `gs_Domain`.
# 2. `Region` — single closed contour with corner detection and resampling.
# 3. `MeshBuilder` — assembles unique points, dedups segments across regions,
#    then calls Triangulate.jl to produce the final triangulation.
#
# HDF5 schema for `save_gs_mesh`/`load_gs_mesh` matches Python exactly so the
# files are interchangeable between the two language fronts.

using HDF5
using JSON3
using LinearAlgebra: norm
using OrderedCollections: OrderedDict
using Triangulate: TriangulateIO, triangulate

export save_gs_mesh, load_gs_mesh, GsDomain, define_region!,
       add_polygon!, add_rectangle!, add_annulus!, add_enclosed!,
       build_mesh!, get_coils, get_conductors

# ----------------------------------------------------------------------------
# HDF5 / JSON I/O (unchanged from earlier)

function save_gs_mesh(pts::AbstractMatrix, tris::AbstractMatrix,
                      regions::AbstractVector, coil_dict::AbstractDict,
                      cond_dict::AbstractDict, filename::AbstractString;
                      use_hdf5::Bool=true)
    if use_hdf5
        coil_json = JSON3.write(coil_dict)
        cond_json = JSON3.write(cond_dict)
        h5open(filename, "w") do f
            f["mesh/r"] = Matrix{Float64}(pts)
            f["mesh/lc"] = Matrix{Int32}(tris)
            f["mesh/reg"] = Vector{Int32}(regions)
            f["mesh/coil_dict"] = String(coil_json)
            f["mesh/cond_dict"] = String(cond_json)
        end
    else
        d = Dict(
            "mesh" => Dict(
                "r" => collect(eachrow(pts)),
                "lc" => collect(eachrow(tris)),
                "reg" => collect(regions),
                "coil_dict" => coil_dict,
                "cond_dict" => cond_dict,
            )
        )
        open(filename, "w") do io
            JSON3.write(io, d)
        end
    end
    return filename
end

function load_gs_mesh(filename::AbstractString; use_hdf5::Bool=true)
    if use_hdf5
        h5open(filename, "r") do f
            pts_raw = read(f["mesh/r"])
            tris_raw = read(f["mesh/lc"])
            pts = ndims(pts_raw) == 2 && size(pts_raw, 1) == 2 ?
                  Matrix{Float64}(pts_raw') : Matrix{Float64}(pts_raw)
            tris = ndims(tris_raw) == 2 && size(tris_raw, 1) == 3 ?
                   Matrix{Int32}(tris_raw') : Matrix{Int32}(tris_raw)
            regions = Vector{Int32}(read(f["mesh/reg"]))
            coil_str = read(f["mesh/coil_dict"])
            cond_str = read(f["mesh/cond_dict"])
            coil_dict = _json_to_orderdict(coil_str)
            cond_dict = _json_to_orderdict(cond_str)
            return pts, tris, regions, coil_dict, cond_dict
        end
    else
        d = open(filename, "r") do io
            JSON3.read(read(io, String), Dict{String,Any})
        end
        m = d["mesh"]
        pts = Matrix{Float64}(reduce(hcat, m["r"])')
        tris = Matrix{Int32}(reduce(hcat, m["lc"])')
        regions = Vector{Int32}(m["reg"])
        return pts, tris, regions, m["coil_dict"], m["cond_dict"]
    end
end

function _json_to_orderdict(s::AbstractString)
    raw = JSON3.read(s)
    return _to_orderdict(raw)
end
_to_orderdict(x::JSON3.Object) =
    OrderedDict{String,Any}(String(k) => _to_orderdict(v) for (k, v) in pairs(x))
_to_orderdict(x::JSON3.Array) = [_to_orderdict(v) for v in x]
_to_orderdict(x) = x

# ----------------------------------------------------------------------------
# Region — single closed contour

mutable struct Region
    points::Matrix{Float64}              # [nv, 2] vertex list
    segments::Vector{Vector{Int}}        # 1-based local indices into `points`
    id::Int
    dx_vol::Float64
    dx_curve::Float64
    small_thresh::Float64
    resampled_points::Matrix{Float64}    # filled by MeshBuilder
end

function Region(points::AbstractMatrix; dx::Real, dx_curve::Union{Nothing,Real}=nothing,
                angle_tol::Real=30.0, sliver_tol::Real=120.0,
                small_thresh::Union{Nothing,Real}=nothing, id::Integer=0)
    pts = Matrix{Float64}(points)
    any(@view(pts[:, 1]) .< 0.0) &&
        error("Region contour has a negative radial coordinate")
    dx > 0 || error("dx must be positive")
    dxc = dx_curve === nothing ? Float64(dx) : Float64(dx_curve)
    # Drop duplicated last point (closed contour)
    if norm(pts[1, :] .- pts[end, :]) < dx / 1e3
        pts = pts[1:end-1, :]
    end
    nv = size(pts, 1)
    # Corner detection — splits the contour into segments at any vertex
    # where the tangent change exceeds `angle_tol`. Mirrors meshing.py:1119-1142.
    keep_tol = cos(π * angle_tol / 180)
    sliver = cos(π * sliver_tol / 180)
    keep = Int[1]
    tangp = pts[2, :] .- pts[1, :]
    tangp ./= norm(tangp)
    for i in 1:nv
        tang = if i == nv
            pts[1, :] .- pts[i, :]
        else
            pts[i+1, :] .- pts[i, :]
        end
        tn = norm(tang)
        if tn < dxc / 1e3
            @warn "Repeated points detected at vertex $i ($(pts[i, 1]), $(pts[i, 2]))"
            continue
        end
        tang ./= tn
        if i > 1
            ang = sum(tang .* tangp)
            if ang < keep_tol
                push!(keep, i)
                ang < sliver && @warn "Sliver (angle=$(round(180 - acos(ang)*180/π; digits=1))°) at vertex $i"
            end
        end
        tangp = tang
    end
    push!(keep, nv + 2)  # sentinel (matches Python's `nv+1` in 0-based)
    # Build segment list. Python uses 0-based indices and tests `i >= keep[k]`;
    # Julia is 1-based throughout, so the same comparison works after we
    # initialized `keep` and the loop variable both as 1-based.
    segments = Vector{Vector{Int}}()
    seg_tmp = Int[]
    k = 2
    for i in 1:nv
        if i >= keep[k]
            push!(segments, vcat(seg_tmp, [i]))
            seg_tmp = [i]
            k += 1
        else
            push!(seg_tmp, i)
        end
    end
    push!(segments, vcat(seg_tmp, [1]))   # close contour
    st = small_thresh === nothing ? dxc / 2 : Float64(small_thresh)
    return Region(pts, segments, Int(id), Float64(dx), dxc, st, zeros(Float64, 0, 2))
end

# Polygon point-in-region test (cast horizontal ray, count crossings).
# Mirrors `Region.check_in_poly`.
function check_in_poly(rp::AbstractMatrix, pt::AbstractVector)
    n = size(rp, 1)
    ncuts = 0
    for j in 1:n-1
        y0 = rp[j, 2]; y1 = rp[j+1, 2]
        if (pt[2] - y0) * (pt[2] - y1) <= 0.0
            (pt[2] - y0) == 0.0 && continue
            xint = (rp[j+1, 1] - rp[j, 1]) * (pt[2] - y0) / (y1 - y0) + rp[j, 1]
            pt[1] <= xint && (ncuts += 1)
        end
    end
    y0 = rp[end, 2]; y1 = rp[1, 2]
    if (pt[2] - y0) * (pt[2] - y1) <= 0.0
        (pt[2] - y0) == 0.0 && return isodd(ncuts)
        xint = (rp[1, 1] - rp[end, 1]) * (pt[2] - y0) / (y1 - y0) + rp[end, 1]
        pt[1] <= xint && (ncuts += 1)
    end
    return isodd(ncuts)
end

"Find an interior reference point near boundary point `i`. Mirrors
`Region.get_in_point` in meshing.py:1212."
function get_in_point(reg::Region, i::Int, dx::Float64)
    rp = reg.resampled_points
    n = size(rp, 1)
    dx_use = min(dx, reg.dx_curve) / 4.0
    that2 = if i == n
        rp[1, :] .- rp[i, :]
    else
        rp[i+1, :] .- rp[i, :]
    end
    that1 = rp[i, :] .- rp[i == 1 ? n : i-1, :]
    that1 ./= norm(that1); that2 ./= norm(that2)
    nhat = dx_use .* ([-that1[2], that1[1]] .+ [-that2[2], that2[1]])
    pt = rp[i, :] .+ nhat
    if !check_in_poly(rp, pt)
        pt = rp[i, :] .- nhat
    end
    return pt
end

# ----------------------------------------------------------------------------
# GsDomain

mutable struct GsDomain
    rextent::Union{Nothing,Float64}
    zextents::Vector{Union{Nothing,Float64}}
    rpad::Union{Nothing,Float64}
    zpad::Vector{Union{Nothing,Float64}}
    rmax::Float64
    zmin::Float64
    zmax::Float64
    boundary_reg::Union{Nothing,String}
    regions::Vector{Region}
    reg_type_counts::Dict{String,Int}
    region_info::OrderedDict{String,Dict{String,Any}}
    extra_reg_defs::Vector{Vector{Float64}}  # [r, z, region_id, max_area]
    r::Union{Nothing,Matrix{Float64}}
    lc::Union{Nothing,Matrix{Int32}}
    reg::Union{Nothing,Vector{Int32}}
end

function GsDomain(; rextent::Union{Nothing,Real}=nothing,
                    zextents::AbstractVector=Union{Nothing,Float64}[nothing, nothing],
                    rpad::Real=1.2, zpad::AbstractVector=[1.2, 1.2])
    rp = rextent === nothing ? Float64(rpad) : nothing
    zp_in = collect(Union{Nothing,Float64}, zextents)
    zp = Union{Nothing,Float64}[
        zp_in[1] === nothing ? Float64(zpad[1]) : nothing,
        zp_in[2] === nothing ? Float64(zpad[2]) : nothing,
    ]
    return GsDomain(
        rextent === nothing ? nothing : Float64(rextent),
        zp_in, rp, zp,
        0.0, 0.0, 0.0, nothing,
        Region[],
        Dict{String,Int}("plasma" => 0, "vacuum" => 0, "boundary" => 0,
                          "conductor" => 0, "coil" => 0),
        OrderedDict{String,Dict{String,Any}}(),
        Vector{Vector{Float64}}(),
        nothing, nothing, nothing,
    )
end

"""
    define_region!(dom, name, dx, reg_type; eta=nothing, noncontinuous=nothing,
                   nTurns=nothing, coil_set=nothing, allow_xpoints=false,
                   inner_limiter=false)

See Python `gs_Domain.define_region` for full semantics. Region names are
upper-cased on entry (matching Python).
"""
function define_region!(dom::GsDomain, name::AbstractString, dx::Real,
                        reg_type::AbstractString;
                        eta::Union{Nothing,Real}=nothing,
                        noncontinuous::Union{Nothing,Bool}=nothing,
                        nTurns::Union{Nothing,Real}=nothing,
                        coil_set::Union{Nothing,AbstractString}=nothing,
                        allow_xpoints::Bool=false,
                        inner_limiter::Bool=false)
    dx > 0 || error("dx must be positive")
    nm = uppercase(name)
    startswith(nm, "#") && error("Invalid region name (cannot start with '#')")
    haskey(dom.region_info, nm) && error("Region $nm already exists")
    next_id = -1
    if reg_type == "plasma"
        next_id = 1
        allow_xpoints = true
    elseif reg_type == "boundary"
        dom.boundary_reg = nm
    elseif reg_type ∉ ("vacuum", "conductor", "coil")
        error("Unknown region type \"$reg_type\"")
    end
    if next_id < 0
        next_id = length(dom.region_info) - dom.reg_type_counts["plasma"] + 2
    end
    dom.reg_type_counts[reg_type] += 1
    info = Dict{String,Any}(
        "id" => next_id, "dx" => Float64(dx), "count" => 0,
        "type" => reg_type, "allow_xpoints" => allow_xpoints,
    )
    if inner_limiter
        reg_type == "plasma" && error("Plasma cannot be inner limiter")
        info["inner_limiter"] = true
    end
    if eta !== nothing
        reg_type == "conductor" || error("eta only valid for conductor regions")
        info["eta"] = Float64(eta)
    elseif reg_type == "conductor"
        error("Resistivity not specified for conductor region $nm")
    end
    if noncontinuous !== nothing
        reg_type == "conductor" || error("noncontinuous only valid for conductor regions")
        info["noncontinuous"] = noncontinuous
    end
    if nTurns !== nothing
        reg_type == "coil" || error("nTurns only valid for coil regions")
        info["nturns"] = Float64(nTurns)
    end
    if coil_set !== nothing
        reg_type == "coil" || error("coil_set only valid for coil regions")
        startswith(coil_set, "#") && error("Invalid coil_set name (cannot start with '#')")
        info["coil_set"] = String(coil_set)
    end
    dom.region_info[nm] = info
    return dom
end

function _check_dx(info::Dict, name::String)
    haskey(info, "dx") || error("Resolution for region $name not defined")
    return Float64(info["dx"])
end

function _bump_extents!(dom::GsDomain, contour::AbstractMatrix)
    dom.rmax = max(dom.rmax, maximum(@view contour[:, 1]))
    dom.zmax = max(dom.zmax, maximum(@view contour[:, 2]))
    dom.zmin = min(dom.zmin, minimum(@view contour[:, 2]))
end

function add_annulus!(dom::GsDomain, inner_contour::AbstractMatrix,
                      inner_name::AbstractString, outer_contour::AbstractMatrix,
                      annulus_name::AbstractString;
                      parent_name::Union{Nothing,AbstractString}=nothing,
                      angle_tol::Real=30.0, sliver_tol::Real=120.0,
                      small_thresh::Union{Nothing,Real}=nothing)
    inner_name_u = uppercase(inner_name); annulus_name_u = uppercase(annulus_name)
    inner_info = get(dom.region_info, inner_name_u, nothing)
    annulus_info = get(dom.region_info, annulus_name_u, nothing)
    inner_info === nothing && error("Region $inner_name_u not defined")
    annulus_info === nothing && error("Region $annulus_name_u not defined")
    inner_dx = _check_dx(inner_info, inner_name_u)
    annulus_dx = _check_dx(annulus_info, annulus_name_u)
    inner_dx_curve = min(inner_dx, annulus_dx)
    outer_dx_curve = annulus_dx
    if parent_name !== nothing
        parent_u = uppercase(parent_name)
        parent_info = get(dom.region_info, parent_u, nothing)
        parent_info === nothing && error("Region $parent_u not defined")
        outer_dx_curve = min(outer_dx_curve, _check_dx(parent_info, parent_u))
    end
    inner_pts = Matrix{Float64}(inner_contour)
    outer_pts = Matrix{Float64}(outer_contour)
    _bump_extents!(dom, inner_pts); _bump_extents!(dom, outer_pts)
    push!(dom.regions, Region(inner_pts; dx=inner_dx, dx_curve=inner_dx_curve,
                              angle_tol=angle_tol, sliver_tol=sliver_tol,
                              small_thresh=small_thresh, id=inner_info["id"]))
    inner_info["count"] += 1
    push!(dom.regions, Region(outer_pts; dx=annulus_dx, dx_curve=outer_dx_curve,
                              angle_tol=angle_tol, sliver_tol=sliver_tol,
                              small_thresh=small_thresh, id=annulus_info["id"]))
    annulus_info["count"] += 1
    return dom
end

function add_polygon!(dom::GsDomain, contour::AbstractMatrix,
                     name::AbstractString;
                     parent_name::Union{Nothing,AbstractString}=nothing,
                     angle_tol::Real=30.0, sliver_tol::Real=120.0,
                     small_thresh::Union{Nothing,Real}=nothing)
    nm = uppercase(name)
    info = get(dom.region_info, nm, nothing)
    info === nothing && error("Region $nm not defined")
    dx = _check_dx(info, nm)
    dx_curve = dx
    if parent_name !== nothing
        parent_u = uppercase(parent_name)
        parent_info = get(dom.region_info, parent_u, nothing)
        parent_info === nothing && error("Region $parent_u not defined")
        dx_curve = min(dx, _check_dx(parent_info, parent_u))
    end
    pts = Matrix{Float64}(contour)
    _bump_extents!(dom, pts)
    push!(dom.regions, Region(pts; dx=dx, dx_curve=dx_curve,
                              angle_tol=angle_tol, sliver_tol=sliver_tol,
                              small_thresh=small_thresh, id=info["id"]))
    info["count"] += 1
    return dom
end

function add_rectangle!(dom::GsDomain, rc::Real, zc::Real, w::Real, h::Real,
                       name::AbstractString;
                       parent_name::Union{Nothing,AbstractString}=nothing,
                       rot::Union{Nothing,Real}=nothing)
    contour = Float64[
        -w/2 -h/2
        +w/2 -h/2
        +w/2 +h/2
        -w/2 +h/2
    ]
    if rot !== nothing
        θ = deg2rad(Float64(rot))
        R = Float64[cos(θ) -sin(θ); sin(θ) cos(θ)]
        contour = contour * R'
    end
    contour[:, 1] .+= rc
    contour[:, 2] .+= zc
    return add_polygon!(dom, contour, name; parent_name=parent_name)
end

function add_enclosed!(dom::GsDomain, in_point::AbstractVector,
                      name::AbstractString)
    nm = uppercase(name)
    info = get(dom.region_info, nm, nothing)
    info === nothing && error("Region $nm not defined")
    dx = _check_dx(info, nm)
    push!(dom.extra_reg_defs,
          [Float64(in_point[1]), Float64(in_point[2]),
           Float64(info["id"]), dx * dx / 2.0])
    info["count"] += 1
    return dom
end

function get_coils(dom::GsDomain)
    out = OrderedDict{String,Dict{String,Any}}()
    cid = 0
    for (key, info) in dom.region_info
        if info["type"] == "coil"
            out[key] = Dict{String,Any}(
                "reg_id" => info["id"],
                "coil_id" => cid,
                "nturns" => get(info, "nturns", 1.0),
                "coil_set" => get(info, "coil_set", key),
                "allow_xpoints" => get(info, "allow_xpoints", false),
            )
            cid += 1
        end
    end
    return out
end

function get_conductors(dom::GsDomain)
    out = OrderedDict{String,Dict{String,Any}}()
    cid = 0; vid = 0
    for (key, info) in dom.region_info
        t = info["type"]
        if t == "conductor"
            out[key] = Dict{String,Any}(
                "reg_id" => info["id"],
                "cond_id" => cid,
                "eta" => info["eta"],
                "noncontinuous" => get(info, "noncontinuous", false),
                "allow_xpoints" => get(info, "allow_xpoints", false),
            )
            haskey(info, "inner_limiter") &&
                (out[key]["inner_limiter"] = info["inner_limiter"])
            cid += 1
        elseif t in ("vacuum", "boundary")
            entry = Dict{String,Any}(
                "reg_id" => info["id"],
                "vac_id" => vid,
                "allow_xpoints" => get(info, "allow_xpoints", false),
            )
            haskey(info, "inner_limiter") &&
                (entry["inner_limiter"] = info["inner_limiter"])
            out[key] = entry
            vid += 1
        end
    end
    return out
end

# ----------------------------------------------------------------------------
# MeshBuilder — segment dedup + Triangulate.jl invocation
#
# Critical implementation note: the Python merging pass uses 0-based indices
# throughout (`reg_pt_map`, `tmp_pt_map`, `segment[0]`, `_resampled_segments`).
# Translating to Julia, we keep the exact same scheme — 0-based — and only
# convert at the boundary with Triangulate.jl (which expects 1-based vertex
# indices in segmentlist).

mutable struct MeshBuilder
    merge_thresh::Float64
    unique_points::Vector{Vector{Float64}}
    segments::Vector{Tuple{Vector{Int},Float64,Float64}}   # (point-id list, dx, small_thresh)
    reg_seg_map::Vector{Vector{Int}}                       # 0-based (matches Python sign convention)
    resampled_segments::Vector{Vector{Int}}
    resampled_points::Matrix{Float64}
    reg_defs::Vector{Vector{Float64}}
end

function MeshBuilder(regions::Vector{Region};
                     merge_thresh::Float64=1e-4,
                     debug::Bool=false,
                     extra_reg_defs::Vector{Vector{Float64}}=Vector{Vector{Float64}}())
    mb = MeshBuilder(merge_thresh,
                     Vector{Vector{Float64}}(),
                     Vector{Tuple{Vector{Int},Float64,Float64}}(),
                     Vector{Vector{Int}}(),
                     Vector{Vector{Int}}(),
                     zeros(Float64, 0, 2),
                     copy(extra_reg_defs))

    println("Assembling regions:")

    # Pass 1: build unique-point list and dedup segments across regions.
    # Faithful to meshing.py:861-994.
    for (ireg, region) in enumerate(regions)
        local_seg_map = Int[]
        # 0-based ids: each region's local point i maps to a global id.
        reg_pt_map = collect((length(mb.unique_points)):(length(mb.unique_points) + size(region.points, 1) - 1))
        ilocal = length(mb.unique_points) - 1
        seg_start = length(mb.segments)
        for (reg_seg_idx, tmp_pts_local) in enumerate(region.segments)
            tmp_pt_map = [reg_pt_map[i] for i in tmp_pts_local]
            for (ipt_local, reg_id) in enumerate(tmp_pts_local)
                reg_pt = region.points[reg_id, :]
                if reg_pt_map[reg_id] > length(mb.unique_points) - 1
                    matched = false
                    for (jpt, unique_pt) in enumerate(mb.unique_points)
                        if norm(reg_pt .- unique_pt) < merge_thresh
                            debug && println("  Merging point: ireg=$ireg map=$(reg_pt_map[reg_id]) -> $(jpt - 1) ($(reg_pt))")
                            tmp_pt_map[ipt_local] = jpt - 1
                            reg_pt_map[reg_id] = jpt - 1
                            matched = true
                            break
                        end
                    end
                    if !matched
                        ilocal += 1
                        tmp_pt_map[ipt_local] = ilocal
                        reg_pt_map[reg_id] = ilocal
                        push!(mb.unique_points, reg_pt)
                    end
                end
            end
            # Compare against existing segments to detect overlaps.
            matched_existing = false
            for (iseg, segment) in enumerate(mb.segments)
                seg_pts = segment[1]
                # Count adjacent matches between tmp_pt_map and seg_pts
                nOverlap = 0
                for i in 1:length(seg_pts)
                    test_pt = seg_pts[i == 1 ? end : i - 1]
                    matches = findall(==(test_pt), tmp_pt_map)
                    if !isempty(matches)
                        test_pt2 = seg_pts[i]
                        for m in findall(==(test_pt2), tmp_pt_map)
                            if (m - 1 in matches) || (m + 1 in matches)
                                nOverlap += 1
                            end
                        end
                    end
                end
                if nOverlap == length(tmp_pts_local) - 1 && length(tmp_pts_local) == length(seg_pts)
                    # Full overlap. Try forward then backward orientation.
                    fwd_match = all(tmp_pt_map[i] == seg_pts[i] for i in 1:length(seg_pts))
                    if fwd_match
                        debug && println("  Merging curve segments (fwd): ireg=$ireg iseg=$(iseg - 1)")
                        new_dx = min(segment[2], region.dx_curve)
                        new_st = min(segment[3], region.small_thresh)
                        mb.segments[iseg] = (seg_pts, new_dx, new_st)
                        push!(local_seg_map, iseg - 1)
                        matched_existing = true
                        break
                    end
                    bwd_match = all(tmp_pt_map[end - i + 1] == seg_pts[i] for i in 1:length(seg_pts))
                    if bwd_match
                        debug && println("  Merging curve segments (bwd): ireg=$ireg iseg=$(iseg - 1)")
                        new_dx = min(segment[2], region.dx_curve)
                        new_st = min(segment[3], region.small_thresh)
                        mb.segments[iseg] = (seg_pts, new_dx, new_st)
                        push!(local_seg_map, -(iseg - 1))
                        matched_existing = true
                        break
                    end
                elseif nOverlap >= 1 && (iseg - 1) < seg_start
                    # Partial match (only across regions). The Python branch
                    # for partial matches splits segments to share the overlap.
                    # In practice, the test geometries (Solov'ev, ITER) hit
                    # only the full-overlap branch above; partial matches
                    # arise when two non-trivially-shaped regions share an
                    # edge of differing length. We log this case so it's
                    # visible if it does come up; users should split their
                    # contours so segments align exactly until the partial
                    # branch is implemented.
                    @warn "Partial segment overlap not handled (ireg=$ireg, iseg=$(iseg - 1), nOverlap=$nOverlap)"
                end
            end
            if !matched_existing
                push!(mb.segments, (tmp_pt_map, region.dx_curve, region.small_thresh))
                push!(local_seg_map, length(mb.segments) - 1)
            end
        end
        push!(mb.reg_seg_map, local_seg_map)
    end

    # Pass 2: resample each segment to target dx, building a flat point list.
    pts_out = copy(mb.unique_points)
    upoints = Matrix{Float64}(undef, length(mb.unique_points), 2)
    for (i, p) in enumerate(mb.unique_points)
        upoints[i, 1] = p[1]; upoints[i, 2] = p[2]
    end
    resampled = Vector{Vector{Int}}()
    for segment in mb.segments
        seg_pts = segment[1]   # 0-based ids into upoints
        dx = segment[2]
        st = segment[3]
        # Map to coordinates
        coords = upoints[seg_pts .+ 1, :]
        dl = zeros(Float64, length(seg_pts))
        for i in 1:length(seg_pts) - 1
            dl[i + 1] = dl[i] + norm(coords[i + 1, :] .- coords[i, :])
        end
        seg_tmp = Int[seg_pts[1]]
        if floor(Int, dl[end] / dx) >= 2
            n_inner = floor(Int, dl[end] / dx)
            for dl_samp in collect(range(0.0, dl[end]; length=n_inner + 1))[2:end-1]
                # Linear interp along the polyline
                idx = clamp(searchsortedlast(dl, dl_samp), 1, length(dl) - 1)
                t = (dl_samp - dl[idx]) / (dl[idx + 1] - dl[idx])
                interp_pt = coords[idx, :] .* (1 - t) .+ coords[idx + 1, :] .* t
                push!(pts_out, interp_pt)
                push!(seg_tmp, length(pts_out) - 1)   # 0-based id
            end
        elseif dl[end] < st
            @warn "Small feature (dl=$(round(dl[end]; sigdigits=3))) at $(coords[1, :])"
        end
        push!(seg_tmp, seg_pts[end])
        push!(resampled, seg_tmp)
    end
    # Reindex points so only those used by segments end up in the final list.
    npts_total = length(pts_out)
    reindex = fill(-1, npts_total)
    pt_count = 0
    for segment in resampled
        for i in 1:length(segment) - 1
            if reindex[segment[i] + 1] < 0
                pt_count += 1
                reindex[segment[i] + 1] = pt_count - 1
            end
            if reindex[segment[i + 1] + 1] < 0
                pt_count += 1
                reindex[segment[i + 1] + 1] = pt_count - 1
            end
            segment[i] = reindex[segment[i] + 1]
        end
        segment[end] = reindex[segment[end] + 1]
    end
    resampled_pts = zeros(Float64, pt_count, 2)
    for i in 1:npts_total
        if reindex[i] >= 0
            resampled_pts[reindex[i] + 1, 1] = pts_out[i][1]
            resampled_pts[reindex[i] + 1, 2] = pts_out[i][2]
        end
    end
    mb.resampled_segments = resampled
    mb.resampled_points = resampled_pts

    println("  # of unique points    = $(size(resampled_pts, 1))")
    println("  # of unique segments  = $(length(mb.segments))")

    # Pass 3: build region "in-point" list for Triangle's region attribute pass.
    # Each region needs an interior reference point + max-area constraint.
    for (ireg, region) in enumerate(regions)
        imin = -1
        dmin = 1e99
        dmax = -1.0
        reg_points = Vector{Vector{Float64}}()
        for seg_id in mb.reg_seg_map[ireg]
            tmp_pt_map = if seg_id < 0
                reverse(resampled[-seg_id + 1])
            else
                resampled[seg_id + 1]
            end
            for (ipt, unique_id) in enumerate(tmp_pt_map)
                ipt == length(tmp_pt_map) && continue
                reg_pt = resampled_pts[unique_id + 1, :]
                push!(reg_points, reg_pt)
                dmin_loc = 1e99
                for jpt in 1:size(resampled_pts, 1)
                    (jpt - 1) == unique_id && continue
                    d = norm(reg_pt .- resampled_pts[jpt, :])
                    dmin_loc = min(dmin_loc, d)
                end
                dmin = min(dmin, dmin_loc)
                if dmin_loc > dmax
                    dmax = dmin_loc
                    imin = length(reg_points)
                end
            end
        end
        rp_mat = Matrix{Float64}(undef, length(reg_points), 2)
        for (i, p) in enumerate(reg_points)
            rp_mat[i, :] = p
        end
        region.resampled_points = rp_mat
        in_pt = get_in_point(region, imin, dmin)
        push!(mb.reg_defs,
              [in_pt[1], in_pt[2], Float64(region.id), region.dx_vol^2 / 2.0])
    end

    return mb
end

"""
    get_mesh!(mb) -> (pts::Matrix, tris::Matrix{Int32}, regions::Vector{Int32})

Run Triangle on the assembled segment soup and return the final mesh.
Triangle indices are returned 0-based (matching Python `gs_Domain.build_mesh`).
"""
function get_mesh!(mb::MeshBuilder)
    println("Generating mesh with Triangle:")
    # Flatten resampled_segments into pairwise edges.
    seg_pairs = Vector{Tuple{Int,Int}}()
    for segment in mb.resampled_segments
        for i in 1:length(segment) - 1
            push!(seg_pairs, (segment[i], segment[i + 1]))
        end
    end

    tio_in = TriangulateIO()
    # Triangulate.jl expects pointlist as (2, n) Float64.
    tio_in.pointlist = Matrix{Float64}(mb.resampled_points')
    # segmentlist (2, m) Int32, 1-based.
    seg_mat = Matrix{Cint}(undef, 2, length(seg_pairs))
    for (j, (a, b)) in enumerate(seg_pairs)
        seg_mat[1, j] = a + 1
        seg_mat[2, j] = b + 1
    end
    tio_in.segmentlist = seg_mat
    # regionlist (4, k) Float64: r, z, attribute, max_area.
    reg_mat = Matrix{Float64}(undef, 4, length(mb.reg_defs))
    for (j, rd) in enumerate(mb.reg_defs)
        reg_mat[1, j] = rd[1]; reg_mat[2, j] = rd[2]
        reg_mat[3, j] = rd[3]; reg_mat[4, j] = rd[4]
    end
    tio_in.regionlist = reg_mat

    # Switches: p = read PSLG, q = quality (default 20°), a = use max area
    # constraints from regionlist, A = output region attributes, Q = quiet.
    # We use 1-based output (omit `z`) since that matches Triangulate.jl's
    # default; convert to 0-based at return for compatibility with Python's
    # `gs_Domain.build_mesh` (which subtracts 1 inside `setup_mesh`).
    tio_out, _ = triangulate("pqaAQ", tio_in)

    # tio_out.pointlist is (2, np); transpose to (np, 2). Triangle list (3, nc).
    pts_mat = Matrix{Float64}(tio_out.pointlist')
    tris_raw = tio_out.trianglelist
    # Convert 1-based → 0-based to match Python convention.
    tris_mat = Matrix{Int32}(tris_raw' .- Int32(1))
    attrs = tio_out.triangleattributelist  # (1, nc)
    regs = Vector{Int32}(undef, size(tris_mat, 1))
    for i in 1:length(regs)
        regs[i] = Int32(round(attrs[1, i]))
    end

    println("  # of points  = $(size(pts_mat, 1))")
    println("  # of cells   = $(size(tris_mat, 1))")
    println("  # of regions = $(maximum(regs))")
    return pts_mat, tris_mat, regs
end

"""
    build_mesh!(dom; debug=false, merge_thresh=1e-4, require_boundary=true)

Generate the final triangulation. Returns `(pts[np, 2], lc[nc, 3], reg[nc])`.
`lc` is 0-based to match Python's `gs_Domain.build_mesh`; the `setup_mesh!`
helper auto-detects this.
"""
function build_mesh!(dom::GsDomain; debug::Bool=false,
                     merge_thresh::Real=1e-4,
                     require_boundary::Bool=true)
    dom.reg_type_counts["plasma"] == 0 && error("No plasma region defined")
    has_other = (dom.reg_type_counts["vacuum"] > 0 ||
                 dom.reg_type_counts["coil"] > 0 ||
                 dom.reg_type_counts["conductor"] > 0)
    if has_other && require_boundary
        dom.boundary_reg === nothing && error("No boundary region defined")
        if dom.region_info[dom.boundary_reg]["count"] == 0
            # Auto-create boundary contour from extents/padding.
            if dom.rextent === nothing
                dom.rextent = dom.rpad * dom.rmax
            elseif dom.rmax > dom.rextent
                error("rextent does not enclose all regions")
            end
            if dom.zextents[1] === nothing
                dom.zextents[1] = dom.zpad[1] * dom.zmin
            elseif dom.zmin < dom.zextents[1]
                error("zextents[1] does not enclose all regions")
            end
            if dom.zextents[2] === nothing
                dom.zextents[2] = dom.zpad[2] * dom.zmax
            elseif dom.zmax > dom.zextents[2]
                error("zextents[2] does not enclose all regions")
            end
            vac_dx = dom.region_info[dom.boundary_reg]["dx"]
            contour = Float64[
                0.0           dom.zextents[1]
                dom.rextent   dom.zextents[1]
                dom.rextent   dom.zextents[2]
                0.0           dom.zextents[2]
            ]
            push!(dom.regions,
                  Region(contour; dx=vac_dx, id=dom.region_info[dom.boundary_reg]["id"]))
            dom.region_info[dom.boundary_reg]["count"] += 1
        end
    end
    # Sanity: every region must have been used at least once.
    for (key, info) in dom.region_info
        info["count"] == 0 && error("Region $key defined but never created")
    end
    # Re-index regions in canonical order: plasma=1, then boundary, vacuum,
    # conductor, coil. Mirrors meshing.py:592-604.
    n_regs = length(dom.region_info)
    reg_reorder = fill(-1, n_regs)
    reg_reorder[1] = 1
    reg_id = 1
    for reg_type in ("boundary", "vacuum", "conductor", "coil")
        for (key, info) in dom.region_info
            if info["type"] == reg_type
                reg_id += 1
                reg_reorder[info["id"]] = reg_id
                info["id"] = reg_id
            end
        end
    end
    for region in dom.regions
        region.id = reg_reorder[region.id]
    end
    for pdef in dom.extra_reg_defs
        pdef[3] = Float64(reg_reorder[Int(pdef[3])])
    end

    mb = MeshBuilder(dom.regions; merge_thresh=Float64(merge_thresh),
                     debug=debug, extra_reg_defs=dom.extra_reg_defs)
    pts, lc, reg = get_mesh!(mb)
    minimum(reg) <= 0 && error("Meshing error: unclaimed region detected")
    dom.r = pts; dom.lc = lc; dom.reg = reg
    return pts, lc, reg
end

end # module

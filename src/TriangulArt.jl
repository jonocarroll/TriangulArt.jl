module TriangulArt

using Images
using FileIO
using StatsBase
using DelaunayTriangulation
using StableRNGs
using Plots

using PyCall

const skimage = PyNULL()

function __init__()
    copy!(skimage, pyimport("skimage"))
end

function edge_points(image::Array{RGB{N0f8},3},
    length_scale::Integer=200,
    n_horizontal_points::Integer=0,
    n_vertical_points::Integer=0,
)
    ymax, xmax = size(image)[1:2]
    if n_horizontal_points == 0
        n_horizontal_points = round(Int, xmax / length_scale)
    end
    if n_vertical_points == 0
        n_vertical_points = round(Int, ymax / length_scale)
    end

    Δx = xmax / n_horizontal_points
    Δy = ymax / n_vertical_points

    res = [0 0; 0 ymax; xmax 0; xmax ymax]
    tmp = zeros(n_horizontal_points, 2)
    tmp[:, 1] = Δx * (1:n_horizontal_points)
    res = vcat(res, tmp)
    tmp[:, 2] .= ymax
    res = vcat(res, tmp)
    tmp = zeros(n_vertical_points, 2)
    tmp[:, 2] = Δy * (1:n_vertical_points)
    res = vcat(res, tmp)
    tmp[:, 1] .= xmax
    res = vcat(res, tmp)
    unique(res, dims=1)
end

function generate_uniform_random_points(image::Matrix{RGB{N0f8}}, n_points::Integer=100)
    ymax, xmax = size(image)[1:2]
    rand(n_points, 2) .* [xmax ymax]
end

function generate_weighted_points(image::Matrix{RGB{N0f8}}, n_points::Integer=100)
    im1 = Gray.(image) * 255 |> x -> x .|> UInt8
    xsize, ysize = size(im1)
    length_scale = sqrt(xsize * ysize / n_points)
    sigma = 0.3 * length_scale
    radius = floor(Int, length_scale)

    amp = 3
    entropy_width = 0.2 * length_scale
    im2 = skimage.filters.rank.entropy(im1, skimage.morphology.disk(entropy_width))

    epts = Array{Int}(undef, (n_points, 2))
    for n in 1:n_points
        maxpos = last(findmax(im2))
        x, y = getindex(maxpos, 1), getindex(maxpos, 2)
        epts[n, :] = [y x]
        for i in (x-2radius):(x+2radius)
            for j in (y-2radius):(y+2radius)
                if !checkbounds(Bool, im2, i, j)
                    continue
                end
                im2[i, j] = im2[i, j] - amp * exp(-((i - x)^2 + (j - y)^2) / 2sigma^2)
            end
        end
    end
    unique(epts, dims=1)
end

function which_tri(tri::Triangulation, x, y, k=0)
    q = (float(x), float(y))
    V = iszero(k) ? jump_and_march(tri, q) : jump_and_march(tri, q; k)
    return DelaunayTriangulation.sort_triangle(V)
end

"""
    triangulArt(
        image::Matrix{RGB{N0f8}};
        npts::Integer=100,
        fast::Bool=false,
        refine::Bool=false,
        debug::Bool=false,
        showimage::Bool=false,
        showimagecol::Bool=false
    )

Create an artistic Delaunay Triangulation from an image

# Arguments
- `image::Matrix{RGB{N0f8}}`: an image loaded via `FileIO::load()`
- `npts::Integer=100`: the number of triangle vertices to create 
- `fast::Bool=false`: if `true`, use random sampling, otherwise use the local entropy to select the vertices
- `refine::Bool=false`: use `DelaunayTriangulation::refine()` to refine the vertices
- `debug::Bool=false`: use debug mode; plot the triangle vertices and edges
- `showimage::Bool=false`: plot a grayscale version of the original image beneath the points
- `showimagecol::Bool=false`: plot a colour version of the original image beneath the points. Setting this helps hide the triangle edges.

"""
function triangulArt(
    image::Matrix{RGB{N0f8}};
    npts::Integer=100,
    fast::Bool=false,
    refine::Bool=false,
    debug::Bool=false,
    showimage::Bool=false,
    showimagecol::Bool=true
)
    orig_image = copy(image)
    # use a 3D representation
    image = image[:, :, :]
    if fast
        points = generate_uniform_random_points(orig_image, npts)
    else
        points = generate_weighted_points(orig_image, npts)
    end
    boundary_points = edge_points(image)
    all_points = vcat(points, boundary_points)
    all_points = @views unique!(tuple.(all_points[:, 1], all_points[:, 2]))
    rng = StableRNG(2)
    tri = triangulate(all_points; rng)
    if refine
        refine!(tri, max_area=prod(size(image)) / 200, min_area=prod(size(image)) / (2000); rng)
    end

    # change from color, width, height to height, width, color
    img_CHW = channelview(image)
    img_HWC = permutedims(img_CHW, (3, 2, 1, 4)) # 2 * 2 * 3
    # correctly order height, width
    img_HWC = permutedims(img_HWC[:, :, :, 1], (2, 1, 3))

    locs = Vector{Tuple{Int,Int,Vector{N0f8},NTuple{3,Int}}}()
    sizey, sizex = size(img_HWC)
    k = 0
    for x in 1:10:sizex
        for y in 1:10:sizey
            V = which_tri(tri, x, y, k)
            k = rand(rng, triangle_vertices(V)) # just take some random vertex to start from 
            push!(locs, (x, y, img_HWC[y, x, :], V))
        end
    end
    # average colours
    at = getindex.(locs, 4)
    ind = map(elem -> findall(isequal(elem), at), unique(at))
    cols = Dict{NTuple{3,Int},NTuple{3,Float32}}()
    for i in ind
        tr = locs[i][1][4]
        r = map(e -> e[3][1], locs[i])
        g = map(e -> e[3][2], locs[i])
        b = map(e -> e[3][3], locs[i])
        cols[tr] = (median(r), median(g), median(b))
    end

    if showimage
        @views pl = plot(Gray.(img_HWC[:, :, 1, 1]), showaxis=false, xlim=(1, sizex), ylim=(1, sizey), aspect_ratio=:auto)
    elseif showimagecol
        pl = plot(orig_image, showaxis=false, xlim=(1, sizex), ylim=(1, sizey), aspect_ratio=:auto)
    else
        pl = plot(showaxis=false, yflip=true, xlim=(1, sizex), ylim=(1, sizey))
    end
    @views if debug
        scatter!(pl, boundary_points[:, 1], boundary_points[:, 2], markercolor=:blue, legend=false)
        scatter!(pl, points[:, 1], points[:, 2], markercolor=:red, legend=false)
    end

    for V in each_solid_triangle(tri)
        V = DelaunayTriangulation.sort_triangle(V) # note that we sorted earlier. This ensures that (1) all triangulations are represented in the same order (e.g. (3, 1, 2) and (1, 2, 3) both go to (1, 2, 3)), and that (2) the keys in the Dict only store the triangle once
        u, v, w = triangle_vertices(V)
        p, q, r = get_point(tri, u, v, w)
        xs = [getx(p), getx(q), getx(r), getx(p)]
        ys = [gety(p), gety(q), gety(r), gety(p)]
        cr, cg, cb = get(cols, V, (0.0f0, 0.0f0, 0.0f0))
        tricol = RGBA(cr, cg, cb, 1.0)
        if debug
            plot!(pl, xs, ys, linecolor=:black, legend=false)
        else
            plot!(pl, xs, ys, linewidth=0, linealpha=0, seriestype=:shape, fillcolor=tricol, legend=false)
        end
    end

    pl

end

export triangulArt

end

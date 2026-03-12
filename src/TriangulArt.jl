module TriangulArt

using Images
using FileIO
using StatsBase
using DelaunayTriangulation
using StableRNGs
using Plots
import Colors: clamp01

using PyCall

const skimage = PyNULL()

function __init__()
    copy!(skimage, pyimport("skimage"))
end

function edge_points(image::Matrix{RGB{N0f8}},
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
    ysize, xsize = size(im1) 
    length_scale = sqrt(xsize * ysize / n_points)
    sigma = 0.3 * length_scale
    radius = floor(Int, length_scale)

    amp = 3
    entropy_width = 0.2 * length_scale
    im2 = Float64.(skimage.filters.rank.entropy(im1, skimage.morphology.disk(entropy_width)))

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
        debug::Bool=false
    )

Create an artistic Delaunay Triangulation from an image

# Arguments
- `image::Matrix{RGB{N0f8}}`: an image loaded via `FileIO::load()`
- `npts::Integer=100`: the number of triangle vertices to create 
- `fast::Bool=false`: if `true`, use random sampling, otherwise use the local entropy to select the vertices
- `refine::Bool=false`: use `DelaunayTriangulation::refine()` to refine the vertices
- `debug::Bool=false`: use debug mode; overlay the triangle vertices and edges on the output

"""
function triangulArt(
    image::Matrix{RGB{N0f8}};
    npts::Integer=100,
    fast::Bool=false,
    refine::Bool=false,
    debug::Bool=false
)

    @info "triangulArt called" npts=npts

    orig_image = copy(image)
    orig_image = map(Images.clamp01nan, orig_image)

    # use a 3D representation
    # image = image[:, :, :]
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

    # prepare image array used for colouring
    img_CHW = channelview(image)
    img_HWC = permutedims(img_CHW, (2, 3, 1))
    img_HWC = clamp01.(img_HWC)
    sizey, sizex = size(img_HWC)          # height, width

    # sample the image colour at the centroid of triangle (p, q, r)
    function _centroid_colour(p, q, r)
        cx = clamp(round(Int, (getx(p) + getx(q) + getx(r)) / 3), 1, sizex)
        cy = clamp(round(Int, (gety(p) + gety(q) + gety(r)) / 3), 1, sizey)
        RGB{Float32}(
            clamp01(Float32(img_HWC[cy, cx, 1])),
            clamp01(Float32(img_HWC[cy, cx, 2])),
            clamp01(Float32(img_HWC[cy, cx, 3]))
        )
    end

    # build triangle list and pre-compute colours
    tri_list = [DelaunayTriangulation.sort_triangle(V) for V in each_solid_triangle(tri)]
    cols_list = [begin
        u, v, w = triangle_vertices(V)
        p, q, r = get_point(tri, u, v, w)
        _centroid_colour(p, q, r)
    end for V in tri_list]

    # build triangle → colour lookup (keyed by sorted triangle)
    tri_colors = Dict(V => c for (V, c) in zip(tri_list, cols_list))

    # rasterize: for each output pixel find its containing triangle via
    # jump-and-march, then look up the pre-computed colour
    output = Matrix{RGB{Float32}}(undef, sizey, sizex)
    for row in 1:sizey, col in 1:sizex
        V = which_tri(tri, col, row)
        output[row, col] = get(tri_colors, V,
            RGB{Float32}(clamp01(Float32(img_HWC[row, col, 1])),
                         clamp01(Float32(img_HWC[row, col, 2])),
                         clamp01(Float32(img_HWC[row, col, 3]))))
    end

    pl = plot(output, showaxis=false, aspect_ratio=:auto)

    if debug
        scatter!(pl, boundary_points[:, 1], boundary_points[:, 2],
                 markercolor=:blue, legend=false)
        scatter!(pl, points[:, 1], points[:, 2], markercolor=:red,
                 legend=false)
        for V in tri_list
            u, v, w = triangle_vertices(V)
            p, q, r = get_point(tri, u, v, w)
            xs = [getx(p), getx(q), getx(r), getx(p)]
            ys = [gety(p), gety(q), gety(r), gety(p)]
            plot!(pl, xs, ys, linecolor=:black, legend=false)
        end
    end

    pl
end

export triangulArt

end

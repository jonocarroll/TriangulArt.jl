using TriangulArt
using Test

using FileIO
using ReferenceTests

refimg_path = joinpath(dirname(dirname(pathof(TriangulArt))), "test", "readme_images")
readmeimg_path = joinpath(dirname(dirname(pathof(TriangulArt))), "examples")

@testset "TriangulArt" begin
    bird = load(joinpath(readmeimg_path, "bird_small.jpg"))
    readme1 = triangulArt(bird, npts=50, fast=true) # don't test this: rng is not controlled 
    # Another problem with the above is that, due to the random sampling, 
    # not all of the image will be covered by a color.
    readme2 = triangulArt(bird)
    @test_reference joinpath(refimg_path, "readme2.png") readme2
    everest = load(joinpath(readmeimg_path, "everest.jpeg"))
    readme3 = triangulArt(everest)
    @test_reference joinpath(refimg_path, "readme3.png") readme3
    readme4 = triangulArt(bird, npts=500)
    @test_reference joinpath(refimg_path, "readme4.png") readme4
    readme5 = triangulArt(bird, npts=100, refine=true)
    @test_reference joinpath(refimg_path, "readme5.png") readme5
    readme6 = triangulArt(bird, npts=100, debug=true)
    @test_reference joinpath(refimg_path, "readme6.png") readme6
    readme7 = triangulArt(bird, npts=100, showimagecol=false, debug=true)
    @test_reference joinpath(refimg_path, "readme7.png") readme7
end
# TriangulArt

[![Build Status](https://github.com/jonocarroll/TriangulArt.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jonocarroll/TriangulArt.jl/actions/workflows/CI.yml?query=branch%3Amain)

```
using FileIO
using TriangulArt

bird = load("examples/bird_small.jpg");
triangulArt(bird)

triangulArt(bird, 500)
triangulArt(bird, 500, refine=true)
```
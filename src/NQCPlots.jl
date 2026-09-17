"""
A module containing Makie recipes for NQCBase structures.

I wish this could've been an extension of NQCBase, but the Julia extension system doesn't allow for it.
"""
module NQCPlots

using Makie: register_computation!
using Makie
using PeriodicTable
using UnitfulAtomic
using Unitful
using NQCBase

# Default atom colours and radii:
include("defaults.jl")

export atomic_structures_theme

const atomic_structures_theme = Theme(
    Axis=( # Fit axes to data by default, so atoms are circular.
        aspect=Makie.DataAspect(),
    ),
    Axis3=( # Fit axes to data by default, so atoms are spherical.
        aspect=:data,
    ),
)

bettersphere = Makie.GeometryBasics.mesh(
    Makie.GeometryBasics.Tesselation(
        Makie.GeometryBasics.Sphere(Point3f(0.0, 0.0, 0.0), 1.0),
        128
    )
)

include("atoms3d.jl")

include("atoms2d.jl")

end

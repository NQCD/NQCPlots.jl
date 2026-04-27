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

@recipe Atoms3D (structure,) begin
    """
    Dictionary{Symbol, Float64} mapping element symbols to atomic radii in Angstroms for plotting. Default values are based on empirical atomic radii in picometers from Wikipedia data page.

    e.g. `:C => 0.7` for carbon atoms with a radius of 0.7 Å.
    """
    atomicradii = Dict([el => ustrip(uconvert(u"Å", rad * u"pm")) for (el, rad) in default_atomic_radii_pm]...)
    """
    Dictionary{Symbol, Colorant} mapping element symbols to colors for plotting. Default values are based on CPK coloring scheme taken from PeriodicTable.jl.

        e.g. `:C => colorant"black"` for carbon atoms colored black.
    """
    atomcolors = default_atom_fills
    """
    Colour multiplier used to darken the atom edge colour relative to the fill colour. Default value is 0.5, which means the edge colour will be half as bright as the fill colour.
    """
    atomsstrokecolormultiplier = 0.5
    """
    The marker to use for plotting atoms. Default is a high-resolution sphere, which might make vector graphics larger.
    Consider using `rasterize=true` in the plot attributes to mitigate this.
    """
    marker = bettersphere
    "Atom edge width in pt, default is 1.5 pt."
    strokewidth = 1.5
    "Whether to show cell boundaries. Default is true."
    show_cell = true
    "Colour for the cell boundaries. Default is black."
    cellcolor = colorant"black"
    "Line width for the cell boundaries in pt. Default is 1.0 pt."
    celllinewidth = 1.0
    "Alpha transparency for the cell boundaries. Default is 0.5."
    cellalpha = 1.0
end

# Plot in 3D by default
# Makie.args_preferred_axis(::Type{<: Atoms3D}) = Makie.Axis3

function Makie.plot!(plot::Atoms3D{<:Tuple{<:NQCBase.Structure}},)
    input_nodes = [
        :structure,
        :atomicradii,
        :atomcolors,
        :show_cell,
    ] # Take a single structure in
    output_nodes = [
        :positions, # Convert structure positions to plottable positions in Å.
        :markersize, # Generate marker sizes based on conversion dict
        :color, # Generate colors based on conversion dict
    ]
    register_computation!(plot.attributes, input_nodes, output_nodes) do inputs, changed, cached
        # Unpack inputs
        atomic_structure = inputs[:structure]
        atomicradii = inputs[:atomicradii]
        atomcolors = inputs[:atomcolors]

        # Convert positions to Å and then to Point3f format for plotting
        positions_conv = auconvert.(u"Å", atomic_structure.positions) .|> ustrip
        positions_AA::Vector{Point3f} = [Point3f(pos...) for pos in eachcol(positions_conv)]
        # Generate marker sizes based on conversion dict
        marker_sizes = [atomicradii[el] for el in atomic_structure.atoms.types]
        # Generate colors based on conversion dict
        color = [atomcolors[el] for el in atomic_structure.atoms.types]
        return (positions_AA, marker_sizes, color)
    end
    register_computation!(
        plot.attributes,
        [:structure, :show_cell],
        [:cell_poly],
    ) do inputs, changed, cached
        show_cell = inputs[:show_cell]
        atomic_structure = inputs[:structure]
        if show_cell && isa(atomic_structure.cell, PeriodicCell)
            # Generate Cell polygon from structure.cell.vectors if show_cell is true
            cell_poly = Point3f[]
            cell_vectors = auconvert.(u"Å", atomic_structure.cell.vectors) |> ustrip # Convert to Å
            # Convert unit cell into vertices and indices for mesh plotting. Each column of cell vectors is the unit vector in one dimension.
            poly_points = [
                zero(cell_vectors[:, 1]), # Origin
                cell_vectors[:, 1],
                cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 2],
                zero(cell_vectors[:, 1]), # Bottom face done
                cell_vectors[:, 3],
                cell_vectors[:, 3] + cell_vectors[:, 1],
                cell_vectors[:, 3] + cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 3] + cell_vectors[:, 2],
                cell_vectors[:, 3], # Top face done
                cell_vectors[:, 3] + cell_vectors[:, 1],
                cell_vectors[:, 1],
                cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 3] + cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 3] + cell_vectors[:, 2],
                cell_vectors[:, 2],
            ]
            for vertex in poly_points
                push!(cell_poly, Point3f(vertex...))
            end
            return (cell_poly,)
        end
        return (Point3f[],)
    end

    # map!(plot.attributes, input_nodes, output_nodes) do (atomic_structure, atomicradii, atomcolors)
    #     # Convert positions to Å and then to Point3f format for plotting
    #     positions_conv = auconvert.(u"Å", atomic_structure.positions) .|> ustrip
    #     positions_AA::Vector{Point3f} = [Point3f(pos...) for pos in eachcol(positions_conv)]
    #     # Generate marker sizes based on conversion dict
    #     marker_sizes = [atomicradii[el] for el in atomic_structure.atoms.types]
    #     # Generate colors based on conversion dict
    #     color = [atomcolors[el] for el in atomic_structure.atoms.types]
    #     return (positions_AA, marker_sizes, color)
    # end

    # Generate positions in Angstroms and convert to Point3f format for plotting
    meshscatter!(
        plot,
        plot.positions,
        marker=plot[:marker],
        markersize=plot.markersize,
        color=plot.color,
        # markerspace = :data,
        #ssao = true,
        # strokewidth = atom_strokewidth,
        # strokecolor = [get(atom_edges, i, colorant"black") for i in trj_postprocessed[:atoms].types[atoms_to_draw]],
        rasterize=true,
    )
    poly!(
        plot,
        plot.cell_poly,
        strokecolor=plot.cellcolor,
        strokewidth=plot.celllinewidth,
        alpha=plot.cellalpha,
        color=colorant"transparent", # Transparent fill
    )
end

Makie.preferred_axis_type(::Atoms3D) = Makie.Axis3

@recipe Atoms2D (structure,) begin
    """
    Dictionary{Symbol, Float64} mapping element symbols to atomic radii in Angstroms for plotting. Default values are based on empirical atomic radii in picometers from Wikipedia data page.

    e.g. `:C => 0.7` for carbon atoms with a radius of 0.7 Å.
    """
    atomicradii = Dict([el => ustrip(uconvert(u"Å", rad * u"pm")) for (el, rad) in default_atomic_radii_pm]...)
    """
    Dictionary{Symbol, Colorant} mapping element symbols to colors for plotting. Default values are based on CPK coloring scheme taken from PeriodicTable.jl.

        e.g. `:C => colorant"black"` for carbon atoms colored black.
    """
    atomcolors = default_atom_fills
    """
    Colour multiplier used to darken the atom edge colour relative to the fill colour. Default value is 0.5, which means the edge colour will be half as bright as the fill colour.
    """
    atomsstrokecolormultiplier = 0.5
    """
    The marker to use for plotting atoms. Default is a high-resolution sphere, which might make vector graphics larger.
    Consider using `rasterize=true` in the plot attributes to mitigate this.
    """
    marker = Circle
    "Atom edge width in pt, default is 1.5 pt."
    strokewidth = 1.5
    "Whether to show cell boundaries. Default is true."
    show_cell = true
    "Colour for the cell boundaries. Default is black."
    cellcolor = colorant"black"
    "Line width for the cell boundaries in pt. Default is 1.0 pt."
    celllinewidth = 1.0
    "Alpha transparency for the cell boundaries. Default is 0.5."
    cellalpha = 1.0
end

function Makie.plot!(plot::Atoms2D{<:Tuple{<:NQCBase.Structure}},)
    input_nodes = [
        :structure,
        :atomicradii,
        :atomcolors,
        :show_cell,
        :atomsstrokecolormultiplier,
    ] # Take a single structure in
    output_nodes = [
        :positions, # Convert structure positions to plottable positions in Å.
        :markersize, # Generate marker sizes based on conversion dict
        :color, # Generate colors based on conversion dict
        :strokecolor, # Generate stroke colors based on conversion dict and multiplier
    ]
    register_computation!(plot.attributes, input_nodes, output_nodes) do inputs, changed, cached
        # Unpack inputs
        atomic_structure = inputs[:structure]
        atomicradii = inputs[:atomicradii]
        atomcolors = inputs[:atomcolors]
        multiplier = inputs[:atomsstrokecolormultiplier]

        # Convert positions to Å and then to Point3f format for plotting
        positions_conv = auconvert.(u"Å", atomic_structure.positions) .|> ustrip
        positions_AA::Vector{Point3f} = [Point3f(pos...) for pos in eachcol(positions_conv)]
        # Generate marker sizes based on conversion dict
        marker_sizes = [atomicradii[el] for el in atomic_structure.atoms.types]
        # Generate colors based on conversion dict
        color = [atomcolors[el] for el in atomic_structure.atoms.types]
        edge_colors = [atomcolors[el] * multiplier for el in atomic_structure.atoms.types]
        return (positions_AA, marker_sizes, color, edge_colors)
    end
    register_computation!(
        plot.attributes,
        [:structure, :show_cell],
        [:cell_poly],
    ) do inputs, changed, cached
        show_cell = inputs[:show_cell]
        atomic_structure = inputs[:structure]
        if show_cell || isa(atomic_structure.cell, PeriodicCell)
            # Generate Cell polygon from structure.cell.vectors if show_cell is true
            cell_poly = Point3f[]
            cell_vectors = auconvert.(u"Å", atomic_structure.cell.vectors) |> ustrip # Convert to Å
            # Convert unit cell into vertices and indices for mesh plotting. Each column of cell vectors is the unit vector in one dimension.
            poly_points = [
                zero(cell_vectors[:, 1]), # Origin
                cell_vectors[:, 1],
                cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 2],
                zero(cell_vectors[:, 1]), # Bottom face done
                cell_vectors[:, 3],
                cell_vectors[:, 3] + cell_vectors[:, 1],
                cell_vectors[:, 3] + cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 3] + cell_vectors[:, 2],
                cell_vectors[:, 3], # Top face done
                cell_vectors[:, 3] + cell_vectors[:, 1],
                cell_vectors[:, 1],
                cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 3] + cell_vectors[:, 1] + cell_vectors[:, 2],
                cell_vectors[:, 3] + cell_vectors[:, 2],
                cell_vectors[:, 2],
            ]
            for vertex in poly_points
                push!(cell_poly, Point3f(vertex...))
            end
            return (cell_poly,)
        end
        return (Point3f[],)
    end

    # Generate positions in Angstroms and convert to Point3f format for plotting
    scatter!(
        plot,
        plot.positions,
        marker=plot[:marker],
        markersize=plot.markersize,
        color=plot.color,
        markerspace=:data,
        #ssao = true,
        strokewidth=plot.strokewidth,
        strokecolor=plot.strokecolor,
        rasterize=true,
    )
    poly!(
        plot,
        plot.cell_poly,
        strokecolor=plot.cellcolor,
        strokewidth=plot.celllinewidth,
        alpha=plot.cellalpha,
        color=colorant"transparent", # Transparent fill
    )
end

Makie.preferred_axis_type(::Atoms2D) = Makie.Axis

end

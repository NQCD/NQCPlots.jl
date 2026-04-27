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


#atom_colors = [Symbol(el.symbol) => parse(Makie.Colors.Colorant, el.cpk_hex) for el in PeriodicTable.elements]
atom_colors = []
for el in PeriodicTable.elements
    col = colorant"#6e6e6e"
    try
        col = parse(Makie.Colors.Colorant, el.cpk_hex)
    catch e
    end
    push!(atom_colors, Symbol(el.symbol) => col)
end

export atomic_structures_theme

const atomic_structures_theme = Theme(
    Axis=( # Fit axes to data by default, so atoms are circular.
        aspect=Makie.DataAspect(),
    ),
    Axis3=( # Fit axes to data by default, so atoms are spherical.
        aspect=:data,
    ),
)

# Default atom colours
const default_atom_fills = Dict(atom_colors...)
# Empirical atomic radii (pm) from Wikipedia data page:
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
const default_atomic_radii_pm = Dict{Symbol,Float64}(
    :H => 25,
    :He => 120,
    :Li => 145,
    :Be => 105,
    :B => 85,
    :C => 70,
    :N => 65,
    :O => 60,
    :F => 50,
    :Ne => 160,
    :Na => 180,
    :Mg => 150,
    :Al => 125,
    :Si => 110,
    :P => 100,
    :S => 100,
    :Cl => 100,
    :Ar => 71,
    :K => 220,
    :Ca => 180,
    :Sc => 160,
    :Ti => 140,
    :V => 135,
    :Cr => 140,
    :Mn => 140,
    :Fe => 140,
    :Co => 135,
    :Ni => 135,
    :Cu => 135,
    :Zn => 135,
    :Ga => 130,
    :Ge => 125,
    :As => 115,
    :Se => 115,
    :Br => 115,
    :Rb => 235,
    :Sr => 200,
    :Y => 180,
    :Zr => 155,
    :Nb => 145,
    :Mo => 145,
    :Tc => 135,
    :Ru => 130,
    :Rh => 135,
    :Pd => 140,
    :Ag => 160,
    :Cd => 155,
    :In => 155,
    :Sn => 145,
    :Sb => 145,
    :Te => 140,
    :I => 140,
    :Cs => 260,
    :Ba => 215,
    :La => 195,
    :Ce => 185,
    :Pr => 185,
    :Nd => 185,
    :Pm => 185,
    :Sm => 185,
    :Eu => 185,
    :Gd => 180,
    :Tb => 175,
    :Dy => 175,
    :Ho => 175,
    :Er => 175,
    :Tm => 175,
    :Yb => 175,
    :Lu => 175,
    :Hf => 155,
    :Ta => 145,
    :W => 135,
    :Re => 135,
    :Os => 130,
    :Ir => 135,
    :Pt => 135,
    :Au => 135,
    :Hg => 150,
    :Tl => 190,
    :Pb => 180,
    :Bi => 160,
    :Po => 190,
    :Ra => 215,
    :Ac => 195,
    :Th => 180,
    :Pa => 180,
    :U => 175,
    :Np => 175,
    :Pu => 175,
    :Am => 175,
    :Cm => 176,
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

end

"""
A module containing Makie recipes for NQCBase structures. 

I wish this could've been an extension of NQCBase, but the Julia extension system doesn't allow for it. 
"""
module NQCPlots

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
    push!(atom_colors, Symbol(el) => col)
end

# Default atom colours
const default_atom_fills = Dict(atom_colors...)
# Empirical atomic radii (pm) from Wikipedia data page:
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
const default_atomic_radii_pm = Dict{Symbol, Float64}(
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
         Makie.GeometryBasics.Sphere(Point3f(0.0,0.0,0.0), 1.0),
         128
     )
 )

@recipe AtomicStructure3DAtoms (structure, ) begin
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
# Makie.args_preferred_axis(::Type{<: AtomicStructure3DAtoms}) = Makie.Axis3

function Makie.plot!(plot::AtomicStructure3DAtoms, )
    structure = plot[:structure]    
    # Generate positions in Angstroms and convert to Point3f format for plotting
    positions_AA = auconvert.(u"Å", structure.positions) .|> ustrip
    positions_AA = [Point3f(pos...) for pos in eachcol(positions_AA)]
    # Generate marker sizes based on conversion dict
    marker_sizes = [plot[:atomicradii][el] for el in structure.atoms.symbols]
    # Generate colors based on conversion dict
    colors = [plot[:atomcolors][el] for el in structure.atoms.symbols]
    meshscatter!(
          plot,
          positions_AA,
          color = colors,
          marker = plot[:marker],
          markersize = marker_sizes,
          #ssao = true,
          # strokewidth = atom_strokewidth,
          # strokecolor = [get(atom_edges, i, colorant"black") for i in trj_postprocessed[:atoms].types[atoms_to_draw]],
          rasterize = true,
        )
    return plot
end

end

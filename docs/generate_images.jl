# Script to generate example images for README documentation.
# Run from the NQCPlots.jl root directory:
#   julia --project=test docs/generate_images.jl
# GLMakie is used instead of CairoMakie to get correct draw order for 3D meshscatter.

using GLMakie
GLMakie.activate!(; visible=false)  # Headless/offscreen rendering
using NQCBase
import NQCBase: Structure, Atoms, InfiniteCell, PeriodicCell
using NQCPlots
using Unitful
using UnitfulAtomic

# Output directory for images
img_dir = joinpath(@__DIR__, "img")
mkpath(img_dir)

# ─── Helper: save a figure ────────────────────────────────────────────────────
function save_fig(name, fig)
    path = joinpath(img_dir, name)
    save(path, fig; px_per_unit=2)
    println("Saved: $path")
end

# ─── 1. Build a simple water molecule (no periodic cell) ─────────────────────
# Positions are in atomic units (Bohr).
# O–H bond length ≈ 1.81 Bohr, H–O–H angle ≈ 104.5°
function water_structure()
    atoms = Atoms([:O, :H, :H])
    # O at origin, H atoms placed at ~0.96 Å = 1.81 Bohr
    bond_au = austrip(0.96u"Å")
    angle = deg2rad(104.5 / 2)
    positions = [
        0.0 bond_au*sin(angle) -bond_au*sin(angle);
        0.0 bond_au*cos(angle) bond_au*cos(angle);
        0.0 0.0 0.0
    ]
    return Structure(atoms, positions, InfiniteCell())
end

# ─── 2. Build a small FCC Cu(111) slab with one CO adsorbate ─────────────────
function cu_co_slab()
    # FCC Cu lattice constant ≈ 3.61 Å
    a = austrip(3.61u"Å")
    # (111) surface vectors
    a1 = [a, 0.0, 0.0]
    a2 = [a / 2, a * sqrt(3) / 2, 0.0]
    a3 = [0.0, 0.0, a * sqrt(6) / 3 + 15.0]  # Add vacuum of ~15 Å above slab
    cell_mat = hcat(a1, a2, a3)
    cell = PeriodicCell(cell_mat)

    layer_z_spacing = austrip(2.08u"Å")   # ≈ (111) interlayer distance

    # 3-layer 2×2 slab
    cu_types = fill(:Cu, 12)
    atoms = Atoms(vcat(cu_types, [:C, :O]))

    xs_base = [0.0, a / 2, a / 4, 3a / 4]
    ys_base = [0.0, 0.0, a * sqrt(3) / 4, a * sqrt(3) / 4]

    positions = Matrix{Float64}(undef, 3, 14)
    idx = 1
    for layer in 0:2
        for k in 1:4
            positions[:, idx] = [xs_base[k], ys_base[k], layer * layer_z_spacing]
            idx += 1
        end
    end
    # CO on top site above first Cu of top layer
    c_z = 3 * layer_z_spacing + austrip(1.3u"Å")
    o_z = c_z + austrip(1.15u"Å")
    positions[:, 13] = [xs_base[1], ys_base[1], c_z]
    positions[:, 14] = [xs_base[1], ys_base[1], o_z]

    return Structure(atoms, positions, cell)
end

# ─── 3. Build a methane molecule ─────────────────────────────────────────────
function methane_structure()
    atoms = Atoms([:C, :H, :H, :H, :H])
    d = austrip(1.089u"Å")   # C–H bond length
    s = d / sqrt(3)
    positions = [
        0.0 s -s s -s;
        0.0 s s -s -s;
        0.0 s -s -s s
    ]
    return Structure(atoms, positions, InfiniteCell())
end

# ─── 4. Build a NaCl unit cell ───────────────────────────────────────────────
function nacl_structure()
    a = austrip(5.64u"Å")   # NaCl lattice constant
    cell = PeriodicCell([a 0 0; 0 a 0; 0 0 a] .* 1.0)
    atoms = Atoms([:Na, :Cl, :Na, :Cl, :Na, :Cl, :Na, :Cl])
    h = a / 2
    positions = [
        0.0 h h 0.0 0.0 h h 0.0;
        0.0 0.0 h h 0.0 0.0 h h;
        0.0 0.0 0.0 0.0 h h h h
    ]
    return Structure(atoms, positions, cell)
end

# ════════════════════════════════════════════════════════════════════════════
# Image 1 — Default plot of a water molecule
# ════════════════════════════════════════════════════════════════════════════
let
    water = water_structure()
    with_theme(atomic_structures_theme) do
        fig, ax, plt = atoms3d(water)
        ax.title = "Water molecule (default)"
        save_fig("water_default.png", fig)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# Image 2 — Methane with default settings
# ════════════════════════════════════════════════════════════════════════════
let
    ch4 = methane_structure()
    with_theme(atomic_structures_theme) do
        fig, ax, plt = atoms3d(ch4)
        ax.title = "Methane (default)"
        save_fig("methane_default.png", fig)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# Image 3 — NaCl with cell visible (default) vs. hidden
# ════════════════════════════════════════════════════════════════════════════
let
    nacl = nacl_structure()
    with_theme(atomic_structures_theme) do
        fig = Figure(size=(900, 400))
        ax1 = Axis3(fig[1, 1]; aspect=:data, title="show_cell = true (default)")
        ax2 = Axis3(fig[1, 2]; aspect=:data, title="show_cell = false")
        atoms3d!(ax1, nacl)
        atoms3d!(ax2, nacl; show_cell=false)
        save_fig("nacl_cell_options.png", fig)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# Image 4 — Cu/CO slab: default, custom cell colour, custom cell linewidth
# ════════════════════════════════════════════════════════════════════════════
let
    slab = cu_co_slab()
    with_theme(atomic_structures_theme) do
        fig = Figure(size=(1200, 400))
        ax1 = Axis3(fig[1, 1]; aspect=:data, title="Default")
        ax2 = Axis3(fig[1, 2]; aspect=:data, title="cellcolor=:royalblue,\ncellalpha=0.8")
        ax3 = Axis3(fig[1, 3]; aspect=:data, title="show_cell=false")
        atoms3d!(ax1, slab)
        atoms3d!(ax2, slab; cellcolor=:royalblue, cellalpha=0.8, celllinewidth=2.0)
        atoms3d!(ax3, slab; show_cell=false)
        save_fig("slab_cell_options.png", fig)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# Image 5 — Custom atomic radii and colours on methane
# ════════════════════════════════════════════════════════════════════════════
let
    ch4 = methane_structure()
    with_theme(atomic_structures_theme) do
        fig = Figure(size=(900, 400))
        ax1 = Axis3(fig[1, 1]; aspect=:data, title="Default radii & colours")
        ax2 = Axis3(fig[1, 2]; aspect=:data, title="Custom radii & colours")
        atoms3d!(ax1, ch4)
        atoms3d!(ax2, ch4;
            atomicradii=Dict(:C => 0.40, :H => 0.25),
            atomcolors=Dict(:C => colorant"#b4da55", :H => colorant"#0065bd"),
        )
        save_fig("methane_custom_radii_colors.png", fig)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# Image 6 — strokewidth / atomsstrokecolormultiplier effect
# ════════════════════════════════════════════════════════════════════════════
let
    water = water_structure()
    with_theme(atomic_structures_theme) do
        fig = Figure(size=(900, 400))
        ax1 = Axis3(fig[1, 1]; aspect=:data, title="Default stroke (1.5 pt)")
        ax2 = Axis3(fig[1, 2]; aspect=:data, title="strokewidth = 4")
        atoms3d!(ax1, water)
        atoms3d!(ax2, water; strokewidth=4)
        save_fig("water_strokewidth.png", fig)
    end
end

println("\nAll images written to: $img_dir")

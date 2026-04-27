using NQCBase, NQCPlots, CairoMakie
structure = structure = read_extxyz("test/Cu111_T-200K_3x3x6_highH_relax.xyz") |> first
f, ax, pl = atoms3d(structure, axis=(title="Atoms3D Test", aspect=:data))

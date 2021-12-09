from mesh_parametrization import mesh_deformation as md

# import mesh
mesh = md.read_mesh("LProfilTetMsh.xdmf")

# import nodes and displacements
nodes, displacements = md.read_input_file("example_input.txt")

md.deform_mesh(mesh, nodes, displacements)

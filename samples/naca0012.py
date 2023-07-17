import hexed

solver = hexed.cpp.make_solver(2, 6, 1.)
solver.mesh().add_tree([hexed.cpp.new_copy(hexed.cpp.Nonpenetration()) for i in range(4)])
solver.visualize_field_tecplot("hexed_out/new_mesh")

"""
solver = hexed.Solver(2, 6)
solver.generate_mesh([-10, -10], [10, 10], geometries = [hexed.naca("0012", int(1e4))], surf_resolution = .1, refine_sweeps = 8)
solver.visualize_field("mesh")
"""

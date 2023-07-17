import hexed

solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], geometries = ['hexed.naca("0012", int(1e4))'], surf_resolution = .1, refine_sweeps = 8)
solver.visualize_field_tecplot("hexed_out/new_mesh")

"""
solver = hexed.Solver(2, 6)
solver.generate_mesh([-10, -10], [10, 10], geometries = [hexed.naca("0012", int(1e4))], surf_resolution = .1, refine_sweeps = 8)
solver.visualize_field("mesh")
"""

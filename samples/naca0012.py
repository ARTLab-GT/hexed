import hexed

solver = hexed.Solver(2, 6)
solver.generate_mesh([-10, -10], [10, 10], geometries = [hexed.naca("0012")], refine_sweeps = 8)
solver.visualize_field("mesh")

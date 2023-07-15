import hexed

solver = hexed.Solver(2, 6)
solver.generate_mesh([-10, -10], [10, 10])
solver.visualize_field("mesh")

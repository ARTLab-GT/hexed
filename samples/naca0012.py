import hexed

solver = hexed.Solver(2, 6, 1.)
solver.mesh().add_tree()
solver.visualize("test")

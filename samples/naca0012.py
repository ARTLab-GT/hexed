import hexed

solver = hexed.Solver(2, 6, 1.)
solver.mesh().add_tree()
#solver.visualize("test")
print(hexed.matrix_shape.__doc__)

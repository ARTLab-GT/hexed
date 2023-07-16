import hexed
import matplotlib.pyplot as plt

geom = hexed.naca("0012", int(1e4))
plt.plot(geom[:, 0], geom[:, 1])
plt.axis("equal")
plt.show()
solver = hexed.Solver(2, 6)
solver.generate_mesh([-10, -10], [10, 10], geometries = [geom], refine_sweeps = 8)
solver.visualize_field("mesh")

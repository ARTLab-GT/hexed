import hexed
import numpy as np

solver = hexed.Solver(2, 6)
angle = np.linspace(0, 2*np.pi, 100)
solver.generate_mesh([-10, -10], [10, 10], geometries = [np.array([np.cos(angle), np.sin(angle)]).transpose()])
solver.visualize_field("mesh")

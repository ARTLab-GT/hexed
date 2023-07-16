import hexed
import numpy as np
import matplotlib.pyplot as plt

solver = hexed.Solver(2, 6)
solver.generate_mesh([-10, -10], [10, 10], geometries = [hexed.naca("0012")])
solver.visualize_field("mesh")

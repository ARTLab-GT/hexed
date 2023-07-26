import hexed
import os

solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], freestream = hexed.flow_state(altitude = 0, velocity = [0., 0.]),
                             geom = hexed.make_geom("samples/crmhl-2dcut.igs", n_dim = 2),
                             final_resolution = "res_bad > 3e-3 && ref_level < 13", final_max_iters = 20)
print(f"n elements: {solver.mesh().n_elements()}")
solver.visualize()

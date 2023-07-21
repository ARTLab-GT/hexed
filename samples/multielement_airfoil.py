import hexed
import os

solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], freestream = hexed.flow_state(altitude = 0, velocity = [0., 0.]),
                             #geometries = ["samples/crmhl-2dcut.igs"],
                             geometries = ["samples/" + file_name for file_name in os.listdir("samples") if ".txt" in file_name],
                             final_resolution = "res_bad > 3e-3 && ref_level < 13", final_max_iters = 30)
solver.visualize()

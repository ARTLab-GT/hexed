import hexed

"""
solver = hexed.create_solver(3, 3, [-.1, 0., -.3], [.5, .6, .3], freestream = hexed.flow_state(altitude = 0., velocity = [0., 0., 0.]),
                             geom = hexed.make_geom("samples/wing_store.IGS", n_dim = 3),
                             final_resolution = "res_bad > 3e-2 && ref_level < 8", final_max_iters = 8, resolution_metric = "continuity",
                             n_smooth = 10)
"""
solver = hexed.create_solver(3, 5, [-.1, 0., -.3], [.5, .6, .3], freestream = hexed.flow_state(altitude = 0., velocity = [0., 0., 0.]),
                             geom = hexed.make_geom("samples/wing_store.IGS", n_dim = 3),
                             final_resolution = "res_bad > 1e-2 && ref_level < 8", final_max_iters = 8,
                             n_smooth = 15)
print(f"n elements: {solver.mesh().n_elements()}")
print("visualizing")
solver.visualize_surface_tecplot(solver.mesh().surface_bc_sn(), "hexed_out/mesh", 10)
print("Geometry processing performance: " + hexed.cpp.Simplex_geom[3].performance_report())

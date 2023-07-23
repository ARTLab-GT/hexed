import hexed

solver = hexed.create_solver(3, 3, [-.1, 0., -.3], [.5, .6, .3], freestream = hexed.flow_state(altitude = 0., velocity = [0., 0., 0.]),
                             #geom = hexed.make_geom("samples/wing_store.IGS", n_dim = 3),
                             geom = hexed.make_geom("samples/wing_store.STL"),
                             final_resolution = "res_bad > 3e-2 && ref_level < 8", final_max_iters = 3, resolution_metric = "continuity",
                             n_smooth = 20)
solver.visualize_surface_tecplot(solver.mesh().surface_bc_sn(), "hexed_out/mesh", 3)
print("Geometry processing performance: " + hexed.cpp.Simplex_geom[3].performance_report())

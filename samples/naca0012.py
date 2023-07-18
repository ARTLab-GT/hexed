import hexed

solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], hexed.flow_state(altitude = 0, mach = .2, attack = 3*hexed.cpp.degree),
                             geometries = [hexed.naca("0012", int(1e4))], final_resolution = "is_extruded && ref_level < 8")
solver.visualize_field_tecplot("hexed_out/new_mesh")

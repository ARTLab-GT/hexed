import hexed

print("go")
solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], freestream = hexed.flow_state(altitude = 0, mach = 0.2, attack = 3*hexed.cpp.degree),
                             geometries = [hexed.naca("0012", int(1e4))], final_resolution = "is_extruded && ref_level < 8")
solver.visualize_field_tecplot("hexed_out/initial")

@solver.method
def done(self):
    return self.iteration_status().wall_time() > 600

solver.run()

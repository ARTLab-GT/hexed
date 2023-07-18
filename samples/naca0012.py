import hexed

solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], hexed.flow_state(altitude = 0, mach = .2, attack = 3*hexed.cpp.degree),
                             geometries = [hexed.naca("0012", int(1e4))], final_resolution = "is_extruded && ref_level < 8")
solver.visualize_field_tecplot("hexed_out/initial")
print(solver.iteration_status().header())
for i in range(10):
    for j in range(10):
        for k in range(10):
            solver.update()
        status = solver.iteration_status()
        print(status.report())
    solver.visualize_field_tecplot(f"hexed_out/iter{status.iteration}")

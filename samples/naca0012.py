import hexed

solver = hexed.create_solver(2, 6, [-10, -10], [10, 10], freestream = hexed.flow_state(altitude = 0, mach = 0.2, attack = 3*hexed.cpp.degree),
                             geometries = [hexed.naca("0012", int(1e4))], final_resolution = "is_extruded && ref_level < 9")
solver.n_step = 1000
solver.vis_freq = None

@solver.method
def done(self):
    return self.iteration == 10000

stopwatch = hexed.cpp.Stopwatch()
stopwatch.start()
solver.run()
stopwatch.pause()
print(stopwatch.time())
print(solver.stopwatch_tree().report())

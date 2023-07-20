import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import hexed

fig = plt.figure()
xdata = [0]
ydata = [1]
lines = plt.plot(xdata, ydata)
def init():
    plt.xlim(0, 1)
    plt.ylim(1, 2)
    plt.yscale("log")
def update(frame):
    time.sleep(2)
    xdata.append(xdata[-1] + 1)
    ydata.append((xdata[-1] + 1)**2)
    lines[0].set_data(xdata, ydata)
    if plt.xlim()[1] < xdata[-1] or plt.ylim()[1] < ydata[-1]:
        plt.xlim(0, 2*xdata[-1])
        plt.ylim(1, ydata[-1]**2)
    return lines
animation = ani.FuncAnimation(fig, update, init_func = init, cache_frame_data = False)
"""
mon = hexed.Convergence_monitor()
mon.push(0, 1.)
time.sleep(2)
mon.push(1, .5)
time.sleep(2)
mon.push(2, .25)
time.sleep(2)
"""
plt.show()
exit()

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

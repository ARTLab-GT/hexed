import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

exper = np.loadtxt("experiment.csv", delimiter = ",")
#plt.scatter(exper[:, 0], exper[:, 1])
exper = np.loadtxt("tewfik.csv", delimiter = ",")
plt.scatter(exper[:, 0], exper[:, 1]/exper[0, 1])
analy = np.loadtxt("analytic.csv", delimiter = ",")
#plt.plot(analy[:, 0], analy[:, 1], linestyle = "--", color = "k")
for suffix in ["8"]:
    sim = pd.read_csv(f"hexed_out_isothermal{suffix}/simulation.csv")
    sim.sort_values("angle", inplace = True)
    max_angle = 80
    sim = sim.loc[sim["angle"] < max_angle*np.pi/180]
    plt.plot(sim["angle"]*180/np.pi, sim["heat_flux"]/sim["heat_flux"][0])
plt.show()

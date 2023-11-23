import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

max_angle = 80
exper = np.loadtxt("tewfik.csv", delimiter = ",")
exper = exper[exper[:, 0] < max_angle + 3]
plt.scatter(exper[:, 0], exper[:, 1]/exper[0, 1], label = "experiment")
sim = pd.read_csv(f"hexed_out_isothermal8/simulation.csv")
sim.sort_values("angle", inplace = True)
sim = sim.loc[sim["angle"] < max_angle*np.pi/180]
ad = pd.read_csv(f"hexed_out_ad8/simulation.csv")
ad.sort_values("angle", inplace = True)
ad = ad.loc[ad["angle"] < max_angle*np.pi/180]
htc = np.array(sim["heat_flux"])/(np.array(ad["temperature"]) - np.array(sim["temperature"]))
plt.plot(sim["angle"]*180/np.pi, htc/htc[0], label = "simulation")
plt.xlabel(r"$\theta$ ($\degree$)")
plt.ylabel(r"$c_h(\theta)/c_h(0)$")
plt.grid(True)
plt.legend()
plt.show()

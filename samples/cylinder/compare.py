import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

max_angle = 70
experiment = np.loadtxt("experiment_nusselt.csv", delimiter = ",")
recovery = np.loadtxt("experiment_recovery.csv", delimiter = ",")
experiment = experiment[experiment[:, 0] <= max_angle, :]
plt.scatter(experiment[:, 0], experiment[:, 1]*1e3, label = "experiment")

isothermal = pd.read_csv("hexed_out_isothermal/surface_data.csv")
isothermal.sort_values("angle", inplace = True)
isothermal = isothermal.loc[isothermal["angle"] < max_angle*np.pi/180]
adiabatic = pd.read_csv("hexed_out_adiabatic/surface_data.csv")
adiabatic.sort_values("angle", inplace = True)
adiabatic = adiabatic.loc[adiabatic["angle"] < max_angle*np.pi/180]
angle = isothermal["angle"]*180/np.pi
nusselt = np.array(isothermal["heat_flux"])/(np.array(adiabatic["temperature"]) - 329.631)*2/(0.00776826*0.82645/0.7)
plt.plot(angle, nusselt, label = "simulation (constant Pr)")
nusselt = np.array(isothermal["heat_flux"])/(np.array(adiabatic["temperature"]) - 329.631)*2/0.00776826
plt.plot(angle, nusselt, label = "simulation (non-constant Pr)")

plt.xlabel(r"angle ($\degree$)")
plt.ylabel("Nusselt number")
plt.legend()
plt.grid(True)
plt.show()

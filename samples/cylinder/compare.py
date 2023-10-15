import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

max_angle = 80
experiment = np.loadtxt("experiment_nusselt.csv", delimiter = ",")
recovery = np.loadtxt("experiment_recovery.csv", delimiter = ",")
experiment = experiment[experiment[:, 0] <= max_angle, :]
recovery = recovery[recovery[:, 0] <= max_angle, :]

stag_temp = 394.261
cases = []
temps = []
fluxes = []
for i in range(2):
    cases.append(pd.read_csv(f"hexed_out_temperature{i}/surface_data.csv"))
    cases[-1].sort_values("angle", inplace = True)
    cases[-1] = cases[-1].loc[cases[-1]["angle"] < max_angle*np.pi/180]
    temps.append(np.array(cases[-1]["temperature"]))
    fluxes.append(np.array(cases[-1]["heat_flux"]))
angle = cases[0]["angle"]*180/np.pi
recovery_temp = temps[0] - fluxes[0]*(temps[0] - temps[1])/(fluxes[0] - fluxes[1])

plt.scatter(recovery[:, 0], recovery[:, 1], label = "experiment")
plt.plot(angle, recovery_temp/stag_temp, label = "simulation")
plt.xlabel(r"angle ($\degree$)")
plt.ylabel("$T_r/T_0$")
plt.legend()
plt.grid(True)
plt.figure()

plt.scatter(experiment[:, 0], experiment[:, 1]*1e3, label = "experiment")
nusselt = fluxes[0]/(recovery_temp - temps[0])*2/0.00776826
plt.plot(angle, nusselt, label = "simulation")
plt.xlabel(r"angle ($\degree$)")
plt.ylabel("Nusselt number")
plt.legend()
plt.grid(True)
plt.show()

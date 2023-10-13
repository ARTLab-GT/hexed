import matplotlib.pyplot as plt
import numpy as np

experiment = np.loadtxt("experiment_nusselt.csv", delimiter = ",")
recovery = np.loadtxt("experiment_recovery.csv", delimiter = ",")
plt.scatter(experiment[:, 0], experiment[:, 1], label="experiment")
plt.xlabel(r"angle ($\degree$)")
plt.ylabel("Nusselt number")
plt.legend()
plt.grid(True)
plt.show()

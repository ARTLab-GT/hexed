import numpy as np

angle = np.linspace(0, 2*np.pi, 10**4)
coords = np.array([np.sin(angle), np.cos(angle)]).transpose()*.127/2
np.savetxt("auto.csv", coords, delimiter=", ")

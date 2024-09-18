import numpy as np

angle = np.linspace(.5*np.pi, np.pi, 10**4)
coords = np.array([np.cos(angle), np.sin(angle)]).transpose()*.02
np.savetxt("auto.csv", coords, delimiter=", ")

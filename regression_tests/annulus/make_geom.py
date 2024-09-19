import numpy as np

angle = np.linspace(0., 2*np.pi, 10**4)
radius = [1., .25]
for i in range(2):
    coords = radius[i]*np.array([np.cos(angle), np.sin(angle)]).transpose()
    np.savetxt(f"auto{i}.csv", coords, delimiter=", ")

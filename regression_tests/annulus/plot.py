import matplotlib.pyplot as plt
import numpy as np

def plot(size, error, order):
    plt.scatter(size, error, label = f"p = {order}")
    coefs = np.polyfit(np.log(2.**np.array(size)), np.log(error), 1)
    print(f"observed order (p = {order}): {-coefs[0]}")

plot([3, 4, 5], [0.000286724, 3.52875e-05, 3.41315e-06], 3)
plot([3, 4], [6.55922e-06, 2.90258e-07], 5)
plt.xlabel("refinement level")
plt.ylabel("$L_2$ error")
plt.legend()
plt.yscale("log")
plt.grid(True)
plt.show()

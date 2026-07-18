import numpy as np
import matplotlib.pyplot as plt

x, y, r = np.loadtxt("points.dat", unpack=True)

plt.scatter(x, y, c=r, s=1, cmap="tab10")
plt.axis("equal")
plt.show()

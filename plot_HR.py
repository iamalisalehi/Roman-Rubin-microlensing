import numpy as np 
import matplotlib.pyplot as plt 

x, thin, bulge, thick, halo, rho, N = np.loadtxt("./files/density/D1.dat", unpack="True")

plt.plot(x, thin, label="Thin Disk")
plt.plot(x, bulge, label="Bulge")
plt.plot(x, thick, label="Thick Disk")
plt.plot(x, halo, label="Halo")
plt.legend()
plt.show()


plt.plot(x, rho, label=r"$\rho_\star$")
plt.legend()
plt.show()

plt.plot(x, N, label=r"$N_\star$")
plt.legend()
plt.show()

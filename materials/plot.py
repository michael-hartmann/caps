import numpy as np
import matplotlib.pyplot as plt

hbar = 6.582119514e-16 # eV*s/rad

gold   = np.loadtxt("gold.csv")

plt.xlabel("ξ in rad/s")
plt.ylabel("ε(iξ)")
plt.loglog(gold[:,0], gold[:,1])

plt.show()

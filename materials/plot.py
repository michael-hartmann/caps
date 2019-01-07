import numpy as np
import matplotlib.pyplot as plt

hbar = 6.582119514e-16 # eV*s/rad

dalvit = np.loadtxt("GoldDalvit.dat")
gold   = np.loadtxt("GoldPalik.dat")
palik  = np.loadtxt("GoldPalik.dat")

plt.xlabel("ξ in rad/s")
plt.ylabel("ε(iξ)")
plt.loglog(gold[:,0], gold[:,1])
plt.loglog(dalvit[:,0], dalvit[:,1])
plt.loglog(palik[:,0], palik[:,1])

plt.show()

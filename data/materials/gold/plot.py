from math import log,e
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

hbar = 6.582119514e-16 # eV/s

def epsilon_drude(xi, omegap, gamma):
    """Fit function for Drude model"""
    return 1+omegap**2/(xi*(xi+gamma))


# read data of gold
data = np.loadtxt("GoldEpsIm.dat")
xi      = data[:,0]
epsilon = data[:,1]

# fit small
npts = 50 # use first npts points for fit
xdata = xi[:npts]
ydata = epsilon[:npts]
popt, pcov = curve_fit(epsilon_drude, xdata, ydata, p0=(9/hbar, 0.035/hbar))
omega_p, gamma = popt

# fit small
npts = 50 # use first npts points for fit
xdata = xi[-npts:]
ydata = epsilon[-npts:]
popt, pcov = curve_fit(epsilon_drude, xdata, ydata, p0=(9/hbar, 0.035/hbar))
omega_p2, gamma2 = popt


print("xi small")
print("ωp = %.8g eV" % (omega_p*hbar))
print("γ  = %.8g eV" % (gamma*hbar))
print()

print("xi large")
print("ωp = %.8g eV" % (omega_p2*hbar))
print("γ  = %.8g eV" % (gamma2*hbar))


# plot
xi = np.logspace(log(1e11), log(1e19), 1000, base=e)
drude = (omega_p)**2/(xi*(xi+gamma))

drude2 = (omega_p2)**2/(xi*(xi+gamma2))

plt.plot(data[:,0], data[:,1]-1)
plt.plot(xi, drude)
plt.plot(xi, drude2)
plt.xlabel("xi")
plt.ylabel("epsilon(i*xi)")
plt.xscale('log')
plt.yscale('log')
plt.show()

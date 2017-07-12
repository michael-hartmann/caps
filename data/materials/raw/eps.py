import numpy as np
from math import atan, pi, log, e
from scipy.integrate import simps
import matplotlib.pyplot as plt

class Gold:
    def __init__(self, omega, eps2, omegap=8.45, gamma=0.047):
        self.omega  = omega
        self.eps2   = eps2
        self.omegap = omegap
        self.gamma  = gamma
        self.omegac = omega[-1]

    def epsm1(self, xi):
        xi  = np.array(xi)
        eps = np.empty(xi.shape)

        gamma = self.gamma
        omegap = self.omegap
        omegac = self.omegac

        for i,xi_i in enumerate(xi):
            # drude
            I1 = omegap**2/(xi_i*(xi_i**2-gamma**2))*( xi_i*atan(omegac/gamma) - gamma*atan(omegac/xi_i) )
            
            integrand = self.omega*self.eps2/(self.omega**2+xi_i**2)
            I2 = -simps(integrand, self.omega)

            eps[i] = 2/pi*(I1+I2)

        return eps


if __name__ == "__main__":
    # olmon_EV
    data_olmonEV = np.loadtxt("Olmon_PRB2012_EV.dat", delimiter=",", usecols=(0,3))
    omega_olmonEV = data_olmonEV[:,0]
    eps2_olmonEV = data_olmonEV[:,1]
    gold_olmonEV = Gold(omega_olmonEV, eps2_olmonEV, omegap=8.45, gamma=0.047)

    # olmon_SC
    data_olmonSC = np.loadtxt("Olmon_PRB2012_SC.dat", delimiter=",", usecols=(0,3))
    omega_olmonSC = data_olmonSC[:,0]
    eps2_olmonSC = data_olmonSC[:,1]
    gold_olmonSC = Gold(omega_olmonSC, eps2_olmonSC, omegap=8.45, gamma=0.047)

    # olmon_TS
    data_olmonTS = np.loadtxt("Olmon_PRB2012_TS.dat", delimiter=",", usecols=(0,3))
    omega_olmonTS = data_olmonTS[:,0]
    eps2_olmonTS = data_olmonTS[:,1]
    gold_olmonTS = Gold(omega_olmonTS, eps2_olmonTS, omegap=8.45, gamma=0.047)

    # palik
    data_palik = np.loadtxt("palik.csv", delimiter=",", usecols=(0,3))
    omega_palik = data_palik[:,0]
    eps2_palik = data_palik[:,1]
    gold_palik = Gold(omega_palik, eps2_palik, omegap=9, gamma=0.0342)

    # brandli
    data_brandli = np.loadtxt("brandli.csv", delimiter=",", usecols=(0,3))
    omega_brandli = data_brandli[:,0]
    eps2_brandli = data_brandli[:,1]

    # dft
    data_dft = np.loadtxt("Werner.csv", delimiter=",", usecols=(0,3))
    omega_dft = data_dft[:,0]
    eps2_dft = data_dft[:,1]

    # Hagemann
    data_hagemann = np.loadtxt("Hagemann.csv", delimiter=",", usecols=(0,3))
    omega_hagemann = data_hagemann[:,0]*1e6
    eps2_hagemann = data_hagemann[:,1]

    # merged
    data_merged = np.loadtxt("mixed.csv", delimiter=",", usecols=(0,3))
    omega_merged = data_merged[:,0]
    eps2_merged = data_merged[:,1]
    gold_merged = Gold(omega_merged, eps2_merged, omegap=8.45, gamma=0.047)


    # olmon + palik
    omega = np.logspace(log(0.0001), log(1000), 1000, base=e)
    eps_olmon  = gold_olmonEV.epsm1(omega)
    eps_palik  = gold_palik.epsm1(omega)
    eps_merged = gold_merged.epsm1(omega)


    #plt.loglog(omega_palik, eps2_palik, "x", label="Palik")
    plt.loglog(omega_olmonEV, eps2_olmonEV, "o", label="Olmon EV")
    plt.loglog(omega_olmonSC, eps2_olmonSC, "x", label="Olmon SC")
    plt.loglog(omega_olmonTS, eps2_olmonTS, "d", label="Olmon TS")
    plt.loglog(omega_dft, eps2_dft, "+", label="DFT")
    plt.loglog(omega_palik, eps2_palik, "+", label="Palik")
    #plt.loglog(omega_hagemann, eps2_hagemann, "d", label="Hagemann")
    plt.xlabel("omega [eV]")
    plt.ylabel("eps''(omega)")
    plt.legend()
    plt.show()

    dalvit = np.loadtxt("../GoldEpsIm.dat")
    dalvit[:,0] /= 1.519e15
    dalvit[:,1] -= 1

    plt.loglog(omega, eps_palik, label="Palik")
    plt.loglog(omega, eps_merged, label="Olmon EV + Palik")
    plt.loglog(dalvit[:,0], dalvit[:,1], label="Diego")
    plt.xlabel("omega [eV]")
    plt.ylabel("eps(i*omega)-1")
    plt.legend()
    plt.show()

import numpy as np
from math import atan, pi, log10
from scipy.integrate import simps

hbar_eV = 6.582119514e-16 # eV s

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
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Rotate optical data to imaginary axis")

    # required argument
    parser.add_argument("filename", type=str, help="input filename")

    # optional arguments
    parser.add_argument("--start",  action="store", dest="start", type=float, default=1e11, help="omega start (in rad/s)")
    parser.add_argument("--stop",   action="store", dest="stop",  type=float, default=1e19, help="omega stop (in rad/s)")
    parser.add_argument("--points", action="store", dest="N",     type=int,   default=800,  help="number of points")

    parser.add_argument("--omegap_low",  action="store", dest="omegap_low", type=float, default=8.45,  help="omegap low (in eV)")
    parser.add_argument("--gamma_low",   action="store", dest="gamma_low",  type=float, default=0.047, help="gamma low (in eV)")

    parser.add_argument("--omegap_high", action="store", dest="omegap_high", type=float, default=0, help="omegap high (in eV)")
    parser.add_argument("--gamma_high",  action="store", dest="gamma_high",  type=float, default=0, help="gamma high (in eV)")

    args = parser.parse_args()

    data = np.loadtxt(args.filename, delimiter=",", usecols=(0,3))
    gold = Gold(data[:,0], data[:,1], omegap=args.omegap_low, gamma=args.gamma_low)

    omega = np.logspace(log10(args.start), log10(args.stop), args.N)
    eps2  = 1+gold.epsm1(omega*hbar_eV)

    print("# This file contains the dielectric function epsilon(i*xi) for Gold")
    print("# input filename: %s" % args.filename)
    print("#")
    print("# Drude parameters for low frequencies:")
    print("# omegap_low = %geV" % args.omegap_low)
    print("# gamma_low = %geV"  % args.gamma_low)
    print("#")
    print("# Drude parameters for high frequencies:")
    print("# omegap_high = %geV" % args.omegap_high)
    print("# gamma_high = %geV"  % args.gamma_high)
    print("#")
    print("# xi in rad/s            epsilon(i*xi)")

    for i,omega_ in enumerate(omega):
        print("%.6e\t%.6e" % (omega_, eps2[i]))

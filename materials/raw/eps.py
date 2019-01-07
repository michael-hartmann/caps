import numpy as np
from math import atan, pi, log10
from scipy.integrate import simps

class Material:
    """Compute dielectric function at the imaginary axis from optical data

    Compute the dielectric function at imaginary frequencies, i.e.
    epsilon(i*xi), from (measured) the imaginary part of the dielectric
    function at real frequencies.

    At real frequencies omega, the dielectric function is a complex valued function
        epsilon(omega) = eps1 + i*eps2
    where eps1 is the real part and eps2 is the imaginary part. Using the
    imaginary part eps2 of the dielectric function, it is possible to compute
    the dielectric function at imaginary frequencies, i.e., to compute
    epsilon(i*xi) where xi is real. For small frequencies there might be not
    enough data points. In this case an extrapolation assuming the Drude model
    is used.

    For integration the Simpson's rule is used.

    References:
    [1] Lambrecht, Reynaud, "Casimir force between metallic mirros", Eur. Phys.
        J. D 8, 309 (2000), https://doi.org/10.1007/s100530050041
    [2] Hartmann, "Casimir effect in the plane-sphere geometry: Beyond the
        proximity force approximation", phd thesis (2018),
        https://opus.bibliothek.uni-augsburg.de/opus4/44798
    """
    def __init__(self, omega, eps2, omegap=9, gamma=0.035):
        """Initialize the Material class

        Parameters:
            omega:  numpy array with (real) frequencies in rad/s
            eps2:   numpy array with imaginary parts of the dielectric function
                    corresponding to the array omega
            omegap: plasma frequency used for extrapolation for small
                    frequencies; in eV/hbar
            gamma:  relaxation frequency used for extrapolation for small
                    frequencies; in eV/hbar
        """
        self.omega  = omega
        self.eps2   = eps2
        self.omegap = omegap
        self.gamma  = gamma
        self.omegac = omega[-1]

    def epsm1(self, xi):
        """Compute epsilon(i*xi)-1"""
        xi  = np.array(xi)
        eps = np.empty(xi.shape)

        gamma = self.gamma
        omegap = self.omegap
        omegac = self.omegac

        for i,xi_i in enumerate(xi):
            # extrapolation using Drude model, (2.55) of [2]
            I1 = omegap**2/(xi_i*(xi_i**2-gamma**2))*( xi_i*atan(omegac/gamma) - gamma*atan(omegac/xi_i) )
            
            # (2.53) of [2]
            integrand = self.omega*self.eps2/(self.omega**2+xi_i**2)
            I2 = -simps(integrand, self.omega)

            # (2.53) of [2]
            eps[i] = 2/pi*(I1+I2)

        return eps


if __name__ == "__main__":
    from argparse import ArgumentParser
    from sys import argv

    command = " ".join(argv)

    parser = ArgumentParser(description="Rotate optical data to imaginary axis")

    # required argument
    parser.add_argument("filename", type=str, help="input filename")

    # optional arguments
    parser.add_argument("--start",  action="store", dest="start", type=float, default=1e11, help="omega start (in rad/s)")
    parser.add_argument("--stop",   action="store", dest="stop",  type=float, default=1e19, help="omega stop (in rad/s)")
    parser.add_argument("--points", action="store", dest="N",     type=int,   default=1024, help="number of points")

    parser.add_argument("--omegap_low",  action="store", dest="omegap_low", type=float, default=9,     help="omegap low (in eV)")
    parser.add_argument("--gamma_low",   action="store", dest="gamma_low",  type=float, default=0.035, help="gamma low (in eV)")

    parser.add_argument("--omegap_high", action="store", dest="omegap_high", type=float, default=0, help="omegap high (in eV)")
    parser.add_argument("--gamma_high",  action="store", dest="gamma_high",  type=float, default=0, help="gamma high (in eV)")

    parser.add_argument("--description", action="store", dest="desc", type=str, default="unknown", help="short description of the material")

    args = parser.parse_args()

    data = np.loadtxt(args.filename, delimiter=",", usecols=(0,3))
    material = Material(data[:,0], data[:,1], omegap=args.omegap_low, gamma=args.gamma_low)

    hbar_eV = 6.582119514e-16 # eV s

    omega = np.logspace(log10(args.start), log10(args.stop), args.N)
    eps2  = 1+material.epsm1(omega*hbar_eV)

    print("# %s" % command)
    print("#")
    print("# This file contains the dielectric function epsilon(i*xi) for %s" % args.desc)
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
    print("# xi in rad/s\tepsilon(i*xi)")

    for i,omega_ in enumerate(omega):
        print("%.8e\t%.8e" % (omega_, eps2[i]))

import numpy as np
from glob import glob
from math import pi
from pfa import pressure,force

# constants
kB   = 1.38064852e-23 # m² kg 1/s² 1/K
hbar = 1.0545718e-34  # m² kg / s
hbar_eV = 6.582119514e-16 # hbar [eV s]
c    = 299792458      # m/s

def slurp(filenames):
    data = []

    for i,filename in enumerate(filenames):
        with open(filename, "r") as f:
            drude = plasma = 0
            for line in f:
                line = line.strip()
                if "# plasma" in line:
                    _,line = line.split("=", 1)

                    line = line.strip()
                    plasma = float(line[:line.find(" ")])
                    continue
                if "# xi=0" in line:
                    line = line[16:]
                    line = line[:line.find(",")]
                    drude = float(line)
                    continue
                if line == "" or line[0] == "#":
                    continue

                # L/R, L, R, T, ldim, F*(L+R)/(ħc)
                LbyR, L, R, T, ldim, F_drude = map(float, line.split(","))
                T_scaled = 2*pi*kB*(L+R)/(hbar*c)*T
                F_plasma = F_drude + (plasma-drude)/2*T_scaled/pi

                data.append((L,R,T,ldim,F_drude,F_plasma))

    return np.array(sorted(data))


def deriv(x,y, deriv=1, accuracy=2):
    coefficients = {
        (1,2): (-1/2, 0, 1/2),
        (1,4): (1/12, -2/3, 0, 2/3, -1/12),
        (1,6): (-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60),
        (1,8): (1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280),

        (2,2): (1, -2, 1),
        (2,4): (-1/12, 4/3, -5/2, 4/3, -1/12),
        (2,6): (1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90),
        (2,8): (-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560)
    }

    npts = len(x)
    h = (x[-1]-x[0])/(npts-1)

    c = coefficients[(deriv, accuracy)]
    p = len(c)//2

    dy = np.zeros(npts-2*p)
    for i in range(npts-2*p):
        for j in range(len(c)):
            dy[i] += c[j]*y[i+j]

    dy /= h**deriv
    return np.copy(x[p:-p]), dy



if __name__ == "__main__":
    from sys import argv

    accuracy = 6

    omegap = 9/hbar_eV    # plasma frequency 9eV
    gamma  = 0.03/hbar_eV # dissipation 0.03
    filename = "../../materials/GoldEpsIm.dat"

    data = slurp(argv[1:])
    R = data[0,1]
    T = data[0,2]

    L = data[:,0]
    E_drude  = data[:,4]/(L+R)*(hbar*c)
    E_plasma = data[:,5]/(L+R)*(hbar*c)

    dx, F_drude  = deriv(L, E_drude,  accuracy=accuracy)
    dx, F_plasma = deriv(L, E_plasma, accuracy=accuracy)

    F_drude *= -1
    F_plasma *= -1

    print("# L (m), R (m), T (K), F_Drude, F_Plasma, F_PFA_Drude, F_PFA_Plasma")
    for i,L in enumerate(dx):
        pfa_drude, pfa_plasma = force(L,R,T)
        ratio_drude  = F_drude[i]/pfa_drude
        ratio_plasma = F_plasma[i]/pfa_plasma

        print("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %g, %g" % (L, R, T, F_drude[i], F_plasma[i], pfa_drude, pfa_plasma, ratio_drude, ratio_plasma))

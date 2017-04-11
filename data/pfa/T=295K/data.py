import numpy as np
from glob import glob
from math import pi

# constants
kB   = 1.38064852e-23 # m² kg 1/s² 1/K
hbar = 1.0545718e-34  # m² kg / s
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
    try:
        accuracy = int(argv[1])
    except:
        pass

    R = 151.3e-6

    d = {}
    pfa = np.loadtxt("data_pfa", delimiter=",")
    for i,L in enumerate(pfa[:,0]):
        d[L] = (pfa[i,1], pfa[i,2])

    data = slurp(sorted(glob("gold_eta7/gold*.out"), reverse=True))
    L = data[:,0]
    E_drude  = data[:,4]/(L+R)*(hbar*c)
    E_plasma = data[:,5]/(L+R)*(hbar*c)

    dx, dF_drude  = deriv(L, E_drude,  deriv=2, accuracy=accuracy)
    dx, dF_plasma = deriv(L, E_plasma, deriv=2, accuracy=accuracy)

    dx *= 1e9
    P_drude  = 1/(2*pi*R)*dF_drude *1000
    P_plasma = 1/(2*pi*R)*dF_plasma*1000

    print("# accuracy=%d" % accuracy)
    print("# L (nm), P (Drude, mPa), P (Plasma, mPa), P (Drude, PFA), P (Plasma, PFA), ratio (Drude), ratio (Plasma)")
    for i,x in enumerate(dx):
        index = int(round(x,8))
        print("%.15g, %.15g, %.15g, %.15g, %.15g, %.15g, %.15g" % (x, P_drude[i], P_plasma[i], d[index][0], d[index][1], P_drude[i]/d[index][0], P_plasma[i]/d[index][1]))

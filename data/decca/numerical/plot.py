import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def get_radius(filename):
    with open(filename, "r") as f:
        for line in f:
            if "# R =" in line:
                _, R = line.split("=")
                R = R.strip()
                return float(R[:-2])
    raise ValueError("Could not find radius in file")


for filename in argv[1:]:
    R = get_radius(filename)

    data = np.loadtxt(filename, delimiter=",")

    L = data[:,0]
    plasma     = data[:,3]
    plasma_pfa = data[:,4]

    drude     = data[:,1]
    drude_pfa = data[:,2]

    plt.plot(L, (drude/drude_pfa-1)*R/L)
    plt.plot(L, (plasma/plasma_pfa-1)*R/L)

plt.show()

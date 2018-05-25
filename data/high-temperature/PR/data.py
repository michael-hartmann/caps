import numpy as np
from glob import glob

data = []
for filename in glob("raw/*.out"):
    LbyR, L, R, ldim, F_drude, F_PR = np.loadtxt(filename, delimiter=",")
    if LbyR < 3e-5:
        print(filename)
    data.append((1/LbyR, F_drude, F_PR))

print("# high-temperature limit for perfect reflectors")
print("#")
print("# R/L, F_Drude/(kB*T), F_PR/(kB*T)")

for RbyL, drude, pr in sorted(data):
    print("%.8g, %.14g, %.14g" % (RbyL, drude, pr))

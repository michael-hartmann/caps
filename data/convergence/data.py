import numpy as np
from glob import glob

data = {}
for fname in glob("raw/slurm-*.out"):
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) and line[0] != "#":
                LbyR, L, R, T, ldim, E = map(float, line.split(","))
                RbyL = int(round(1/LbyR,0))
                eta = ldim/RbyL

                data.setdefault((R,RbyL,int(T)),[]).append([ldim,eta,E])


for key in data:
    data[key] = sorted(data[key])
    Eref = data[key][-1][-1]
    del data[key][-1]
    for i,x in enumerate(data[key]):
        E = data[key][i][-1]
        data[key][i].append(1-E/Eref)


for key in data:
    R,RbyL,T = key

    fname = "conv_%d_%d.csv" % (RbyL,T)
    with open(fname, "w") as f:
        print("# R=%gµm, R/L=%g, T=%dK" % (R*1e6,RbyL,T), file=f)
        print("# eta, ldim, E, (E-E_exact)/E_exac", file=f)

        for ldim,eta,E,ratio in data[key]:
            print("%g, %d, %g, %g" % (eta,ldim,E,ratio), file=f)

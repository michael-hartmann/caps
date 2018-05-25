import numpy as np
from math import log10
from os import system, unlink
import random, string
from sys import argv

# script to measure and compare the run-time using HODLR and Cholesky

def runtime(LbyR, m, nT, detalg, repeat=3):
    time = 0

    for i in range(repeat):
        filename = "/tmp/" + ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(10)) + ".csv"
        cmd = "./casimir_logdetD -L %g -R 1 --nT %g -m %d -d %s > %s" % (LbyR, nT, m, detalg, filename)
        system(cmd)
        
        with open(filename) as f:
            for line in f:
                if line == "" or line[0] == "#":
                    continue
                time += float(line.split(",")[-1])
                v = float(line.split(",")[4])

        unlink(filename)

    return v,time/repeat


m = int(argv[1])
nT = float(argv[2])

print("# R/L, m, nT, CHOLESKY, HODLR")
for LbyR in np.logspace(log10(0.1), log10(0.0005), 25):
    vc, cholesky  = runtime(LbyR, m, nT, "CHOLESKY", repeat=1)
    vh, hodlr     = runtime(LbyR, m, nT, "HODLR",    repeat=1)

    print("%.8g, %d, %g, %g, %g, %e" % (1/LbyR, m, nT, cholesky, hodlr, 1-abs(vh/vc)))

#LbyR = float(argv[1])
#m = int(argv[2])
#
#for xiRbyc in np.logspace(log10(0.01), log10(1000), 200):
#    nT = xiRbyc*(1+LbyR)
#
#    vmax = maximum(LbyR, m, nT, "LU")
#    print(LbyR, m, nT, vmax)

import numpy as np
from glob import glob

from sys import path, argv, stdout, stderr
path.append("../../../src/python")
import PFA
from deriv import deriv

hbarc = PFA.hbarc

def write(s, file):
    print(s)
    print(s, file=file)

def slurp(filenames):
    d = []
    for filename in filenames:
        with open(filename) as f:
            for line in f:
                if line and line[0] != "#":
                    LbyR, L, R, T, ldim, E = map(float, line.split(","))

                    d.append((L,R,T,E))

    return np.array(sorted(d))


def usage(f=stdout):
    print("data.py DIRECTORIES", file=f)
    print("  Generate csv with comparison to PFA", file=f)
    print("  This script generates a E.csv and a F.csv file in", file=f)
    print("  each directory from the raw data.", file=f)


if len(argv) < 2:
    usage(stderr)
    exit(1)

directories = argv[1:]
for directory in directories:
    data = slurp(glob(directory + "/slurm-*.out"))
    R = data[0,1]
    T = data[0,2]
    L = data[:,0]
    E = data[:,3]*(hbarc/(L+R))
    L_F, F = deriv(L,E)
    F *= -1

    epsm1 = PFA.epsilonm1_from_file("../../materials/GoldDalvit.dat")
    pfa = PFA.PFA(R,T,epsm1)

    with open(directory + "/E.csv", "w") as f:
        write("# R=%g" % R, file=f)
        write("# T=%g" % T, file=f)
        write("# R/L, L (in m), E (in J), E_pfa (in J), 1-E/E_PFA", file=f)
        for L_,E_ in zip(L,E):
            E_pfa = pfa.E(L_)
            E_ratio = E_/E_pfa
            write("%.12g, %.12g, %.12g, %.12g, %.12g" % (R/L_, L_, E_, E_pfa, 1-E_ratio), file=f)


    with open(directory + "/F.csv", "w") as f:
        write("# R=%g" % R, file=f)
        write("# T=%g" % T, file=f)
        write("# R/L, L (in m), F (in N), F_pfa (in N), 1-F/F_PFA", file=f)
        for L_,F_ in zip(L_F,F):
            F_pfa = pfa.F(L_)
            F_ratio = F_/F_pfa
            write("%.12g, %.12g, %.12g, %.12g, %.12g" % (R/L_, L_, F_, F_pfa, 1-F_ratio), file=f)

#!/usr/bin/python3

import numpy as np
from glob import glob
from math import pi,log,exp

def pfa(LbyR):
    # PFA formula
    return -pi**3/720.*(LbyR+1)/LbyR**2

def slurp(filenames):
    data = []
    for filename in filenames:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                empty = line == ""
                comment = line.startswith("#")
                if not(empty or comment):
                    # support old and new format
                    try:
                        # L/R, lmax, order, alpha, F(T=0)
                        LbyR,lmax,order,alpha,F = map(float, line.split(","))
                    except ValueError:
                        # L/R, L, R, T, ldim, F*(L+R)/(ħc)
                        LbyR,L,R,T,ldim,F = map(float, line.split(","))

                    data.append((LbyR, F, pfa(LbyR), F/pfa(LbyR)))

    return np.array(sorted(data, key=lambda x: x[0]))


if __name__ == "__main__":

    #bimonte = lambda x: 1+theta1*x+theta2*x**2*log(x)

    data = slurp(glob("pc_gk/slurm-*.out"))

    print("# L/R, F*(L+R)/(hbar*c), F_PFA*(L+R)/(hbar*c), F/F_PFA")
    for x in data:
        print("%.15g, %.15g, %.15g, %.15g" % tuple(x))

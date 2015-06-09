#!/ur/bin/python

from __future__ import division
from math import *
from numpy import logspace
from interpolate import Interpolate

kb    = 1.3806488e-23

listR = [ 0.1e-6, 0.2e-6, 0.5e-6, 1e-6, 2e-6]
T_SI   = 300
Lstart = 0.2e-6
Lstop  = 11e-6
N      = 50

inter = Interpolate("data2")

for R in listR:
    for L in logspace(log(Lstart),log(Lstop),N,base=e):
        LbyR, T_SI, S_SI = inter.interpolate_SI(L/R, R, T_SI)
        if LbyR >= 5:
            continue

        print R,L*1e6, S_SI*1e-15/(kb*R**3)
    print

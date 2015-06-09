#!/usr/bin/python

from __future__ import division
from math import *
import numpy as np
import pfa

N = 50
LbyR_start = 0.01
LbyR_stop  = 0.1

l = list(np.logspace(log(LbyR_start), log(LbyR_stop), N, base=e))
l.reverse()

print "# LbyR, a,b (by PFA)"
for i,LbyR in enumerate(l):
    x = []
    y = []
    for T in np.linspace(0.05, 0.2, 20):
        F = pfa.pfa(LbyR,T)
        x.append(T**4)
        y.append(F)
        print "# ", LbyR,T,y[-1]
    b,a = np.polyfit(x, y, 1)
    print "%.15g, %.10g, %.10g" % (LbyR,a,-b)

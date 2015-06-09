#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from numpy import polyfit
from glob import glob

directories = ["0.1_lo", "0.1_hi", "0.07_lo", "0.07_hi"]
for directory in directories:
    fit_x = []
    fit_y = []
    for filename in glob(directory + "/slurm-*.out"):
        fh = open(filename, "r")
        for line in fh:
            line = line.strip()
            if line == "" or line[0] == "#":
                continue
            Q,T,F,lmax,nmax,time = map(float, line.split(","))
            LbyR = 1/Q-1
            fit_x.append(T**4)
            fit_y.append(F)
        fh.close()

    b,a = polyfit(fit_x, fit_y, 1)
    print "directory=%s,\tlmax=%d,\tLbyR=%.10g,\ta=%.10g,\tb=%.10e" % (directory, lmax, LbyR, a,b)

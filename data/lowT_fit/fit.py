#!/usr/bin/python

from __future__ import division
from math import ceil
from numpy import linspace

precision  = 1e-8
lscale     = 7
lmin       = 15
T_min      = 0.05
T_max      = 0.2
N_T        = 20
LbyR_start = 4
LbyR_stop  = 0.1
N_LbyR     = 100

for i,LbyR in enumerate(linspace(LbyR_start, LbyR_stop, N_LbyR)):
    Q = 1/(1+LbyR)
    lmax = max(lmin, ceil(lscale/LbyR))

    cmd = "casimir -Q %.10f -T %f,%f,%d -p %g -L %d > raw/slurm-%05d.out" % (Q, T_min, T_max, N_T, precision, lmax,i)
    print cmd

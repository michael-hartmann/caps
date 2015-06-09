#!/usr/bin/python

from __future__ import division
import numpy as np

eps  = 1e-8
lfac = 6
lmin = 20

LbyR_start = 5
LbyR_stop  = 0.03
LbyR_N     = 1200

T_start    = 0.2
T_stop     = 1.6
T_N        = 450

for LbyR in np.linspace(LbyR_start, LbyR_stop, LbyR_N):
    Q = 1/(1+LbyR)
    lmax = max(lfac/LbyR, lmin)

    if lmax <= 60:
        slurmd = "sbatch ~/run"
        cmd = "casimir -T %f,%f,%d -Q %.15g -p %g -L %d # L/R=%g" % (T_start, T_stop, T_N, Q, eps, lmax, LbyR)

        print slurmd + " " + cmd
    else:
        cores = 5

        for T in np.linspace(T_start, T_stop, T_N):
            slurmd = "sbatch -c %d ~/run" % cores
            cmd = "casimir -T %.15g -Q %.15g -p %g -L %d -c %d # L/R=%g, multicore" % (T, Q, eps, lmax, cores, LbyR)

            print slurmd + " " + cmd

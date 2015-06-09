#!/usr/bin/python

from __future__ import division
import numpy as np

for i,LbyR in enumerate(np.linspace(5,20,200)):
    Q = 1/(1+LbyR)
    print "casimir -Q %.15g -T 0.05,0.25,20 -p 1e-10 -L 10 > slurm_addita_%05d.out" % (Q,i)

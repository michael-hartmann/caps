#!/usr/bin/python

from __future__ import division
import sys
from glob import glob
from scipy.interpolate import interp1d

F0 = []
fh = open("F0", "r")
for line in fh:
    line = line.strip()
    if line == "" or line[0] == "#":
        continue

    Q, LbyR, F, S = map(float, line.split(","))
    F0.append((LbyR, F, S))
fh.close()

F0.sort()
inter = interp1d([x[0] for x in F0], [x[1] for x in F0])

data = []
for filename in glob("density/slurm-*.out"):
    fh = open(filename, "r")
    for line in fh:
        line = line.strip()
        if len(line) == 0 or line[0] == '#' or "slurmd" in line:
            continue

        data.append(map(float, line.split(",")))
        Q = data[-1][0]
        LbyR = 1/Q-1
        try:
            F0 = inter(LbyR)
            #(Q0, LbyR0, F0, S0)
            F = data[-1][2]
            data[-1].append(F/F0)
        except:
            #print LbyR
            data[-1].append(-1)
        
    fh.close()

# function that sorts a list of lists by first and second item
def f_cmp(x,y):
    if x[0] == y[0]:
        return cmp(x[1], y[1])
    else:
        return cmp(x[0], y[0])
data.sort(cmp=f_cmp)

for i in range(len(data)):
    try:
        if data[i][0] == data[i+1][0] and data[i][1] == data[i+1][1]:
            del data[i+1]
    except:
        break

d = {}

for p in data:
    key = p[0]
    if not d.has_key(key):
        d[key] = []

    d[key].append(p)

keys = d.keys()
keys.sort()
keys.reverse()
print "# L/R, R/(L+R), T, F, S, F/F0, lmax, nmax, time"
for key in keys:
    print
    #for i in range(1,len(d[key])-1):
    for i in range(2,len(d[key])-2):
        p = d[key][i]
        #S = - (d[key][i+1][2] - d[key][i-1][2])/(d[key][i+1][1] - d[key][i-1][1])
        S = - (-d[key][i+2][2]+8*d[key][i+1][2]-8*d[key][i-1][2]+d[key][i-2][2])/(12*(d[key][i+1][1] - d[key][i-1][1])/2)
        LbyR = 1/p[0]-1
        RbyScriptL, T, F, lmax, nmax, time, FbyF0 = p
        print "%.10g, %.10g, %.10g, %.10g, %.10g, %.10g, %g, %g, %g" % (LbyR, RbyScriptL, T, F, S, FbyF0, lmax, nmax, time)

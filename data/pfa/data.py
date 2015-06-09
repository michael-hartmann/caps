#!/usr/bin/python

from math import *
from glob import glob
import numpy as np

def slurp(filename):
    """
    Datei filename einlesen und die Daten als dict zurueckgeben. Der key
    entspricht n, die dazugehoerigen values entsprechen einer sortierten Liste
    bestehend aus den Elementen (m,value).

    Return: L/R, T, dict
    """
    data = {}
    LbyR,T = None,None
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if line == "":
                continue
            if line[0] == "#" and "m=" in line:
                line = line[1:]
                line.strip
                n,m,value = line.split(",")
                n     = int(n[3:])
                m     = int(m[3:])
                value = float(value[7:])
                if value == 0:
                    continue

                if n not in data:
                    data[n] = []
                data[n].append((m,value))
            elif "#" not in line:
                LbyR,T,F,lmax,nmax,time = map(float, line.split(","))

    for key in data.keys():
        data[key] = sorted(data[key])

    return LbyR,T,data


def estimate_remainder(data,points=6):
    """
    Die Funktion nimmt die letzten points Punkte, fittet eine
    Exponentialfunktion durch diese Punkte und nutzt die Parameter, um den Rest
    der abgeschnittenen Summe abzuschaetzen.
    """
    k0 = len(data)-1
    x = [    +z[0]  for z in data[-points:]]
    y = [log(-z[1]) for z in data[-points:]]

    m,t = np.polyfit(x,y,1)
    
    a0,q = exp(t),exp(m)
    return -a0*q**k0/(1-q)


def read_file(filename, points=8):
    """
    Die Datei filename einlesen, fuer die letzten 8 Punkte die abgeschnittene Summer fuer F_n und F abschaetzen.

    Return: L/R, T, F, error
    """
    LbyR,T,data = slurp(filename)
    result = []

    for key in sorted(data.keys()):
        remainder_m = estimate_remainder(data[key])

        Fn = sum([z[1] for z in data[key]]) - data[key][0][1]/2

        result.append((key, Fn+remainder_m, Fn, remainder_m))

    remainder_n = estimate_remainder(result)
    F_T = T/pi*(sum([z[1] for z in result]) - result[0][1]/2 + remainder_n)

    error = T/pi*( remainder_n + result[0][3]/2 + sum([z[3] for z in result[1:]]) )

    return LbyR,T,F_T,error


if __name__ == "__main__":
    data = []
    for filename in glob("slurm-*.out"):
        data.append(read_file(filename))
    data = sorted(data)

    print "# LbyR, T, F, error"
    for entry in data:
        print "%.15g, %.15g, %.15g, %g" % (entry)

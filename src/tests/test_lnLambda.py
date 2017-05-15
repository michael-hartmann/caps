#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 100

def prettyprint(x, length=20):
    return str(x)[:length]


def lnLambda(l1,l2,m):
    return log(sqrt( (2*l1+1)*(2*l2+1)*factorial(l1-m)*factorial(l2-m) /( factorial(l1+m)*factorial(l2+m)*l1*(l1+1)*l2*(l2+1) ) ))


if __name__ == "__main__":
    l = []

    for l1 in (1,2,3,5,10,50,100,200,500,1000,2000,5000,7000,10000,20000,50000,75000,100000):
        for l2 in (1,2,3,5,10,50,100,200,500,1000,2000,5000,7000,10000,20000,50000,75000,100000):
            if l1 >= l2:
                for m in (0,1,2,5,10,20,50,100,20,500,1000,2000,5000):
                    if l2 >= m:
                        l.append((l1,l2,m))


    for l1,l2,m in l:
        s = prettyprint(lnLambda(l1,l2,m))
        print("AssertAlmostEqual(&test, casimir_lnLambda(%d,%d,%d), %s);" % (l1,l2,m,s))

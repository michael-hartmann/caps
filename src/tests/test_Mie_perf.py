#!/usr/bin/python

from __future__ import division
from mpmath import *
import numpy as np

mp.dps = 100

def prettyprint(args, length=40):
    lna,sign_a,lnb,sign_b = args
    s_lna    = str(lna)[:length]
    s_lnb    = str(lnb)[:length]
    d_sign_a = int(sign_a)
    d_sign_b = int(sign_b)

    return "%s, %+d, %s, %+d" % (s_lna, d_sign_a, s_lnb, d_sign_b)


def lnab(l,chi):
    al = (-1)**(l+1) * pi/2 * ( l*besseli(l+0.5,chi)-chi*besseli(l-0.5,chi) )/( l*besselk(l+0.5,chi)+chi*besselk(l-0.5,chi) )
    bl = (-1)**(l+1) * pi/2 * besseli(l+0.5,chi)/besselk(l+0.5,chi)

    return log(abs(al)), sign(al), log(abs(bl)), sign(bl)


if __name__ == "__main__":
    for l in (100,200,500,1000,1500,2000,4000,6000,8000):
        print l,prettyprint(lnab(l,1))

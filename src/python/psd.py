#!/usr/bin/python

# Ref: Hu, Xu, Yan, Pade spectrum decomposition of Fermi function and Bose
# function, The Journal of Chemical Physics 133, 101106, 2010

from __future__ import division
import numpy as np
from sympy import *
from mpmath import *

mp.dps = 200

N = 4

z = symbols("z")

Amm,Am = 1/4, 5/4
Bmm,Bm = 3, 15+z/4

for M in range(3,2*N+1):
    A = simplify((2*M+1)*Am + z*Amm/4)
    B = simplify((2*M+1)*Bm + z*Bmm/4)

    Bmm,Bm = Bm,B
    Amm,Am = Am,A


A = Poly(A)
B = Poly(B)
Bp = B.diff(z)
coeffs = map(float, reversed(B.coeffs()))

roots = np.polynomial.polynomial.polyroots(coeffs)
roots = sorted(roots, key=lambda x: -x)

print "double psd%d[][2] = {" % N
for i,mxi2 in enumerate(roots):
    xi_i  = sqrt(-mxi2)/(2*(i+1)*pi)
    eta_i = A.eval(mxi2)/Bp.eval(mxi2)/2
    print "    { %.20g, %.20g },"  % (xi_i, eta_i)
print "    { 0, 0 }"
print "};"

# Test the implementation of the dilogarithm (or spence function) for
# parameters 0<x<=1 against the implementation of mpmath.
#
# The script outputs the maximum relative error found.

import numpy as np
import mpmath as mp
from spence import Li2

mp.dps = 50

def Li2_mpmath(x):
    return mp.polylog(2,mp.mpf(x))

maxerr = 0
x_maxerr = 0
for x in np.linspace(0,1,25000):
    if x == 0:
        continue

    exact = Li2_mpmath(x)
    fast  = Li2(x)
    err   = abs(1-fast/exact) # relative error
    
    if err >= maxerr:
        maxerr = err
        x_maxerr = x

print("Maximum relative error %g at x=%g" % (maxerr, x_maxerr))

from math import ceil
import mpmath
from mpmath import mpf,log,gamma,sqrt

mpmath.mp.dps = 50

def logdetD(LbyR, m, lmax):
    """Calculate logdetD for LbyR = L/R, m and lmax for Matsubara frequency
    xi=0
    """
    y = 1/(1+LbyR)/2

    minimum = max(1,m)
    maximum = lmax
    dimension = maximum-minimum+1

    EE = mpmath.matrix(dimension, dimension)
    MM = mpmath.matrix(dimension, dimension)
    for l1 in range(minimum,maximum+1):
        for l2 in range(l1,maximum+1):
            i,j = l1-minimum,l2-minimum

            # Kronecker delta
            delta = 1 if l1==l2 else 0

            v = y**(1+l1+l2)*gamma(1+l1+l2)/sqrt( gamma(1+l1+m)*gamma(1+l1-m) * gamma(1+l2+m)*gamma(1+l2-m) )

            EE[i,j] = EE[j,i] = delta-v
            MM[i,j] = MM[j,i] = delta-v*sqrt(l1*l2/(l1+1)/(l2+1))

    detEE = mpmath.det(EE)
    detMM = mpmath.det(MM)
    return float(log(detEE)), float(log(detMM))


if __name__ == "__main__":
    for LbyR in (1, 0.5, 0.1, 0.05, 0.01):
        for m in (0,1,2,3,4,5,10,20):
            lmax = int(ceil(max(30, 7/LbyR)))
            logdetEE,logdetMM = logdetD(mpf(LbyR),m,lmax)
            print("    { %.14g, %d, %d, %.14g, %.14g }," % (LbyR, m, lmax, logdetEE, logdetMM))
        print()

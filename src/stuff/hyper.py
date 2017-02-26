from math import *

def lfac(n):
    return lgamma(1+n)

def hyper(a,b,c,z):
    """Evaluate hypergeometric function using the power series"""
    cs = 1/gamma(c)
    v = cs

    for s in range(1,521):
        cs *= (a-1+s)*(b-1+s)/((c+s-1)*s)*z
        v += cs
        print(s,cs,v)
        if cs/v < 1e-16:
            return v

    return v


def Plm(l,m,x):
    prefactor = lfac(l+m)-m*log(2)-lfac(l-m) + m/2*log(x**2-1)
    return prefactor + log(abs(hyper(l+m+1, m-l, m+1, (1-x)/2)))



print(exp(Plm(2500,2,1.01)))

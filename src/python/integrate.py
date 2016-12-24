import numpy as np
from scipy.special import factorial2
import mpmath
from mpmath import exp,quad,inf,log,sqrt,mpf

def Plm(l,m,z):
    array = [0] * (l-m+1)
    if m == 0:
        array[0] = 1;
    else:
        array[0] = factorial2(2*m-1,exact=True)*((z+1)*(z-1))**(m/2)

    if l == m:
        return array[-1];

    array[1] = z*(2*m+1)*array[0]

    for ll in range(m+2, l+1):
        array[ll-m] = ((2*ll-1)*z*array[ll-m-1] - (ll+m-1)*array[ll-m-2])/(ll-m)

    return array[-1]

def Plm_mpmath(l,m,z):
    return mpmath.legenp(l,m,z,type=3)


def integrand(nu,m, tau):
    #return lambda z: exp(-tau*z)/(z**2+2*z)*Plm(nu,2*m,1+z)
    return lambda z: exp(-tau*z)*Plm(nu,2,1+z)


if __name__ == "__main__":
    from sys import argv
    nu  = int(argv[1])
    m   = int(argv[2])
    tau = float(argv[3])

    v = quad(integrand(nu,m,tau), [0,inf])
    print(v)
    print(nu,m,tau,log(v))

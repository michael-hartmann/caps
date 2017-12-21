cimport libc.math
cimport cython

@cython.cdivision(True)
def Li2(double x):
    """Implement the polylogarithm Li_s(x) for s=2 (also known as dilogarithm
    or Spencer function) for 0 <= x <= 1.

    This function uses a series with n^4 convergence, see [1]. In order to
    improve convergence, we make use of the identity:
        Li2(x)+Li2(1-x) = pi^2/6 - log(x)*log(1-x)

    References:
    [1] Robert Morris, The dilogarithm function of a real argument, Math. Comp. 33 (1979)
    """
    cdef double s = 0
    cdef int n
    cdef double xn = 1
    cdef double pi2 = 9.869604401089358 # pi²

    if x == 0:
        return 0
    elif x == 1:
        return pi2/6
    elif x < 0.5:
        for n in range(1,35):
            xn = xn*x
            s += xn/(n*n*(n+1)*(n+1))

        return x/(x+1)*(3+s)-2*(x-1)/(x+1)*libc.math.log1p(-x)
    else:
        return pi2/6-libc.math.log(x)*libc.math.log1p(-x)-Li2(1-x)

cimport libcasimir
from libc.math cimport exp

class constants:
    hbar    = 1.0545718e-34
    hbar_eV = 6.582119514e-16
    kB      = 1.38064852e-23
    c       = 299792458


class sfunc:
    def __check_parameters_plm(int l, int m, double x):
        if l <= 0:
            raise ValueError("l must be positive")
        if m < 0:
            raise ValueError("m must be non-negative")
        if m > l:
            raise ValueError("m must not be large than l")
        if x <= 1:
            raise ValueError("x must be larger than 1")

    def lnPlm(int l, int m, double x):
        """Compute associated Legendre function Plm for x>1."""
        sfunc.__check_parameters_plm(l,m,x)
        return libcasimir.lnPlm(l,m,x)

    def Plm(int l, int m, double x):
        """Compute associated Legendre function Plm for x>1."""
        return exp(lnPlm(l,m,x))

    def Plm_continued_fraction(long l, long m, double x):
        """Calculate fraction P_l^{m-1}(x)/P_l^m(x)"""
        sfunc.__check_parameters_plm(l,m,x)
        return libcasimir.Plm_continued_fraction(l, m, x)

    def dlnPlm(int l, int m, double x):
        cdef double dplm
        cdef double d2plm
        sfunc.__check_parameters_plm(l,m,x)
        dplm = dlnPlm(l, m, x, &d2plm)
        return dplm,d2plm

    def bessel_In(int n, double x):
        return libcasimir.bessel_In(n, x)

    def bessel_Kn(int n, double x):
        return libcasimir.bessel_Kn(n, x)

    def bessel_logIn(int n, double x):
        return libcasimir.bessel_logIn(n, x)

    def bessel_logKn(int n, double x):
        return libcasimir.bessel_logKn(n, x)

    def bessel_ratioI(double nu, double x):
        return libcasimir.bessel_ratioI(nu, x)

    def bessel_logInu_half(int nu, double x):
        return libcasimir.bessel_logInu_half(nu, x)

    def bessel_logKnu_half(int nu, double x):
        return libcasimir.bessel_logKnu_half(nu, x)

    def lfac(unsigned int n):
        return libcasimir.lfac(n)

    def logi(unsigned int n):
        return libcasimir.logi(n)


class material:
    pass

class integration:
    pass

class casimir:
    pass

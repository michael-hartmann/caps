cimport libcasimir

class constants:
    hbar    = 1.0545718e-34
    hbar_eV = 6.582119514e-16
    kB      = 1.38064852e-23
    c       = 299792458


class sfunc:
    def __check_parameters_plm(l, m, x):
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

    def besselI(int n, double x):
        return libcasimir.besselI(n, x)

    def besselI0e(double x):
        return libcasimir.besselI0e(x)

    def besselI1e(double x):
        return libcasimir.besselI1e(x)

    def bessel_lnInu(int nu, double x):
        return libcasimir.bessel_lnInu(nu, x)

    def bessel_lnKnu(int nu, double x):
        return libcasimir.bessel_lnKnu(nu, x)

    def bessel_lnInuKnu(int nu, double x):
        cdef double Inu, Knu
        libcasimir.bessel_lnInuKnu(nu, x, &Inu, &Knu)
        return Inu, Knu

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

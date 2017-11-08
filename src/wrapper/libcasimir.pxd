cdef extern from "plm.h":
    double lnPlm(int l, int m, double x)
    double lnPlm_estimate(int l, int m, double x)
    double Plm_continued_fraction(const long l, const long m, const double x)

cdef extern from "bessel.h":
    double besselI0e(double x)
    double besselI1e(double x)
    double besselI(int n, double x)
    double bessel_continued_fraction(int nu, double x)
    double bessel_lnInu(int nu, double x)
    double bessel_lnKnu(int nu, double x)
    void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p)

cdef extern from "logfac.h":
    double lfac(unsigned int n)
    double logi(unsigned int x)

cdef extern from "plm.h":
    double lnPlm(int l, int m, double x)
    double Plm_continued_fraction(const long l, const long m, const double x)
    double dlnPlm(int l, int m, double x, double *d2lnPlm);

cdef extern from "bessel.h":
    double bessel_In(int n, double x);
    double bessel_Kn(int n, double x);
    double bessel_logIn(int n, double x);
    double bessel_logKn(int n, double x);
    double bessel_ratioI(double nu, double x);
    double bessel_logInu_half(int nu, double x);
    double bessel_logKnu_half(int nu, double x);


cdef extern from "logfac.h":
    double lfac(unsigned int n)
    double logi(unsigned int x)

import cython
from longdouble cimport longdouble_c2py, longdouble_py2c
import math
from libcpp cimport bool

ctypedef long double float80
ctypedef signed char sign_t

cdef extern from "libcasimir.h":
    ctypedef struct casimir_t:
        double LbyR
        int lmax

    int casimir_init(casimir_t *self, double LbyR, double T)
    void casimir_free(casimir_t *self)

    int casimir_set_lmax(casimir_t *self, int lmax)
    int casimir_set_precision(casimir_t *self, double precision)
    void casimir_set_debug(casimir_t *self, bool debug)
    int casimir_set_cores(casimir_t *self, int cores)
    int casimir_set_trace_threshold(casimir_t *self, double threshold)
    void casimir_set_verbose(casimir_t *self, bool verbose)

    float80 casimir_lnLambda(int l1, int l2, int m, sign_t *sign)
    float80 casimir_lnXi(int l1, int l2, int m, sign_t *sign)
    double casimir_epsilon(double xi, double omegap, double gamma_)

    void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0)
    void casimir_lnab(casimir_t *self, const int n_mat, const int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)

    double casimir_F(casimir_t *self, int *nmax) nogil
    double casimir_F_n(casimir_t *self, const int n, int *mmax)
    void casimir_logdetD0(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
    double casimir_logdetD(casimir_t *self, int n, int m)
    double casimir_lnepsilon(double xi, double omegap, double gamma_)


cdef extern from "sfunc.h":
    float80 logadd(const float80 a, const float80 b);
    float80 logadd_s(const float80 a, const sign_t sign_a, const float80 b, const sign_t sign_b, sign_t *sign);

    void bessel_lnInuKnu(int nu, const float80 x, float80 *lnInu_p, float80 *lnKnu_p);
    float80 ln_doublefact(int n);

    void bessel_lnInuKnu(int nu, const float80 x, float80 *lnInu_p, float80 *lnKnu_p);

    float80 gaunt_log_a0(int n, int nu, int m);
    float80 gaunt_a0(int n,int nu,int m);
    void gaunt(const int n, const int nu, const int m, float80 a_tilde[]);
    int gaunt_qmax(const int n, const int nu, const int m);



class sfunc:
    def ln_doublefact(int n):
        """Calculate double factorial n!!."""
        return longdouble_c2py( ln_doublefact(n) )


    def logadd(a, b):
        """Return log(exp(a)+exp(b))."""
        return longdouble_c2py( logadd(longdouble_py2c(a), longdouble_py2c(b)) )


    def logadd_s(a, sign_a, b, sign_b):
        """Calculate v = sign_a*exp(a)+sign_b*exp(b) and return log|v| and sign(v)."""
        cdef result
        cdef sign_t sign
        cdef sign_t sign_a_c = math.copysign(1, sign_a)
        cdef sign_t sign_b_c = math.copysign(1, sign_b)

        result = logadd_s(longdouble_py2c(a), sign_a, longdouble_py2c(b), sign_b_c, &sign)

        return longdouble_c2py(result), <int>sign


    def bessel_lnInuKnu(int nu, x):
        cdef float80 Inu, Knu
        bessel_lnInuKnu(nu, longdouble_py2c(x), &Inu, &Knu)

        return longdouble_c2py(Inu), longdouble_c2py(Knu)


    def bessel_lnInu(int nu, x):
        return sfunc.bessel_lnInuKnu(nu,x)[0]


    def bessel_lnKnu(int nu, x):
        return sfunc.bessel_lnInuKnu(nu,x)[1]



cdef class Casimir:
    cdef casimir_t casimir

    def __init__(self, LbyR=1, T=1, verbose=False, debug=False, omegap_sphere=None, gamma_sphere=None, omegap_plane=None, gamma_plane=None, lmax=None, cores=None, precision=None, trace_threshold=None):
        if LbyR <= 0:
            raise ValueError("LbyR must be positive")
        if T <= 0:
            raise ValueError("T must be positive")
        if lmax != None and lmax <= 0:
            raise ValueError("lmax must be positive")
        if precision != None and precision <= 0:
            raise ValueError("precision must be positive")
        if cores != None and (type(cores) != int or cores < 1):
            raise ValueError("cores must be a positive integer")
        if trace_threshold != None and trace_threshold < 0:
            raise ValueError("trace_threshold must be non negative")

        assert casimir_init(&(self.casimir), LbyR, T) == 0

        if lmax:
            casimir_set_lmax(&self.casimir, lmax)
        if precision:
            casimir_set_precision(&self.casimir, precision)
        if debug:
            casimir_set_debug(&self.casimir, True)
        if cores:
            casimir_set_cores(&self.casimir, cores)
        if trace_threshold:
            casimir_set_trace_threshold(&self.casimir, trace_threshold)
        if verbose:
            casimir_set_verbose(&self.casimir, verbose)



    def lnLambda(Casimir self, int l1, int l2, int m):
        cdef sign_t sign
        cdef float80 res = casimir_lnLambda(l1, l2, m, &sign)
        return longdouble_c2py(res), <int>sign


    def lnXi(Casimir self, int l1, int l2, int m):
        cdef sign_t sign
        cdef float80 res = casimir_lnXi(l1, l2, m, &sign)
        return longdouble_c2py(res), <int>sign


    def epsilon(Casimir self, double xi, double omegap, double gamma_):
        return casimir_epsilon(xi, omegap, gamma_)


    def lnepsilon(Casimir self, double xi, double omegap, double gamma_):
        return casimir_lnepsilon(xi, omegap, gamma_)


    def rp(Casimir self, nT, k):
        raise NotImplementedError()


    def F_SI_to_scaled(Casimir self, double F_SI, double ScriptL):
        raise NotImplementedError()


    def F_scaled_to_SI(Casimir self, double F_Sscaled, double ScriptL):
        raise NotImplementedError()


    def T_SI_to_scaled(Casimir self, double T_SI, double ScriptL):
        raise NotImplementedError()

    
    def T_scaled_to_SI(Casimir self, double T_scaled, double ScriptL):
        raise NotImplementedError()

    
    def lnab0(Casimir self, int l):
        cdef sign_t sign_a0, sign_b0
        cdef double a0,b0
        casimir_lnab0(l, &a0, &sign_a0, &b0, &sign_b0)
        return a0,sign_a0,b0,sign_b0


    def lnab(Casimir self, int n_mat, int l):
        cdef double lna, lnb
        cdef sign_t sign_a, sign_b
        casimir_lnab(&self.casimir, n_mat, l, &lna, &lnb, &sign_a, &sign_b)
        return lna, sign_a, lnb, sign_b


    def F_n(Casimir self, int n):
        cdef int mmax
        cdef double result

        result = casimir_F_n(&self.casimir, n, &mmax)

        return result, mmax


    def F(Casimir self):
        cdef int nmax
        cdef double result
        
        with nogil:
            result = casimir_F(&(self.casimir), &nmax)

        return result, nmax


    def logdetD0(Casimir self, int m):
        cdef double logdet_EE, logdet_MM
        casimir_logdetD0(&self.casimir, m, &logdet_EE, &logdet_MM)
        return logdet_EE, logdet_MM


    def logdetD(Casimir self, int n, int m):
        return casimir_logdetD(&self.casimir, n, m)


    def __dealloc__(self):
        casimir_free(&self.casimir)

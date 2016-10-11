import numpy as np
import cython
import math
from libcpp cimport bool

ctypedef signed char sign_t

cdef extern from "libcasimir.h":
    ctypedef struct casimir_t:
        double LbyR
        double T
        double threshold
        int lmax

    int casimir_init(casimir_t *self, double LbyR, double T)
    void casimir_free(casimir_t *self)

    void casimir_set_debug(casimir_t *self, bool debug)
    bool casimir_get_debug(casimir_t *self)

    void casimir_set_verbose(casimir_t *self, bool verbose)
    bool casimir_get_verbose(casimir_t *self)

    int casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi, void *userdata), void *userdata)

    int casimir_get_lmax(casimir_t *self)
    int casimir_set_lmax(casimir_t *self, int lmax)

    int casimir_get_detalg(casimir_t *self, char detalg[128])
    int casimir_set_detalg(casimir_t *self, const char *detalg)

    int casimir_get_cores(casimir_t *self)
    int casimir_set_cores(casimir_t *self, int cores)

    double casimir_get_precision(casimir_t *self)
    int    casimir_set_precision(casimir_t *self, double precision)

    double casimir_lnLambda(int l1, int l2, int m)

    void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0)
    void casimir_lnab(casimir_t *self, int n, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)
    void casimir_lnab_perf(casimir_t *self, int n, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)

    void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM)

    void casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM)
    double casimir_logdetD(casimir_t *self, int n, int m)


cdef extern from "sfunc.h":
    ctypedef struct plm_combination_t:
        double lnPl1mPl2m;
        int sign_Pl1mPl2m;

        double lndPl1mPl2m;
        int sign_dPl1mPl2m;

        double lnPl1mdPl2m;
        int sign_Pl1mdPl2m;

        double lndPl1mdPl2m;
        int sign_dPl1mdPl2m;


    double ln_doublefact(int n);

    double bessel_lnInu(const int n, const double x);
    double bessel_lnKnu(const int n, const double x);
    void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p);

    double plm_lnPlm (int l, int m, double x, sign_t *sign);
    double plm_lndPlm(int l, int m, double x, sign_t *sign);

    void plm_PlmPlm(int l1, int l2, int m, double x, plm_combination_t *res);


class sfunc:
    """Special functions

    This module provides some special functions that are needed to evaluate the
    Casimir free energy in the plane-sphere geometry. Special functions
    implemented include modified Bessel functions for half-integer orders,
    associated Legendre polynomials, and the double factorial.
    """

    def ln_doublefact(int n):
        """Double factorial

        This is the factorial with every second value skipped, i.e.,
        7!! = 7*5*3*1.
        """
        return ln_doublefact(n)


    def lnInuKnu(int nu, double x):
        """Modified Bessel functions I_{nu+1/2}, K_{nu+1/2}

        Calculate modified Bessel functions I_{nu+1/2} and K_{nu+1/2} for
        the argument x, and return logarithm of value.

        The function returns log(I_{nu+1/2}(x)), log(K_{nu+1/2}(x)).

        Please note that the order is not nu, but nu+1/2.
        """
        cdef double Inu, Knu
        bessel_lnInuKnu(nu, x, &Inu, &Knu)
        return Inu, Knu

    def lnInu(int nu, double x):
        """Modified Bessel function I_{nu+1/2}(x) (see lnInuKnu)"""
        return bessel_lnInu(nu, x)

    def lnKnu(int nu, double x):
        """Modified Bessel function K_{nu+1/2}(x) (see lnInuKnu)"""
        return bessel_lnKnu(nu, x)


    def PlmPlm(int l1, int l2, int m, double x):
        """Product of associated Legendre polynomials

        Calculate the product of the associated Legendre polynomials
            (i)   Pl1m(x)*Pl2m(x)
            (ii)  Pl1m'(x)*Pl2m(x)
            (iii) Pl1m(x)*Pl2m'(x)
            (iv)  Pl1m'(x)*Pl2m'(x)
        for x > 1.

        The function returns the logarithm and the sign for the product (i) to
        (iv).
        """
        cdef plm_combination_t res
        plm_PlmPlm(l1, l2, m, x, &res)
        return res.lnPl1mPl2m, res.sign_Pl1mPl2m, res.lndPl1mPl2m, res.sign_dPl1mPl2m, res.lnPl1mdPl2m, res.sign_Pl1mdPl2m, res.lndPl1mdPl2m, res.sign_dPl1mdPl2m

    def lnPlm(int l, int m, double x):
        """Calculate log(Plm(x)) for x>1 and its sign"""
        cdef sign_t sign
        cdef double v
        v = plm_lnPlm(l, m, x, &sign)
        return v,sign

    def lndPlm(int l, int m, double x):
        """Calculate log(Plm'(x)) for x>1 and its sign"""
        cdef sign_t sign
        cdef double v
        v = plm_lndPlm(l, m, x, &sign)
        return v,sign


cdef class Casimir:
    cdef casimir_t casimir
    cpdef args
    cpdef epsilonm1

    def __init__(self, LbyR=1, T=1, verbose=False, debug=False, lmax=None, cores=None, precision=None, detalg=None):
        """Initialize Casimir object

        Required arguments:
            LbyR: aspect ration L/R, LbyR > 0
            T:    temperature in units of ħc/(2π*kB*(L+R)), T > 0

        Optional arguments:
            lmax:      truncation of vector space
            cores:     number of cores to use (only used for F)
            precision  XXX
            detalg:    LU, LU_LAPACK, QR_GIVENS or QR_LAPACK
            threshold: XXX
            debug:     flag, print debugging information
            verbose:   flag, print some addition information
        """
        cdef char c_detalg[128]

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
        if detalg != None and len(detalg) >= 128:
            raise ValueError("length of detalg must be smaller than 128")

        assert casimir_init(&(self.casimir), LbyR, T) == 0

        if lmax:
            casimir_set_lmax(&self.casimir, lmax)
        if precision:
            casimir_set_precision(&self.casimir, precision)
        if debug:
            casimir_set_debug(&self.casimir, True)
        if cores:
            casimir_set_cores(&self.casimir, cores)
        if verbose:
            casimir_set_verbose(&self.casimir, verbose)
        if detalg:
            # pad to 128 bytes
            detalg = detalg + (128-len(detalg))*"\0"
            detalg_ascii = bytes(detalg, "ascii")
            c_detalg = detalg_ascii
            casimir_set_detalg(&self.casimir, c_detalg)

    def __dealloc__(self):
        casimir_free(&self.casimir)



    def get_debug(self):
        return casimir_get_debug(&self.casimir)

    def get_verbose(self):
        return casimir_get_verbose(&self.casimir)

    def get_lmax(self):
        return casimir_get_lmax(&self.casimir)

    def get_detalg(self):
        cdef char detalg[128]
        casimir_get_detalg(&self.casimir, detalg)
        return str(<bytes>detalg, "ascii")

    def get_cores(self):
        return casimir_get_cores(&self.casimir)

    def get_precision(self):
        return casimir_get_precision(&self.casimir)

    def set_epsilonm1(Casimir self, epsilonm1, args):
        """Set callback function epsilonm1

        epsilonm1 must be of the form:

        def epsilonm1(xi, args):
            ...

        xi corresponds to the (scaled) Matsubara frequency, and the args given
        to set_epsilonm1 will be given to the function epsilonm1.
        """
        cdef void *ptr
        ptr = <void *>self
        self.epsilonm1 = epsilonm1
        self.args = args
        casimir_set_epsilonm1(&self.casimir, &__epsilonm1, ptr)

    def __repr__(self):
        s  = "L/R = %.8g\n" % self.casimir.LbyR
        s += "T   = %.8g\n" % self.casimir.T
        s += "lmax      = %d\n" % self.get_lmax()
        s += "cores     = %d\n" % self.get_cores()
        s += "precision = %g\n" % self.get_precision()
        s += "threshold = %g\n" % self.casimir.threshold
        s += "detalg    = %s" % self.get_detalg()
        return s


    def lnLambda(Casimir self, int l1, int l2, int m):
        """Calculate log(|Λ(l1,l2,m)|)"""
        return casimir_lnLambda(l1, l2, m)

    def lnab0(Casimir self, int l):
        """Calculate prefactors log(a0), log(b0) of Mie coefficients for χ→0"""
        cdef double a0,b0
        cdef sign_t sign_a0, sign_b0
        casimir_lnab0(l, &a0, &sign_a0, &b0, &sign_b0)
        return a0, sign_a0, b0, sign_b0

    def lnab(Casimir self, int n, int l):
        """Calculate Mie coefficients a_l, b_l

        Calculate Mie coefficients for order l and argument χ=nT*R/(R+L). The
        function returns log(a_l), sign(a_l), log(b_l), sign(b_l).
        """
        cdef double lna,lnb
        cdef sign_t sign_a, sign_b
        casimir_lnab(&self.casimir, n, l, &lna, &lnb, &sign_a, &sign_b)
        return lna, sign_a, lnb, sign_b

    def lnab_perf(Casimir self, int n, int l):
        """Calculate Mie coefficients for perfect reflectors"""
        cdef double lna,lnb
        cdef sign_t sign_a, sign_b
        casimir_lnab_perf(&self.casimir, n, l, &lna, &lnb, &sign_a, &sign_b)
        return lna, sign_a, lnb, sign_b


    def rp(Casimir self, double nT, double k):
        """Calculate and return Fresnel coefficients r_TE, r_TM"""
        cdef double r_TE, r_TM
        casimir_rp(&self.casimir, nT, k, &r_TE, &r_TM)
        return r_TE, r_TM


    def logdetD0(Casimir self, int m):
        """Calculate determinant of scattering matrix in the high-temperature limit, i.e., xi=0"""
        cdef double EE, MM
        casimir_logdetD0(&self.casimir, m, &EE, &MM)
        return EE, MM

    def logdetD(Casimir self, int n, int m):
        """Calculate determinant of scattering matrix"""
        return casimir_logdetD(&self.casimir, n, m)


cdef double __epsilonm1(double xi, void *userdata):
    """Wrapper for callbacks"""
    cdef Casimir self

    self = <Casimir> userdata
    return self.epsilonm1(xi, self.args)

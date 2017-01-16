import numpy as np
cimport numpy as np
import cython
import math
from libcpp cimport bool
from libc.stdlib cimport malloc, free

ctypedef signed char sign_t

cdef extern from "matrix.h":
    ctypedef enum detalg_t:
         DETALG_LU, DETALG_QR, DETALG_EIG

    ctypedef struct matrix_t:
        size_t dim,dim2
        size_t lda
        double *M;
        bool free_memory

    matrix_t *matrix_alloc(const size_t dim)
    matrix_t *matrix_view(double *ptr, size_t dim, size_t lda)

    void matrix_free(matrix_t *A)
#    void matrix_setall(matrix_t *A, double z)
#
#    int matrix_save_to_stream(matrix_t *A, FILE *stream)
#    int matrix_save_to_file(matrix_t *A, const char *filename)
#
#    double matrix_trace(matrix_t *A)
#
#    double matrix_logdet_triangular(matrix_t *A)
#    double matrix_logdet(matrix_t *A, double z, detalg_t detalg)
#    double matrix_logdet_lu(matrix_t *A)
#    double matrix_logdet_qr(matrix_t *A)
#    double matrix_logdetIdmM_eig(matrix_t *A, double z);


cdef extern from "libcasimir.h":
    ctypedef struct casimir_t:
        double LbyR
        double threshold
        double tolerance
        int ldim

    ctypedef enum polarization_t:
         TE, TM

    casimir_t *casimir_init(double LbyR)
    void casimir_free(casimir_t *self)

    void casimir_set_debug(casimir_t *self, bool debug)
    bool casimir_get_debug(casimir_t *self)

    void casimir_set_verbose(casimir_t *self, bool verbose)
    bool casimir_get_verbose(casimir_t *self)

    int casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi, void *userdata), void *userdata)

    int casimir_get_ldim(casimir_t *self)
    int casimir_set_ldim(casimir_t *self, int ldim)

    detalg_t casimir_get_detalg(casimir_t *self)
    void casimir_set_detalg(casimir_t *self, detalg_t detalg)

    double casimir_get_tolerance(casimir_t *self);
    int    casimir_set_tolerance(casimir_t *self, double tolerance);

    double casimir_lnLambda(int l1, int l2, int m)

    void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0)
    void casimir_lnab(casimir_t *self, double nT, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)
    void casimir_lnab_perf(casimir_t *self, double nT, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)

    void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM)

    void casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM)
    double casimir_logdetD(casimir_t *self, double nT, int m)

    matrix_t *casimir_M(casimir_t *self, double nT, int m)


cdef extern from "integration.h":
    ctypedef struct integration_t

    integration_t *casimir_integrate_init(casimir_t *self, double xi, int m, double epsrel);
    void casimir_integrate_free(integration_t *integration);

    double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
    double casimir_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign);

    double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
    double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
    double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
    double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);


cdef extern from "sfunc.h":
    double lfac(unsigned int n);
    double logi(unsigned int x);

    double ln_factorial2(unsigned int n);

    void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p);

    double Plm(int l, int m, double x, double factor, int mode);
    double Plm_estimate(int l, int m, double x);


cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

class sfunc:
    """Special functions

    This module provides some special functions that are needed to evaluate the
    Casimir free energy in the plane-sphere geometry. Special functions
    implemented include modified Bessel functions for half-integer orders,
    associated Legendre polynomials, and the double factorial.
    """

    def ln_factorial2(unsigned int n):
        """Logarithm of double factorial

        Logarithm of double factorial log(n!!), where n!! is the factorial with
        every second value skipped, i.e.,
            log(7!! = 7*5*3*1.
        """
        return ln_factorial2(n)

    def lfac(unsigned n):
        """Calculate log(n!)"""
        return lfac(n)

    def logi(unsigned n):
        """Calculate log(n) for n integer"""
        return logi(n)


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
        I,K = sfunc.lnInuKnu(nu,x)
        return I

    def lnKnu(int nu, double x):
        """Modified Bessel function K_{nu+1/2}(x) (see lnInuKnu)"""
        I,K = sfunc.lnInuKnu(nu,x)
        return K

    def Plm_estimate(int l, int m, double x):
        return Plm_estimate(l, m, x)

    def Plm(int l, int m, double x, double factor=1, int mode=0):
        return Plm(l,m,x,factor,mode)


cdef class Casimir:
    cdef casimir_t *casimir
    cpdef args
    cpdef epsilonm1

    def __init__(Casimir self, double LbyR, verbose=False, debug=False, ldim=None, tolerance=None, detalg="LU"):
        """Initialize Casimir object

        Required argument:
            LbyR: aspect ration L/R, LbyR > 0

        Optional arguments:
            ldim:      truncation of vector space
            tolerance  relative error for numerical integration
            detalg:    LU, QR, EIG
            threshold: threshold for matrix elements
            debug:     flag, print debugging information
            verbose:   flag, print some addition information
        """
        if LbyR <= 0:
            self.casimir = NULL
            raise ValueError("invalid value for LbyR")

        self.casimir = casimir_init(LbyR)
        if self.casimir == NULL:
            raise RuntimeError("casimir_init returned NULL.")

        if ldim != None:
            self.set_ldim(ldim)
        if tolerance != None:
            self.set_tolerance(tolerance)
        if detalg != None:
            self.set_detalg(detalg)

        self.set_debug(debug)
        self.set_verbose(verbose)

    def __dealloc__(self):
        if self.casimir != NULL:
            casimir_free(self.casimir)

    def get_debug(self):
        return casimir_get_debug(self.casimir)

    def set_debug(self, bool debug):
        return casimir_set_debug(self.casimir, debug)

    def get_verbose(self):
        return casimir_get_verbose(self.casimir)

    def set_verbose(self, bool debug):
        return casimir_set_verbose(self.casimir, debug)

    def get_ldim(self):
        return casimir_get_ldim(self.casimir)

    def set_ldim(self, int ldim):
        if ldim < 1:
            raise ValueError("invalid value for ldim")
        return casimir_set_ldim(self.casimir, ldim)

    def get_detalg(self):
        return casimir_get_detalg(self.casimir)

    def set_detalg(self, detalg):
        if detalg.upper() == "LU":
            casimir_set_detalg(self.casimir, DETALG_LU)
        elif detalg.upper() == "QR":
            casimir_set_detalg(self.casimir, DETALG_QR)
        elif detalg.upper() == "EIG":
            casimir_set_detalg(self.casimir, DETALG_EIG)
        else:
            raise ValueError("invalid value for detalg")

    def get_tolerance(self):
        return casimir_get_tolerance(self.casimir)

    def set_tolerance(self, double tolerance):
        if tolerance <= 0:
            raise ValueError("invalid value for tolerance")
        return casimir_set_tolerance(self.casimir, tolerance)

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
        casimir_set_epsilonm1(self.casimir, &__epsilonm1, ptr)

    def __repr__(self):
        s  = "L/R       = %.8g\n" % self.casimir.LbyR
        s += "ldim      = %d\n" % self.get_ldim()
        s += "tolerance = %g\n" % self.get_tolerance()
        s += "threshold = %g\n" % self.casimir.threshold
        s += "detalg    = %d\n" % self.get_detalg()
        s += "verbose   = %s\n" % ("True" if self.get_verbose() else "False")
        s += "debug     = %s"   % ("True" if self.get_debug()   else "False")
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

    def lnab(Casimir self, double nT, int l):
        """Calculate Mie coefficients a_l, b_l

        Calculate Mie coefficients for order l and argument χ=nT*R/(R+L). The
        function returns log(a_l), sign(a_l), log(b_l), sign(b_l).
        """
        cdef double lna,lnb
        cdef sign_t sign_a, sign_b
        casimir_lnab(self.casimir, nT, l, &lna, &lnb, &sign_a, &sign_b);
        return lna, sign_a, lnb, sign_b

    def lnab_perf(Casimir self, double nT, int l):
        """Calculate Mie coefficients for perfect reflectors"""
        cdef double lna,lnb
        cdef sign_t sign_a, sign_b
        casimir_lnab_perf(self.casimir, nT, l, &lna, &lnb, &sign_a, &sign_b)
        return lna, sign_a, lnb, sign_b


    def rp(Casimir self, double nT, double k):
        """Calculate and return Fresnel coefficients r_TE, r_TM"""
        cdef double r_TE, r_TM
        casimir_rp(self.casimir, nT, k, &r_TE, &r_TM)
        return r_TE, r_TM


    def logdetD0(Casimir self, int m):
        """Calculate determinant of scattering matrix in the high-temperature limit, i.e., xi=0"""
        cdef double EE, MM
        casimir_logdetD0(self.casimir, m, &EE, &MM)
        return EE, MM

    def logdetD(Casimir self, double nT, int m):
        """Calculate determinant of scattering matrix"""
        return casimir_logdetD(self.casimir, nT, m)

    def M(Casimir self, double nT, int m):
        raise NotImplementedError()

    def get_integration(Casimir self, double nT, int m):
        return Integration(self, nT, m, epsrel=self.get_tolerance())

cdef polarization_t __polarization(p):
    if p.upper() == "TE":
        return TE
    else:
        return TM

cdef class Integration:
    cdef double xi,epsrel
    cdef int m
    cdef integration_t *integration

    def __init__(Integration self, Casimir casimir, xi, m, epsrel=1e-8):
        self.xi = xi
        self.m  = m
        self.epsrel = epsrel
        self.integration = casimir_integrate_init(casimir.casimir, xi, m, epsrel)
        if self.integration == NULL:
            raise RuntimeError("casimir_integrate_init returned NULL.")

    def __dealloc__(Integration self):
         if self.integration != NULL:
             casimir_integrate_free(self.integration)

    def __repr__(Integration self):
        s  = "xi     = %g\n" % self.xi
        s += "m      = %d\n" % self.m
        s += "epsrel = %g" % self.epsrel
        return s

    def I(Integration self, int l1, int l2, p):
        cdef sign_t sign
        v = casimir_integrate_I(self.integration, l1, l2, __polarization(p), &sign)
        return v,sign

    def K(Integration self, int nu, p):
        cdef sign_t sign
        v = casimir_integrate_K(self.integration, nu, __polarization(p), &sign)
        return v,sign

    def A(Integration self, int l1, int l2, p):
        cdef sign_t sign
        v = casimir_integrate_A(self.integration, l1, l2, __polarization(p), &sign)
        return v,sign

    def B(Integration self, int l1, int l2, p):
        cdef sign_t sign
        v = casimir_integrate_B(self.integration, l1, l2, __polarization(p), &sign)
        return v,sign

    def C(Integration self, int l1, int l2, p):
        cdef sign_t sign
        v = casimir_integrate_C(self.integration, l1, l2, __polarization(p), &sign)
        return v,sign

    def D(Integration self, int l1, int l2, p):
        cdef sign_t sign
        v = casimir_integrate_D(self.integration, l1, l2, __polarization(p), &sign)
        return v,sign


cdef double __epsilonm1(double xi, void *userdata):
    """Wrapper for callbacks"""
    cdef Casimir self

    self = <Casimir> userdata
    return self.epsilonm1(xi, self.args)

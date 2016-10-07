from unittest import main, TestCase
import numpy as np
import mpmath as mp
import libcasimir
from math import log,e
from scipy.misc import factorial2
import itertools

# increase accuracy of mpmath
#mp.mp.dps = 500

class TestBessel(TestCase):
    """Test Bessel functions"""

    def _test(self,nu,x, eps=1e-12):
        lnInu, lnKnu = libcasimir.sfunc.lnInuKnu(nu,x)
        lnInu_mp = mp.log(mp.besseli(nu+0.5,x))
        lnKnu_mp = mp.log(mp.besselk(nu+0.5,x))

        errInu = abs((lnInu-lnInu_mp)/lnInu_mp)
        errKnu = abs((lnKnu-lnKnu_mp)/lnKnu_mp)

        if errInu > eps:
            raise AssertionError("lnInu: %.20g != %.20g (mpmath), err=%g>%g, nu=%g, x=%g" % (lnInu, lnInu_mp, errInu, eps, nu, x))
        if errKnu > eps:
            raise AssertionError("lnKnu: %.20g != %.20g (mpmath), err=%g>%g, nu=%g, x=%g" % (lnKnu, lnKnu_mp, errKnu, eps, nu, x))


    def test_bessel_mpmath_small(self):
        lnInuKnu = libcasimir.sfunc.lnInuKnu

        all_upto = 20
        skip     = 71
        maximum  = 1500

        iter_nu = itertools.chain(range(all_upto), range(all_upto, maximum, skip))
        iter_x  = np.logspace(log(1e-6), log(1e6), 100, base=e)

        for nu in iter_nu:
            for x in iter_x:
                self._test(nu,x)


class TestDoublefact(TestCase):
    """Test double factorial n!!"""

    def test_scipy(self):
        ln_doublefact = libcasimir.sfunc.ln_doublefact

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        every = 200     # test all integers up to 200
        skip  = 13      # then only test every 13-th integer
        maximum = 15000 # up to 15000

        for i in itertools.chain(range(0, every), range(every, maximum+1, skip)):
            x = ln_doublefact(i)
            y = mp.log(mp.mpf(factorial2(i, exact=True)))
            assertAlmostEqual(x,y)


class TestLambda(TestCase):
    casimir = libcasimir.Casimir()

    def lnLambda_mpmath(self,l1,l2,m):
        num   = (2*l1+1)*(2*l2+1)*mp.factorial(l1-m)*mp.factorial(l2-m)
        denom = mp.factorial(l1+m)*mp.factorial(l2+m)*l1*(l1+1)*l2*(l2+1)
        return mp.log(mp.sqrt( num/denom ))

        
    def test_mpmath_small(self):
        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=20)

        lnLambda = self.casimir.lnLambda
        lnLambda_mpmath = self.lnLambda_mpmath

        for l1 in range(1,10):
            for l2 in range(1,10):
                lmin = min(l1,l2)
                for m in range(lmin+1):
                    assertAlmostEqual(lnLambda(l1,l2,m), lnLambda_mpmath(l1,l2,m))


    def test_mpmath_arbitrary(self):
        casimir = self.casimir
        lmin = 1      # start from l1,l2 = 1
        lmax = 15000 # up to l1,l2=lmax
        N = 50

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        lnLambda = self.casimir.lnLambda
        lnLambda_mpmath = self.lnLambda_mpmath

        for l1 in map(int, np.logspace(log(lmin), log(lmax), N, base=e)):
            for l2 in map(int, np.logspace(log(lmin), log(lmax), N, base=e)):
                min_l1l2 = min(l1,l2)
                for m in map(int, np.logspace(log(lmin), log(min_l1l2), N, base=e)):
                    assertAlmostEqual(lnLambda(l1,l2,m), lnLambda_mpmath(l1,l2,m))


if __name__ == "__main__":
    main()

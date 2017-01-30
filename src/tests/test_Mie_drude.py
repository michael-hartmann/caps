#!/usr/bin/python3

"""This program creates the test for Mie coefficients with Drude metals."""

from __future__ import division
from mpmath import *
from numpy import linspace
from sys import stderr

mp.dps = 40

def prettyprint(x, length=17):
    """Return mpf float in a nice format"""
    return str(x)[:length]


def epsilon(xi, omegap, gamma):
    """Return dielectric function"""
    return 1 + omegap**2/(xi*(xi+gamma))


def lnab(nT,l,LbyR,omegap,gamma):
    """Calculate Mie coefficients al, bl"""
    nT = mpf(nT)
    LbyR = mpf(LbyR)
    omegap = mpf(omegap)
    gamma = mpf(gamma)

    n2 = epsilon(nT, omegap, gamma)
    n  = sqrt(n2)
    chi  = nT/(1+LbyR)
    nchi = n*chi

    besseli_plus_chi   = besseli(l+0.5, chi)
    besseli_plus_nchi  = besseli(l+0.5, nchi)

    besseli_minus_chi  = besseli(l-0.5, chi)
    besseli_minus_nchi = besseli(l-0.5, nchi)

    besselk_plus_chi  = besselk(l+0.5, chi)
    besselk_minus_chi = besselk(l-0.5, chi)

    sla = besseli_plus_nchi * ( l*besseli_plus_chi  -  chi*besseli_minus_chi )
    slb = besseli_plus_chi  * ( l*besseli_plus_nchi - nchi*besseli_minus_nchi )
    slc = besseli_plus_nchi * ( l*besselk_plus_chi  +  chi*besselk_minus_chi )
    sld = besselk_plus_chi  * ( l*besseli_plus_nchi - nchi*besseli_minus_nchi )

    al = (-1)**(l+1) * pi/2 * (n2*sla-slb)/(n2*slc-sld)
    bl = (-1)**(l+1) * pi/2 * (sla-slb)/(slc-sld)

    return log(abs(al)), sign(al), log(abs(bl)), sign(bl)


if __name__ == "__main__":
    LbyR = 1

    print("/* This code was created by test_Mie_drude.py */")
    print("#include \"sfunc.h\"")
    print("#include \"libcasimir.h\"")
    print("#include \"unittest.h\"")
    print()
    print("#include \"test_mie_drude.h\"")
    print()
    print("int test_mie_drude(void)")
    print("{")
    print("    double lna, lnb;")
    print("    double userdata[2];");
    print("    sign_t sign_a, sign_b;")
    print("    casimir_t *casimir;")
    print("    unittest_t test;")
    print()
    print("    unittest_init(&test, \"Mie (Drude)\", \"Test Mie coefficients for various parameters\");")
    print()
    print("    casimir = casimir_init(%g);" % LbyR)


    for omegap,gamma in ((500,1), (100,1), (50,1), (1,1)):
        print("    userdata[0] = %g; userdata[1] = %g;" % (omegap,gamma))
        print("    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);")
        print()
        for nT in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 10, 100, 1000, 10000, 100000):
            for l in (1, 5, 10, 100, 500, 1000, 5000):
                try:
                    lna, sign_a, lnb, sign_b = lnab(mpf(nT),l,mpf(LbyR),mpf(omegap),mpf(gamma))

                    print("    casimir_lnab(casimir, %g, %d, &lna, &lnb, &sign_a, &sign_b);" % (nT,l))
                    print("    AssertEqual(&test, sign_a, %d);" % sign_a)
                    print("    AssertEqual(&test, sign_b, %d);" % sign_b)
                    print("    AssertAlmostEqual(&test, lna, %s);" % prettyprint(lna))
                    print("    AssertAlmostEqual(&test, lnb, %s);" % prettyprint(lnb))
                    print()
                except (ValueError, libmp.libhyper.NoConvergence):
                    print("/*")
                    print("    casimir_lnab(casimir, %g, %d, &lna, &lnb, &sign_a, &sign_b);" % (nT,l))
                    print("    AssertEqual(&test, sign_a, %d);" % sign_a)
                    print("    AssertEqual(&test, sign_b, %d);" % sign_b)
                    print("    AssertAlmostEqual(&test, lna, %s);" % prettyprint(lna))
                    print("    AssertAlmostEqual(&test, lnb, %s);" % prettyprint(lnb))
                    print("*/")
                    print()
        print()

    print("    casimir_free(casimir);")
    print()
    print("    return test_results(&test, stderr);")
    print("}")

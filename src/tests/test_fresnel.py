#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 50

def fresnel_rp(nT, k, omegap, gamma):
    eps = 1 + omegap**2/(nT*(nT+gamma))
    beta = sqrt(1 + (eps-1)/(1 + (k/nT)**2))

    r_TE = (1-beta)/(1+beta)
    r_TM = (eps-beta)/(eps+beta)

    return r_TE, r_TM


if __name__ == "__main__":

    print("#include \"sfunc.h\"")
    print("#include \"unittest.h\"")
    print()
    print("#include \"test_fresnel.h\"")
    print()
    print("int test_fresnel()")
    print("{")
    print("    float80 r_TE, r_TM;")
    print("    unittest_t test;")
    print("    double omegap, gamma_;")
    print("    casimir_t casimir;")
    print()
    print("    unittest_init(&test, \"Fresnel\", \"Test Fresnel coefficients\");");
    print()
    for omegap,gamma in ((500,1), (250,1), (50,1), (5,1), (1,1)):
        print("    omegap = %g;" % omegap)
        print("    gamma_ = %g;" % gamma)
        print("    casimir_init(&casimir, 1, 1); // L/R, T don't matter")
        print("    casimir_set_drude(&casimir, omegap, gamma_, omegap, gamma_);")
        print()
        for nT in (1e-5, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e5):
            for k in (1e-5, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e5):
                r_TE, r_TM = fresnel_rp(mpf(nT), mpf(k), mpf(omegap), mpf(gamma))
                print("    casimir_rp(&casimir, %g, %g, &r_TE, &r_TM);" % (nT, k))
                print("    AssertAlmostEqual(&test, r_TE, %.16g);" % float(r_TE))
                print("    AssertAlmostEqual(&test, r_TM, %.16g);" % float(r_TM))
                print()
        print()
        print("    casimir_free(&casimir);")
        print()
    print("    return test_results(&test, stderr);")
    print("}")

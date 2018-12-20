from mpmath import *

mp.dps = 50

c = 299792458 # speed of light in m/s
hbar_eV = 6.582119514e-16


def fresnel_rp(xi, k, omegap_eV, gamma_eV):
    omegap = omegap_eV/hbar_eV # in rad/s
    gamma  = gamma_eV/hbar_eV # in rad/s

    kappa = sqrt((xi/c)**2+k**2)

    epsm1 = omegap**2/(xi*(xi+gamma)) # Drude model
    x = (xi/(c*kappa))**2*epsm1
    beta = sqrt(1 + x)

    r_TE = (1-beta)/(1+beta)
    r_TM = (epsm1+1-beta)/(epsm1+1+beta)

    return r_TE, r_TM


if __name__ == "__main__":
    R = 100e-6
    L = 1e-6

    print("#include \"libcasimir.h\"")
    print("#include \"misc.h\"")
    print("#include \"constants.h\"")
    print("#include \"unittest.h\"")
    print()
    print("#include \"test_fresnel.h\"")
    print()
    print("int test_casimir_fresnel()")
    print("{")
    print("    unittest_t test;")
    print("    casimir_t *casimir;")
    print("    double rTE, rTM;")
    print("    double userdata[2];")
    print()
    print("    unittest_init(&test, \"casimir_fresnel\", \"Fresnel coefficients\", 2e-15);")
    print()
    print("    casimir = casimir_init(%e,%e);" % (R,L))
    print()

    for omegap_eV in (0.01, 0.1, 1, 2, 5, 7, 8, 9, 10, 20, 50, 100, 1000):
        for gamma_eV in (1e-5, 1e-4, 0.001, 0.003, 0.035, 0.05, 0.1, 0.5, 1):
            print("    userdata[0] = %g/CASIMIR_hbar_eV; /* %geV */" % (omegap_eV,omegap_eV))
            print("    userdata[1] = %g/CASIMIR_hbar_eV; /* %geV */" % (gamma_eV,gamma_eV))

            print("    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);")
            print()
            for xi_ in (1e-6, 1e-5, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6):
                xi = xi_*c/(R+L)
                for k_ in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6):
                    k = k_/(R+L)

                    r_TE, r_TM = fresnel_rp(mpf(xi), mpf(k), mpf(omegap_eV), mpf(gamma_eV))
                    print("    casimir_fresnel(casimir, %g, %g, &rTE, &rTM);" % (xi_, k_))
                    print("    AssertAlmostEqual(&test, rTE, %.16g);" % float(r_TE))
                    print("    AssertAlmostEqual(&test, rTM, %.16g);" % float(r_TM))
                    print()
    print()
    print("    casimir_free(casimir);")
    print()
    print("    return test_results(&test, stderr);")
    print("}")

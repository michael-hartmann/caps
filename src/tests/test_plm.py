#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 100

def Plm(l,m,x):
    if m % 2 == 0:
        return legenp(l,m,x).real
    else:
        return legenp(l,m,x).imag

def dPlm(l,m,x):
    return ((l-m+1)*Plm(l+1,m,x) - (l+1)*x*Plm(l,m,x) )/(x**2-1)


if __name__ == "__main__":
    print("#include <math.h>")

    print("#include \"unittest.h\"")
    print("#include \"sfunc.h\"")
    print("#include \"libcasimir.h\"")
    print()
    print("#include \"test_Plm.h\"")
    print()
    print("int test_Plm()")
    print("{")
    print("    sign_t sign;")
    print("    unittest_t test;")
    print()
    print("    unittest_init(&test, \"Plm\", \"Test associated Legendre polynomials\");")
    print()

    for l in (1, 5, 10, 50, 100, 500, 1000, 2000, 5000):
        for m in (0,1,5,10,50,100):
            for x in (1.1, 2, 5, 10, 100, 1000, 1e6):
                if l >= m:
                    value = Plm(l,m,mpf(x))
                    print("    AssertAlmostEqual(&test, plm_lnPlm(%d,%d,%g, &sign), %s);" % (l,m,x,float(log(abs(value)))))
                    print("    AssertEqual(&test, sign, %d);" % sign(value))

                    value = dPlm(l,m,mpf(x))
                    print("    AssertAlmostEqual(&test, plm_lndPlm(%d,%d,%g, &sign), %s);" % (l,m,x,float(log(abs(value)))))
                    print("    AssertEqual(&test, sign, %d);" % sign(value))
                    print()

    print()
    print("    return test_results(&test, stderr);")
    print("}")

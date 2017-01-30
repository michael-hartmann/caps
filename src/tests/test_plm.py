#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 50

def Plm(l,m,x):
    if m % 2 == 0:
        return legenp(l,m,x).real
    else:
        return legenp(l,m,x).imag


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
    print("    unittest_t test;")
    print()
    print("    unittest_init(&test, \"Plm\", \"Test associated Legendre polynomials\");")
    print()

    for l in (1, 5, 10, 50, 100, 500, 1000, 2000, 5000, 10000, 15000, 20000, 25000):
        for m in (0,1,5,10,50,100,500,1000):
            for x in (1.01, 1.1, 1.5, 2, 5, 10, 50, 100, 500, 1000, 1e6):
                if l >= m:
                    try:
                        value = Plm(l,m,mpf(x))
                        print("    AssertAlmostEqual(&test, Plm(%d,%d,%g,1,1), %s);" % (l,m,x,float(log(abs(value)))))
                    except:
                        print("    // AssertAlmostEqual(&test, Plm(%d,%d,%g,1,1), ???);" % (l,m,x))
        print()

    print()
    print("    return test_results(&test, stderr);")
    print("}")

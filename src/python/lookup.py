from scipy.special import factorial2
from math import log,lgamma

max_logi = 32768
max_lfac = 32768
max_factorial2 = 256

header = """/* created by lookup.py */
#ifndef LOOKUP_H
#define LOOKUP_H

extern double lookup_logi[];
extern size_t lookup_logi_elems;

extern double lookup_factorial2[];
extern size_t lookup_factorial2_elems;

extern double lookup_lfac[];
extern size_t lookup_lfac_elems;
#endif"""

includes = """/* created by lookup.py */
#include <stdlib.h>
#include <math.h>
"""

with open("../lookup.h", "w") as f:
    print(header, file=f)


with open("../lookup.c", "w") as f:
    print(includes, file=f)
    
    print("""
double lookup_logi[] = { /* %g kb */
    -INFINITY, /* log(0) */""" % (max_logi*8/1024), file=f)

    for i in range(1,max_logi-1):
        print("    %.20g, /* log(%d) */" % (log(i),i), file=f)

    print("    %.20g /* log(%d) */" % (log(max_logi-1),max_logi-1), file=f)
    print("};\n", file=f);

    print("size_t lookup_logi_elems = sizeof(lookup_logi)/sizeof(lookup_logi[0]);", file=f)


    print("\n", file=f)


    print("double lookup_factorial2[] = { /* %g kb */" % (max_factorial2*8/1024), file=f)

    for i in range(max_factorial2-1):
        print("    %.20e, /* %d!! */" % (factorial2(i, exact=True),i), file=f)

    print("    %.20e /* %d!! */" % (factorial2(max_factorial2-1, exact=True),max_factorial2-1), file=f)
    print("};\n", file=f);

    print("size_t lookup_factorial2_elems = sizeof(lookup_factorial2)/sizeof(lookup_factorial2[0]);", file=f)


    print("\n", file=f)


    print("double lookup_lfac[] = { /* %g kb */ " % (max_lfac*8/1024), file=f)

    for i in range(max_lfac-1):
        print("    %.20e, /* lgamma(1+%d) */" % (lgamma(1+i),i), file=f)

    print("    %.20e /* lgamma(1+%d) */" % (lgamma(1+max_lfac-1),max_lfac-1), file=f)
    print("};\n", file=f);

    print("size_t lookup_lfac_elems = sizeof(lookup_lfac)/sizeof(lookup_lfac[0]);", file=f)

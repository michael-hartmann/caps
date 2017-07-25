import mpmath as mp

# set working precision
mp.dps = 30

def log(i,digits=18):
    return mp.nstr(mp.log(i), digits)

def lgamma(i,digits=18):
    return mp.nstr(mp.loggamma(i), digits)

max_logi = 2**16
max_lfac = 2**16

header = """/* created by logfac.py */
#ifndef LOGFAC_H
#define LOGFAC_H

double lfac(unsigned int n) __attribute__ ((pure));
double logi(unsigned int x) __attribute__ ((pure));

#endif"""

includes = """/* created by logfac.py */
#include <stdlib.h>
#include <math.h>

#include "logfac.h"
"""

with open("../logfac.h", "w") as f:
    print(header, file=f)


with open("../logfac.c", "w") as f:
    print(includes, file=f)
    
    print("""
static double lookup_logi[] = { /* %g kb */
    -INFINITY, /* log(0) */""" % (max_logi*8/1024), file=f)

    for i in range(1,max_logi-1):
        print("    %s, /* log(%d) */" % (log(i),i), file=f)

    print("    %s /* log(%d) */" % (log(max_logi-1),max_logi-1), file=f)
    print("};\n", file=f);

    print("static size_t lookup_logi_elems = sizeof(lookup_logi)/sizeof(lookup_logi[0]);", file=f)


    print("\n", file=f)


    print("static double lookup_lfac[] = { /* %g kb */" % (max_lfac*8/1024), file=f)

    for i in range(max_lfac-1):
        print("    %s, /* lgamma(1+%d) */" % (lgamma(1+i),i), file=f)

    print("    %s /* lgamma(1+%d) */" % (lgamma(1+max_lfac-1),max_lfac-1), file=f)
    print("};\n", file=f);

    print("static size_t lookup_lfac_elems = sizeof(lookup_lfac)/sizeof(lookup_lfac[0]);", file=f)


    print(r"""
/** @brief Calculate log(x) for x integer
 *
 * This function uses a lookup table to avoid calling log() for n "small".
 *
 * @param [in] n integer
 * @retval log log(n)
 */
double logi(unsigned int n)
{
    if(n < lookup_logi_elems)
        return lookup_logi[n];
    else
        return log(n);
}


/** @brief Calculate log(n!) = log(Γ(n+1))
 *
 * @param [in] n integer
 * @retval lfac log(n!) = log(Γ(n+1))
 */
double lfac(unsigned int n)
{
    if(n < lookup_lfac_elems)
        return lookup_lfac[n];
    else
        return lgamma(1+n);
}""", file=f)

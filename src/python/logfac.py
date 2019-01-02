import datetime
import mpmath as mp

# script to create logfac.c and logfac.h

# set working precision
mp.dps = 64

def date():
    today = datetime.date.today()
    d = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    return d[today.month-1], today.year


def log(i,digits=19):
    return mp.nstr(mp.log(i), digits)

def lgamma(i,digits=19):
    return mp.nstr(mp.loggamma(i), digits)

if __name__ == "__main__":
    max_logi = 2**16
    max_lfac = 2**10

    assert max_lfac > 1023

    header = """/* created by logfac.py */
#ifndef LOGFAC_H
#define LOGFAC_H

#ifdef __cplusplus
extern "C" {
#endif

double lfac(unsigned int n) __attribute__ ((pure));
double logi(unsigned int x) __attribute__ ((pure));

double lfac2(unsigned int n) __attribute__ ((pure));

#ifdef __cplusplus
}
#endif

#endif"""

    month,year = date()
    includes = """/* created by logfac.py */

/**
 * @file   logfac.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   %s, %d
 * @brief  computation of logarithm and factorial for integer arguments; created by logfac.py
 */

#include <stdlib.h>
#include <math.h>

#include "logfac.h"\n\n""" % (month, year)

    with open("../include/logfac.h", "w") as f:
        print(header, file=f)

    with open("../logfac.c", "w") as f:
        print(includes, file=f)

        print(r"""/** lookup table for \f$\log(n)\f$, see \ref logi */
static double lookup_logi[] = { /* %g kb */
    -INFINITY, /* log(0) */""" % (max_logi*8/1024), file=f)

        for i in range(1,max_logi-1):
            print("    %s, /* log(%d) */" % (log(i),i), file=f)

        print("    %s /* log(%d) */" % (log(max_logi-1),max_logi-1), file=f)
        print("};\n\n", file=f);


        print(r"""/** lookup table for \f$n!\f$, see \ref lfac */
static double lookup_lfac[] = { /* %g kb */""" % (max_lfac*8/1024), file=f)

        for i in range(max_lfac-1):
            print("    %s, /* lgamma(1+%d) */" % (lgamma(1+i),i), file=f)

        print("    %s /* lgamma(1+%d) */" % (lgamma(1+max_lfac-1),max_lfac-1), file=f)
        print("};\n", file=f);

        print(r"""

const size_t __lookup_logi_elems = sizeof(lookup_logi)/sizeof(lookup_logi[0]);
const size_t __lookup_lfac_elems = sizeof(lookup_lfac)/sizeof(lookup_lfac[0]);

/** @brief Calculate \f$\log(n)\f$ for integer \f$n\f$
 *
 * This function uses a lookup table to avoid calling log() for \f$n \le %d\f$
 *
 * @param [in] n integer
 * @retval logn \f$\log(n)\f$
 */
double logi(unsigned int n)
{
    if(n < __lookup_logi_elems)
        /* use lookup table */
        return lookup_logi[n];
    else
        /* use log function */
        return log(n);
}

/** @brief Calculate \f$\log(n!) = \log(\Gamma(n+1))\f$
 *
 * This function computes the logarithm of the factorial \f$n!\f$. This function uses
 * a lookup table for \f$n \le %d\f$
 *
 * @param [in] n integer
 * @retval lfac \f$\log(n!)\f$
 */
double lfac(unsigned int n)
{
    if(n < __lookup_lfac_elems)
        return lookup_lfac[n];

    /* Use Stirling's approximation to compute the logarithm of the factorial
     * function for integer arguments greater or equal than 1024.
     */
    const double C = 0.91893853320467274178; /* log(2*pi)/2 */
    if(n < __lookup_logi_elems)
        return (n+0.5)*lookup_logi[n] - n + C + (1./12)/n;
    else
        return (n+0.5)*log(n) - n + C + (1./12)/n;
}
""" % (max_logi, max_lfac), file=f)

        print(r"""/** @brief Calculate \f$\log(n!!)\f$
 *
 * This function computes the logarithm of the double factorial \f$n!!\f$.
 *
 * @param [in] n argument
 * @retval lfac2 \f$n!!\f$
 */
double lfac2(unsigned int n)
{
    if(n % 2 == 0)
    {
        /* even */
        unsigned int k = n/2;
        return k*log(2)+lfac(k);
    }
    else
    {
        /* odd */
        unsigned int k = (n+1)/2;
        return lfac(2*k)-k*log(2)-lfac(k);
    }
}""", file=f)

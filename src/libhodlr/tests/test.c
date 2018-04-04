/* Test HODLR library
 *
 * Test of the HODLR library and the wrapper. We use the kernel for the MM
 * polarization block in the high-temperature limit. For m=0 the result is
 * known analytically. [1]
 *
 * We compare the numerical result with the exact expression for different
 * aspect ratios R/L. The dimension of the vector space is chosen as
 * ldim=10*R/L.
 *
 * [1] Giuseppe Bimonte, Classical Casimir interaction of a perfectly
 *     conducting sphere and plate, Phys. Rev. D 95, 065004 (2017)
 *     https://doi.org/10.1103/PhysRevD.95.065004
 */

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>

#include <hodlr.h>

/* ANSI codes for colour */
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"

/* return time in seconds since 1970 */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec*1e-6;
}

/* analytical solution of kernel, see [1] */
double analytical(double RbyL)
{
    const int lmax = 20000;
    const double x = 1/RbyL;
    const double Z = 1/(1+x+sqrt(x*(2+x)));

    double sum = 0;
    for(int l = 0; l < lmax; l++)
        sum += log1p(-pow(Z,2*l+3));

    return sum;
}

/* kernel for high-temperature limit, MM block */
double kernel(int i, int j, void *args)
{
    const double l1 = i+1, l2 = j+1;
    const double y = *(double *)args;

    return exp( (l1+l2+1)*y + lgamma(1+l1+l2) - lgamma(1+l1) - lgamma(1+l2)) *sqrt(l1*l2/((l1+1)*(l2+1)));
}

int main(int argc, char *argv[])
{
    char *status[] = { KGRN "ok" KNRM, KRED "failed" KNRM };
    unsigned int nLeaf = 100;
    double tolerance = 1e-10;
    int is_symmetric = 1;
    double eta = 10;
    double RbyL_start = 100;
    double RbyL_stop  = 100000;
    int npts = 20;

    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    printf("# R/L, analytical, HODLR, relative error, status\n");
    for(int i = 0; i < npts; i++)
    {
        double RbyL = RbyL_start*exp(i*log(RbyL_stop/RbyL_start)/(npts-1));

        double dim = RbyL*eta;
        double y = log(0.5/(1+1./RbyL));
        double t = now();
        double logdet = hodlr_logdet(dim, kernel, &y, nLeaf, tolerance, is_symmetric);
        t = now()-t;
        double exact = analytical(RbyL);

        double relerr = 1-fabs(logdet/exact);
        printf("%g, %.14g, %.14g, %.14g, %.5g, %s\n", RbyL, exact, logdet, relerr, t, relerr < 2e-9 ? status[0] : status[1]);
    }

    return 0;
}

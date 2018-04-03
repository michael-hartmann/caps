#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>

#include <hodlr.h>

/* return time in seconds since 1970 */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec*1e-6;
}

/* analytical solution of kernel */
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

    for(int i = 0; i < npts; i++)
    {
        double RbyL = RbyL_start*exp(i*log(RbyL_stop/RbyL_start)/(npts-1));

        double dim = RbyL*eta;
        double y = log(0.5/(1+1./RbyL));
        double t = now();
        double logdet = hodlr_logdet(dim, kernel, &y, nLeaf, tolerance, is_symmetric);
        t = now()-t;
        double exact = analytical(RbyL);

        double error = 1-fabs(logdet/exact);
        printf("%g, %.14g, %g, %s\n", RbyL, error, t, error < 2e-9 ? "ok" : "failed");
    }

    return 0;
}

#include <math.h>

#include "psd.h"
#include "utils.h"

/* prototype for LAPACK routine */
int dstemr_(char *jobz, char *range, int *n, double *
    d__, double *e, double *vl, double *vu, int *il, 
    int *iu, int *m, double *w, double *z__, int *ldz, 
     int *nzc, int *isuppz, int *tryrac, double *work, 
    int *lwork, int *iwork, int *liwork, int *info);

/* Hu, Xu, Yan, J. Chem. Phys. 133, 101106 (2010) */

/* compute expansion coefficients eta according to the paragraph around
 * equations (12) and (13). */
static double _eta(int N, double z)
{
    double Am = 1/4.; /* A1 */
    double A  = 5/4.; /* A2 */

    double Bm = 3;      /* B1 */
    double B  = 15+z/4; /* B2 */

    double dBm = 0;    /* B1'=dB1/dz */
    double dB  = 1/4.; /* B2'=dB2/dz */

    /* Eq. (12) */
    for(int M = 3; M <= N; M++)
    {
        double Amm = Am, Bmm = Bm, dBmm = dBm;

        Am = A;
        Bm = B;
        dBm = dB;

        dB = (2*M+1)*dBm + Bmm/4 + z*dBmm/4; /* derivative of Eq. (12) */
        A = (2*M+1)*Am + z*Amm/4;            /* Eq. (12) */
        B = (2*M+1)*Bm + z*Bmm/4;            /* Eq. (12) */

        /* rescale to prevent overflows */
        Am /= A;
        Bm /= A;
        B  /= A;
        dBm /= A;
        dB  /= A;
        A = 1;
    }

    return A/(2*dB);
}

/** @brief Compute poles and xi_j and expansion coefficients eta_j for PSD
 *
 * This function computes the poles xi_j (at imaginary frequency) and the
 * expansion coefficients eta_j for the Pade spectrum decomposition of order N,
 * see reference [1].  The poles are stored in the array xi, the coefficients
 * are stored in the array eta.
 *
 * References:
 * [1] Hu, Xu, Yan, J. Chem. Phys. 133, 101106 (2010)
 *
 * @param [in]  N   order
 * @param [out] xi  poles
 * @param [out] eta expansion coefficients
 * @retval success 0 if successful
 */
int psd(int N, double xi[N], double eta[N])
{
    char jobz  = 'N';    /* compute eigenvalues only */
    char range = 'I';    /* find IL-th through IU-th eigenvalues */
    int ndim   = 2*N;    /* dimension of matrix */
    double vl  = 0;      /* not referenced */
    double vu  = 0;      /* not referenced */
    int il = 1;          /* index of smallest eigenvalue to compute */
    int iu = N;          /* index of largest eigenvalue to compute */
    int m = 0;           /* total number of eigenvalues found */
    double *z = NULL;    /* not referenced */
    int ldz = 1;         /* probably not referenced or ignored */
    int nzc = ndim;      /* number of eigenvalues to be hold in Z */
    int *isuppz = NULL;  /* not referenced */
    int tryrac = 1;      /* check for high accuracy */
    int lwork = 12*ndim; /* dimension of array work */
    int liwork = 8*ndim; /* dimension of array iwork */
    int info = 0;        /* status */

    /* allocate memory */
    double *d, *e, *w, *work;
    int *iwork;
    d     = xcalloc(ndim,   sizeof(double)); /* diagonal elements are zero, see Eq. (14) */
    e     = xcalloc(ndim,   sizeof(double)); /* off-diagonal elements, set later; see Eq. (14) */
    w     = xcalloc(ndim,   sizeof(double)); /* eigenvalues in ascending order */
    work  = xcalloc(lwork,  sizeof(double)); /* working array */
    iwork = xcalloc(liwork, sizeof(double)); /* working array */

    /* set off-diagonal elements */
    for(int j = 1; j <= ndim; j++)
        e[j-1] = 1/sqrt((2*j+1)*(2*j+3)); /* Eq. (14) */

    /* compute eigenvalues of a real, symmetric, tridiagonal matrix */
    dstemr_(&jobz, &range, &ndim, d, e, &vl, &vu, &il, &iu, &m, w, z, &ldz, &nzc, isuppz, &tryrac, work, &lwork, iwork, &liwork, &info);

    if(info != 0)
        goto out;

    /* iterate over the first N (negative) eigenvalues */
    for(int i = 0; i < N; i++)
    {
        /* compute ξ_j from eigenvalues w_j: ξ_j = ±2/w_j */
        xi[i] = -2/w[i];

        /* compute η_j, see paragraph around Eq. (12) */
        eta[i] = _eta(2*N, -xi[i]*xi[i]);
    }

    out:
    /* free memory */
    xfree(w);
    xfree(e);
    xfree(d);
    xfree(work);
    xfree(iwork);

    return info;
}

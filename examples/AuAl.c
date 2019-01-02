#include <stdio.h>

#include "libcasimir.h"

/* Dielectric function epsilon(xi)-1 for gold: omegap=9eV, gamma=35meV.
 * The frequency xi is in rad/s.
 */
double epsm1_gold(double xi, void *args)
{
    /* plasma frequency and relaxation frequency in rad/s */
    double omegap = 9/CASIMIR_hbar_eV, gamma = 0.035/CASIMIR_hbar_eV;

    /* epsilon-1 = omegap²/(xi*(xi+gamma)) */
    return omegap*omegap/(xi*(xi+gamma));
}

/* Dielectric function epsilon(xi)-1 for aluminium: omegap=11.5eV, gamma=50meV.
 * The frequency xi is in rad/s.
 */
double epsm1_aluminium(double xi, void *args)
{
    /* plasma frequency and relaxation frequency in rad/s */
    double omegap = 11.5/CASIMIR_hbar_eV, gamma = 0.05/CASIMIR_hbar_eV;

    /* epsilon-1 = omegap²/(xi*(xi+gamma)) */
    return omegap*omegap/(xi*(xi+gamma));
}

int main(int argc, char *argv[])
{
    const double R = 50e-6; /* radius R=50µm */
    const double L = 1e-6;  /* separation L=1µm */
    const double T = 300;   /* temperature T=300K */

    /* initialize casimir object */
    casimir_t *casimir = casimir_init(R,L);

    /* Set dimension of vector space to 6.4*R/L; this gives an estimated error
     * of 1e-5 due the truncation of the vector space, see Table 5.1 in
     * [Hartmann, Casimir effect in the plane-sphere geometry: Beyond the
     * proximity force approximation, phd thesis, 2018]
     */
    const int ldim = 6.4*R/L;
    casimir_set_ldim(casimir, ldim);

    /* we assume Drude metals; the sphere is made out of gold, the plate is
     * made out of aluminium
     */
    casimir_set_epsilonm1_plate (casimir, epsm1_aluminium, NULL);
    casimir_set_epsilonm1_sphere(casimir, epsm1_gold, NULL);

    printf("# We assume Drude metals; sphere is gold, plate is aluminium\n#\n");

    /* print information about the casimir object */
    casimir_info(casimir, stdout, "# ");

    /* compute the high-temperature contribution; F is in units of kB*T */
    double F = casimir_ht_drude(casimir);

    /* compute and sum the contributions from the Matsubara frequencies xi_n
     * with n>0
     */
    for(int n = 1; 1; n++)
    {
        /* n-th Matsubara frequency */
        double xi_n = 2*M_PI*n*CASIMIR_kB*T/CASIMIR_hbar;

        double calLbyc = (L+R)/CASIMIR_c; /* (L+R)c */

        /* compute the sum over m */
        double sum_m = 0;
        for(int m = 0; 1; m++)
        {
            double v = casimir_logdetD(casimir, xi_n*calLbyc, m);

            /* the contribution for +m and -m are identical, so for m>0 we have
             * the same contribution from -m
             */
            if(m)
                v *= 2;

            sum_m += v;

            /* stop summation over m if summation converged */
            if(v/sum_m < 1e-8)
                break;
        }

        F += sum_m;

        /* stop summation over Matsubara frequencies if summation converged */
        if(sum_m/F < 1e-7)
            break;
    }

    /* result */
    printf("free energy: %.8g kb*T\n", F);

    /* free object and release allocated memory */
    casimir_free(casimir);

    return 0;
}

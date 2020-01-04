#include <stdio.h>

#include "libcaps.h"

/* Dielectric function epsilon(xi)-1 for gold: omegap=9eV, gamma=35meV.
 * The frequency xi is in rad/s.
 */
static double epsm1_gold(double xi, void *args)
{
    /* plasma frequency and relaxation frequency in rad/s */
    double omegap = 9/CAPS_HBAR_EV, gamma = 0.035/CAPS_HBAR_EV;

    /* epsilon-1 = omegap²/(xi*(xi+gamma)) */
    return omegap*omegap/(xi*(xi+gamma));
}

/* Dielectric function epsilon(xi)-1 for aluminium: omegap=11.5eV, gamma=50meV.
 * The frequency xi is in rad/s.
 */
static double epsm1_aluminium(double xi, void *args)
{
    /* plasma frequency and relaxation frequency in rad/s */
    double omegap = 11.5/CAPS_HBAR_EV, gamma = 0.05/CAPS_HBAR_EV;

    /* epsilon-1 = omegap²/(xi*(xi+gamma)) */
    return omegap*omegap/(xi*(xi+gamma));
}

int main(int argc, char *argv[])
{
    const double R = 50e-6; /* radius R=50µm */
    const double L = 1e-6;  /* separation L=1µm */
    const double T = 300;   /* temperature T=300K */

    /* initialize caps object */
    caps_t *caps = caps_init(R,L);

    /* Set dimension of vector space to 6.4*R/L; this gives an estimated error
     * of 1e-5 due the truncation of the vector space, see Table 5.1 in
     * [Hartmann, Casimir effect in the plane-sphere geometry: Beyond the
     * proximity force approximation, phd thesis, 2018]
     */
    const int ldim = 6.4*R/L;
    caps_set_ldim(caps, ldim);

    /* we assume Drude metals; the sphere is made out of gold, the plate is
     * made out of aluminium
     */
    caps_set_epsilonm1_plate (caps, epsm1_aluminium, NULL);
    caps_set_epsilonm1_sphere(caps, epsm1_gold, NULL);

    printf("# We assume Drude metals; sphere is gold, plate is aluminium\n#\n");

    /* print information about the caps object */
    caps_info(caps, stdout, "# ");

    /* compute the high-temperature contribution; F is in units of kB*T */
    double F = caps_ht_drude(caps);

    /* compute and sum the contributions from the Matsubara frequencies xi_n
     * with n>0
     */
    for(int n = 1; 1; n++)
    {
        /* n-th Matsubara frequency */
        double xi_n = 2*CAPS_PI*n*CAPS_KB*T/CAPS_HBAR;

        double calLbyc = (L+R)/CAPS_C; /* (L+R)c */

        /* compute the sum over m */
        double sum_m = 0;
        for(int m = 0; 1; m++)
        {
            double v = caps_logdetD(caps, xi_n*calLbyc, m);

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
    caps_free(caps);

    return 0;
}

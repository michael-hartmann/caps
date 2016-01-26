#include <stdio.h>
#include <math.h>

#include "libcasimir.h"
#include "integration_drude.h"
#include "sfunc.h"

#define TE 0
#define TM 1

float80 TraceD0(casimir_t *self, int m);
void TraceD(casimir_t *self, int n, int m, float80 Tr_EE[2], float80 Tr_MM[2]);
void usage(FILE *stream, const char *self);

float80 TraceD0(casimir_t *self, int m)
{
    double lnRbyScriptL = log(self->RbyScriptL);
    float80 Tr_EE = 0;

    int min = MAX(m,1);
    int max = self->lmax;

    /* calculate the logarithm of the matrix elements of D */
    for(int l = min; l <= max; l++)
    {
        sign_t sign_a0, sign_b0, sign_xi;
        double lna0, lnb0;
        double lnXiRL = casimir_lnXi(l,l,m,&sign_xi)+(2*l+1)*lnRbyScriptL;
        casimir_lnab0(l, &lna0, &sign_a0, &lnb0, &sign_b0);

        Tr_EE += -sign_xi*sign_a0*exp80(lna0+lnXiRL);
    }

    return Tr_EE;
}

void TraceD(casimir_t *self, int n, int m, float80 Tr_EE[2], float80 Tr_MM[2])
{
    int min,max;
    Tr_EE[TE] = Tr_EE[TM] = 0;
    Tr_MM[TE] = Tr_MM[TM] = 0;

    min = MAX(m,1);
    max = self->lmax;

    if(n == 0)
    {
        Tr_EE[TM] = TraceD0(self, m);
        return;
    }

    for(int l = min; l <= max; l++)
    {
        casimir_integrals_t cint;
        double ln_al, ln_bl;
        sign_t sign_al, sign_bl;

        casimir_mie_cache_get(self, l, n, &ln_al, &sign_al, &ln_bl, &sign_bl);

        casimir_integrate_drude(self, &cint, l, l, m, n, self->T);

        /* EE */
        Tr_EE[TE] += -sign_al*cint.signA_TE*exp80(ln_al+cint.lnA_TE);
        Tr_EE[TM] += -sign_al*cint.signB_TM*exp80(ln_al+cint.lnB_TM);

        /* MM */
        Tr_MM[TE] += -sign_bl*cint.signB_TE*exp80(ln_bl+cint.lnB_TE);
        Tr_MM[TM] += -sign_bl*cint.signA_TM*exp80(ln_bl+cint.lnA_TM);
    }
}

void usage(FILE *stream, const char *self)
{
    fprintf(stream, "Usage: %s D/R T omegap omegap/gamma\n", self);
}


int main(int argc, char *argv[])
{
    casimir_t casimir;
    int lmax = 1;
    float80 Tr = 0;
    double eps = 1e-12;

    /* in beiden Skalierungen identisch */
    double DbyR; /* DbyR = ScriptL/R = 1/RbyScriptL */
    double T;

    /* Skalierung Stefan */
    double omegap_stefan; /* omega_pl * R / 2pi c */
    double omegapbygamma; /* omega_p/gamma */

    if(argc != 5)
    {
        usage(stderr, argv[0]);
        exit(1);
    }
    DbyR          = atof(argv[1]);
    T             = atof(argv[2]);
    omegap_stefan = atof(argv[3]);
    omegapbygamma = atof(argv[4]);

    printf("# DbyR=%g, T*2pi*kb*D/(hbar*c)=%g, omegap*R/(2pi*c) = %g, omegap/gamma = %g\n", DbyR, T, omegap_stefan, omegapbygamma);
    /* Umrechnen in meine Skalierung */
    double omegap = 2*M_PI*omegap_stefan*DbyR;
    double gamma_ = omegap/omegapbygamma;

    {
        casimir_init(&casimir, DbyR-1, T);

        /* set lmax */
        casimir_set_lmax(&casimir, lmax);

        /* set omega_p */
        casimir_set_omegap_sphere(&casimir, omegap);
        casimir_set_omegap_plane (&casimir, omegap);

        /* set gamma */
        casimir_set_gamma_sphere(&casimir, gamma_);
        casimir_set_gamma_plane (&casimir, gamma_);

        /* set order of Gauss-Laguerre integration */
        casimir_set_integration(&casimir, 150);

        casimir_info(&casimir, stdout, "# ");

        printf("#\n");
        printf("# n, m, P, p, Tr(M)\n");

        {
            Tr = 0;
            float80 Tr_0 = 0;

            for(int n = 0; 1; n++)
            {
                float80 Tr_n = 0;

                for(int m = 0; m <= lmax; m++)
                {
                    float80 Tr_EE[2], Tr_MM[2], sum;

                    TraceD(&casimir, n, m, Tr_EE, Tr_MM);
                    sum = Tr_EE[0] + Tr_EE[1] + Tr_MM[0] + Tr_MM[1];

                    printf("%d, %d, EE, TE, %.15Lg\n", n, m, Tr_EE[TE]);
                    printf("%d, %d, EE, TM, %.15Lg\n", n, m, Tr_EE[TM]);
                    printf("%d, %d, MM, TE, %.15Lg\n", n, m, Tr_MM[TE]);
                    printf("%d, %d, MM, TM, %.15Lg\n", n, m, Tr_MM[TM]);

                    if(m == 0)
                        sum /= 2;

                    Tr_n += sum;
                }

                if(n == 0)
                {
                    Tr_0 = Tr_n;
                    Tr = Tr_n/2;
                }
                else
                {
                    Tr += Tr_n;
                    if(Tr_n/Tr_0 < eps)
                        break;
                }
            }
        }
    }

    return 0;
}

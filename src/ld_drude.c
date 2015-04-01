#include <stdio.h>
#include <math.h>

#include "libcasimir.h"
#include "integration_drude.h"
#include "sfunc.h"

#define TE 0
#define TM 1

edouble TraceD0(casimir_t *self, int m);
void TraceD(casimir_t *self, int n, int m, edouble Tr_EE[2], edouble Tr_MM[2]);
void usage(FILE *stream, int argc, char *argv[]);

edouble TraceD0(casimir_t *self, int m)
{
    int l,min,max;
    double lnRbyScriptL = log(self->RbyScriptL);
    edouble Tr_EE = 0;

    min = MAX(m,1);
    max = self->lmax;

    /* calculate the logarithm of the matrix elements of D */
    for(l = min; l <= max; l++)
    {
        sign_t sign_a0, sign_b0, sign_xi;
        double lna0, lnb0;
        double lnXiRL = casimir_lnXi(l,l,m,&sign_xi)+(2*l+1)*lnRbyScriptL;
        casimir_lnab0(l, &lna0, &sign_a0, &lnb0, &sign_b0);

        Tr_EE += -sign_xi*sign_a0*expq(lna0+lnXiRL);
    }

    return Tr_EE;
}

void TraceD(casimir_t *self, int n, int m, edouble Tr_EE[2], edouble Tr_MM[2])
{
    int min,max,l;
    Tr_EE[TE] = Tr_EE[TM] = 0;
    Tr_MM[TE] = Tr_MM[TM] = 0;

    min = MAX(m,1);
    max = self->lmax;

    if(n == 0)
    {
        Tr_EE[TM] = TraceD0(self, m);
        return;
    }

    for(l = min; l <= max; l++)
    {
        casimir_integrals_t cint;
        double ln_al, ln_bl;
        sign_t sign_al, sign_bl;

        casimir_mie_cache_get(self, l, n, &ln_al, &sign_al, &ln_bl, &sign_bl);

        casimir_integrate_drude(self, &cint, l, l, m, n, self->T);

        /* EE */
        Tr_EE[TE] += -sign_al*cint.signA_TE*expq(ln_al+cint.lnA_TE);
        Tr_EE[TM] += -sign_al*cint.signB_TM*expq(ln_al+cint.lnB_TM);

        /* MM */
        Tr_MM[TE] += -sign_bl*cint.signB_TE*expq(ln_bl+cint.lnB_TE);
        Tr_MM[TM] += -sign_bl*cint.signA_TM*expq(ln_bl+cint.lnA_TM);
    }
}

void usage(FILE *stream, int argc, char *argv[])
{
    fprintf(stream, "Usage: %s D/R T omegap omegap/gamma\n", argv[0]);
}


int main(int argc, char *argv[])
{
    casimir_t casimir;
    int n,m;
    int lmax = 1;
    edouble Tr = 0;
    double eps = 1e-12;

    /* in beiden Skalierungen identisch */
    double DbyR; /* DbyR = ScriptL/R = 1/RbyScriptL */
    double T;

    /* Skalierung Stefan */
    double omegap_stefan; /* omega_pl * R / 2pi c */
    double omegapbygamma; /* omega_p/gamma */

    if(argc != 5)
    {
        usage(stderr, argc, argv);
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
        casimir_init(&casimir, 1./DbyR, T);

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
            edouble Tr_0 = 0;

            for(n = 0; 1; n++)
            {
                edouble Tr_n = 0;

                for(m = 0; m <= lmax; m++)
                {
                    edouble Tr_EE[2], Tr_MM[2], sum;

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

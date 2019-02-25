#include <algorithm>
#include <string>

#include <stdio.h>
#include <strings.h>

#include "argparse.h"
#include "bessel.h"
#include "matrix.h"

#include "cylinder.h"

#include "cquadpack/include/cquadpack.h"

static double __kernel(int i, int j, void *args_);
static double __integrand_dirichlet(double x, void *args);
static double __integrand_neumann(double x, void *args);

class CasimirCP {
    double R, d, calL;
    int lmax;
    detalg_t detalg;

    public:
        // constructor
        CasimirCP(double R, double d, detalg_t detalg=DETALG_HODLR) {
            this->R = R;
            this->d = d;
            this->calL = R+d;
            this->lmax = std::max(25, (int)(5*R/d));
            this->detalg = detalg;
        }

        double get_R() const { return R; }
        double get_d() const { return d; }
        double get_calL() const { return calL; }
        int get_lmax() const { return lmax; }
        detalg_t get_detalg() const { return detalg; }

        void set_lmax(int lmax) { this->lmax = lmax; }
        void set_detalg(detalg_t detalg) { this->detalg = detalg; }

        double logdet_dirichlet(double q) {
            return this->logdet(q, 'D');
        }

        double logdet_neumann(double q) {
            return this->logdet(q, 'N');
        }

        double logdet(double q, char DN) {
            const int dim = lmax;
            kernel_args_t args;

            DN = toupper(DN);
            if(DN != 'D')
                DN = 'N';

            // initialize caches
            args.lmax = lmax;
            args.DN = DN;

            args.cache_ratio = new double[lmax+1];
            args.cache_K     = new double[2*(lmax+1)];

            for(int j = 0; j < 2*(lmax+1); j++)
                args.cache_K[j] = bessel_logKn(j, 2*calL*q);

            if(DN == 'D')
            {
                // Dirichlet: I_n(Rq)/K_n(Rq)
                for(int j = 0; j < lmax+1; j++)
                    args.cache_ratio[j] = bessel_logIn(j, R*q)-bessel_logKn(j, R*q);
            }
            else
            {
                // Neumann: I'_n(Rq)/K'_n(Rq)
                double I[3] = { bessel_logI0(R*q), bessel_logI1(R*q), bessel_logIn(2,R*q) };
                double K[3] = { bessel_logK0(R*q), bessel_logK1(R*q), bessel_logKn(2,R*q) };

                // I'_0(Rq)/K'_0(Rq) = -I_1(x)/K_1(x)
                args.cache_ratio[0] = I[1]-K[1];

                for(int j = 1; j < lmax+1; j++)
                {
					/* denom = -2K'_j(x); K'_j(x) = -1/2*[ K_{j+1}(x) + K_{j-1}(x) ] */
					double denom = K[2]+log1p(exp(K[0]-K[2]));

					/* num = 2I'_j(x); I'_j(x) = = 1/2*[ I_{j+1}(x) + I_{j-1}(x) ] = dI */
					double num = I[0]+log1p(exp(I[2]-I[0]));

					args.cache_ratio[j] = num-denom;

                    I[0] = I[1];
                    I[1] = I[2];
                    I[2] = bessel_logIn(j+2,R*q);

                    K[0] = K[1];
                    K[1] = K[2];
                    K[2] = bessel_logKn(j+2,R*q);
                }
            }
            
            // 1- M_00
            double M00 = exp(args.cache_ratio[0]+args.cache_K[0]);
            double log_rho = log1p(-M00);
            args.alpha = 2/(1-M00);

            // compute first determinant
            args.type = 0;
            double logdet1 = kernel_logdet(dim, __kernel, &args, true, detalg);

            // compute second determinant
            args.type = 1;
            double logdet2 = kernel_logdet(dim, __kernel, &args, true, detalg);

            // free memory
            delete [] args.cache_ratio;
            delete [] args.cache_K;

            return log_rho+logdet1+logdet2;
        }

        double energy(char p, double T, double epsrel=1e-8)
        {
            p = toupper(p);
            if(p != 'D')
                p = 'N';

            if(T == 0)
            {
				double integral, abserr;
				int neval, ier;

                if(p == 'D')
				    integral = dqagi(__integrand_dirichlet, 0, 1, 0, epsrel, &abserr, &neval, &ier, this);
                else
				    integral = dqagi(__integrand_neumann, 0, 1, 0, epsrel, &abserr, &neval, &ier, this);

				/* energy in units of hbar*c*L */
				return integral/(4*M_PI*2*d);
            }
            else
                return NAN;
        }
};

static double __kernel(int i, int j, void *args_)
{
    kernel_args_t *args = (kernel_args_t *)args_;

    // µ1,µ2 >= 1
    const int mu1 = i+1, mu2 = j+1;

    double *ratio = args->cache_ratio; // I_n(qR)/K_n(qR)
    double *K     = args->cache_K;     // K_n(2*calL*q)

    double term1 = 0.5*(ratio[mu1]+ratio[mu2]);
    double U = exp(term1 + K[    mu1+mu2 ]);
    double V = exp(term1 + K[abs(mu1-mu2)]);

    if(args->type == 0)
        // log det(1-U+V) = log det(1-(U-V))
        return U-V;

    // log det(1-U-V-2/rho*vv^T) = log det(1-(U+V+α vv^T)) with α=2/ρ
    double vvT = exp(term1+ratio[0]+K[mu1]+K[mu2]);
    return U+V+args->alpha*vvT;
}

static double __integrand_dirichlet(double x, void *args)
{
	CasimirCP *casimir = (CasimirCP *)args;

    double q = x/(2*casimir->get_d());
    return q*casimir->logdet_dirichlet(q);
}

static double __integrand_neumann(double x, void *args)
{
	CasimirCP *casimir = (CasimirCP *)args;

    double q = x/(2*casimir->get_d());
    return q*casimir->logdet_neumann(q);
}

int main(int argc, const char *argv[]) {
    double epsrel = 1e-8, eta = 6, T = 0, R = NAN, d = NAN;
    const char *detalg = NULL;
    int verbose = 0, lmax = 0;

    const char *const usage[] = {
        "casimir_cylinder [options] [[--] args]",
        "casimir_cylinder [options]",
        NULL,
    };

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Mandatory options"),
        OPT_DOUBLE('R', "radius",      &R, "radius of cylinder in m", NULL, 0, 0),
        OPT_DOUBLE('d', "separation",  &d, "seperation between cylinder and plate in m", NULL, 0, 0),
        OPT_GROUP("Further options"),
        OPT_DOUBLE('T', "temperature", &T, "temperature in K (not implemented yet)", NULL, 0, 0),
        OPT_INTEGER('l', "lmax", &lmax, "dimension of vector space", NULL, 0, 0),
        OPT_DOUBLE('e', "epsrel", &epsrel, "relative error for integration", NULL, 0, 0),
        OPT_DOUBLE('n', "eta", &eta, "set eta", NULL, 0, 0),
        OPT_BOOLEAN('v', "verbose", &verbose, "be verbose", NULL, 0, 0),
        OPT_STRING('D', "detalg", &detalg, "algorithm to compute determinants (HODLR, LU, QR or CHOLESKY)", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse, "\nCompute the Casimir interaction in the cylinder-plane geometry.", NULL);
    argc = argparse_parse(&argparse, argc, argv);

    /* check arguments */
    if(isnan(R))
    {
        fprintf(stderr, "missing mandatory option: -R radius of sphere\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(R <= 0)
    {
        fprintf(stderr, "radius of sphere must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(isnan(d))
    {
        fprintf(stderr, "missing mandatory option: -d separation between cylinder and plate\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(d <= 0)
    {
        fprintf(stderr, "separation between cylinder and plate must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(T != 0)
    {
        fprintf(stderr, "Sorry, but finite temperature is not implemented yet.\n");
        return 1;
    }
    /*
    if(T < 0)
    {
        fprintf(stderr, "temperature must be non-negative\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    */
    if(epsrel <= 0)
    {
        fprintf(stderr, "epsrel must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(eta <= 0)
    {
        fprintf(stderr, "eta must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }

    /* PFA for Dirichlet/Neumann in units of hbar*c*L, i.e., E_PFA^DN / (hbar*c*L) */
    double E_PFA_DN = -M_PI*M_PI*M_PI/1920*sqrt(R/(2*d))/(d*d);
    /* PFA for EM in units of hbar*c*L, i.e., E_PFA^DN / (hbar*c*L) */
    double E_PFA = 2*E_PFA_DN;

    CasimirCP casimir(R, d);

    // set lmax
    if(lmax > 0)
        casimir.set_lmax(lmax);
    else
        casimir.set_lmax(std::max(25,(int)ceil(eta/d*R)));

    printf("# R/d = %.15g\n", R/d);
    printf("# d = %.15g\n", d);
    printf("# R = %.15g\n", R);
    printf("# T = %.15g\n", T);
    printf("# lmax = %d\n", casimir.get_lmax());
    printf("# epsrel = %g\n", epsrel);

    // set detalg
    if(detalg != NULL)
    {
        if(strcasecmp(detalg, "LU") == 0)
        {
            printf("# detalg = LU\n");
            casimir.set_detalg(DETALG_LU);
        }
        else if(strcasecmp(detalg, "QR") == 0)
        {
            printf("# detalg = QR\n");
            casimir.set_detalg(DETALG_QR);
        }
        else if(strcasecmp(detalg, "CHOLESKY") == 0)
        {
            printf("# detalg = CHOLESKY\n");
            casimir.set_detalg(DETALG_CHOLESKY);
        }
        else
            printf("# detalg = HODLR\n");
    }

    /* energy Dirichlet in units of hbar*c*L */
    double E_D = casimir.energy('D', 0, epsrel);
    double E_N = casimir.energy('N', 0, epsrel);

    /* energy EM in units of hbar*c*L */
    double E_EM = E_D+E_N;

    printf("#\n");

    printf("# d/R, d, R, T, lmax, E_PFA/(L*hbar*c), E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA\n");
    printf("%.15g, %.15g, %.15g, %.15g, %d, %.15g, %.15g, %.15g, %.15g\n", d/R, d, R, T, casimir.get_lmax(), E_PFA, E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA);

    return 0;
}

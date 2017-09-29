#ifndef CASIMIR_SCALAR
#define CASIMIR_SCALAR

typedef struct {
    double max;
    double xi_;
    int l,m;
} integrand_t;

typedef struct {
    double alpha;
    double LbyR;
    int ldim;
    char X;
    char Y;
    double epsrel;
} integrand_xi_t;

typedef struct {
    double xi_;    /*< xi_ = xi*(L+R)/c */
    double epsrel; /*< relative error for numerical integration */
    int m;         /*< azimutal quantum number m */
    double *K;
} integration_scalar_t;

typedef struct {
    double LbyR;   /*< L/R */
    int ldim;      /*< dimension of vector space */
    int m;         /*< azimuthal quantum number m */
    double xi_;    /*< xi*R/c */
    double epsrel; /*< relative accuracy for integration */
    char X;        /*< boundary condition on plate (D or N) */
    char Y;        /*< boundary condition on sphere (D or N) */
    double *rs;    /*< cache for rs */
    integration_scalar_t *integration;
    casimir_t *casimir;
} args_t;

integration_scalar_t *integration_init(int m, double xi_, int ldim, double epsrel);
void integration_free(integration_scalar_t *self);

double logdetD_m(double LbyR, double xi_, int m, int ldim, char X, char Y, double epsrel);
double logdetD(double LbyR, double xi_, int ldim, char X, char Y, double epsrel);
double casimir_E(double LbyR, char X, char Y, int ldim, double epsrel);;

#endif

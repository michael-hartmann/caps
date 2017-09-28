#ifndef CASIMIR_SCALAR
#define CASIMIR_SCALAR

typedef struct {
    double max;
    double xi_;
    int l,m;
} integrand_t;

typedef struct {
    double xi_;    /*< xi_ = xi*(L+R)/c */
    double epsrel; /*< relative error for numerical integration */
    int m;         /*< azimutal quantum number m */
    double *K;
} integration_t;

typedef struct {
    double LbyR;   /*< L/R */
    int ldim;      /*< dimension of vector space */
    int m;         /*< azimuthal quantum number m */
    double xi_;    /*< xi*R/c */
    double epsrel; /*< relative accuracy for integration */
    char X;        /*< boundary condition on plate (D or N) */
    char Y;        /*< boundary condition on sphere (D or N) */
    double *rs;    /*< cache for rs */
    integration_t *integration;
} args_t;

double logdetD(double LbyR, double xi_, int m, int ldim, char X, char Y, double epsrel);
integration_t *integration_init(int m, double xi_, int ldim, double epsrel);
void integration_free(integration_t *self);

#endif

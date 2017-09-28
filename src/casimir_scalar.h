#ifndef CASIMIR_SCALAR
#define CASIMIR_SCALAR

typedef struct {
    double LbyR;
    double xi_; /* xi_ = xi*(L+R)/c */
    double epsrel;
    int m;
    char X,Y;
    double *cache_rs;
} args_t;

typedef struct {
    double max;
    double tau;
    int l1,l2,m;
} integrand_t;

double logdetD(double LbyR, double xi_, int m, int ldim, char X, char Y, double epsrel);

#endif

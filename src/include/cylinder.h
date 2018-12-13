#ifndef __CASIMIR_CYLINDER
#define __CASIMIR_CYLINDER

typedef struct {
    int lmax, type;
    char DN;
    double alpha; /* 2/rho */
    double *cache_ratio, *cache_K;
} kernel_args_t;

#endif

#ifndef __CAPC
#define __CAPC

typedef struct {
    int lmax, type;
    char DN;
    double alpha; /* 2/rho */
    double *cache_ratio, *cache_K;
} kernel_args_t;

#endif

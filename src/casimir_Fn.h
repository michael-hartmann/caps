#ifndef CASIMIR_FN
#define CASIMIR_FN

#include <pthread.h>

typedef struct {
    double LbyR, T, precision;
    int lmax, n, m;
    double time, logdet;
} param_t;

double sumF(double *values, int lmax);
void usage(FILE *stream, const char *name);
void *logdetD(void *p);
pthread_t *start_thread(double LbyR, double T, int n, int m, int lmax, double precision);

#endif

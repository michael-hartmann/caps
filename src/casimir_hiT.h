#ifndef CASIMIR_HIT
#define CASIMIR_HIT

#include <pthread.h>

typedef struct {
    double LbyR, precision;
    int lmax, m;
    double value, time, logdet_EE, logdet_MM;
} param_t;

double sumF(double *values, int lmax);
void usage(FILE *stream, const char *name);
void *logdetD0(void *p);
pthread_t *start_thread(double LbyR, int m, int lmax, double precision);

#endif

#ifndef CASIMIR_HIT
#define CASIMIR_HIT

#include <pthread.h>

typedef struct {
    double LbyR, threshold;
    int ldim, m;
    double value, time, logdet_EE, logdet_MM;
} param_t;

double sumF(double *values, int ldim);
void usage(FILE *stream, const char *name);
void *logdetD0(void *p);
pthread_t *start_thread(double LbyR, int m, int ldim, double precision);

#endif

#define _GNU_SOURCE

#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "casimir_hiT.h"
#include "libcasimir.h"
#include "utils.h"

#define IDLE 1000 // in µs
#define LSCALE 7.0
#define PRECISION 5e-9

double sumF(double *values, int ldim)
{
    double F = 0;

    for(int i = ldim-1; i >= 0; i--)
        F += values[i];

    return F;
}

void usage(FILE *stream, const char *name)
{
    fprintf(stream, "Usage: %s -x L/R [-l lscale -L ldim -p precision]\n\n", name);
    fprintf(stream, "Options:\n");
    fprintf(stream, "\t-x L/R:    ratio of L and R, L/R > 0\n");
    fprintf(stream, "\t-L ldim:   set ldim\n");
    fprintf(stream, "\t-l lscale: ldim = lscale*R/L (ignored if -L is used), default: %g\n", LSCALE);
    fprintf(stream, "\t-p prec:   precision, default: %g\n", PRECISION);
    fprintf(stream, "\t-c cores:  cores used for computation, default: 1\n");
}

pthread_t *start_thread(double LbyR, int m, int ldim, double precision)
{
    pthread_t *thread = xmalloc(sizeof(pthread_t));
    param_t *p        = xmalloc(sizeof(param_t));

    p->LbyR      = LbyR;
    p->precision = precision;
    p->ldim      = ldim;
    p->m         = m;
    p->value     = -1;
    p->time      = -1;

    pthread_create(thread, NULL, &logdetD0, (void *)p);

    return thread;
}


void *logdetD0(void *p)
{
    casimir_t *casimir;
    double start = now();
    double logdet_EE = 0, logdet_MM = 0;
    param_t *params = p;
    double LbyR      = params->LbyR;
    double precision = params->precision;
    int m            = params->m;
    int ldim         = params->ldim;

    casimir = casimir_init(LbyR);
    //casimir_set_precision(casimir, precision);
    casimir_set_ldim(casimir, ldim);

    casimir_logdetD0(casimir, m, &logdet_EE, &logdet_MM);
    casimir_free(casimir);

    // n == 0
    logdet_EE /= 2;
    logdet_MM /= 2;

    if(m == 0)
    {
        logdet_EE /= 2;
        logdet_MM /= 2;
    }

    params->value     = logdet_EE+logdet_MM;
    params->logdet_EE = logdet_EE;
    params->logdet_MM = logdet_MM;
    params->time      = now()-start;

    return params;
}

int main(int argc, char *argv[])
{
    double precision = PRECISION;
    double lscale = LSCALE;
    int ldim = -1;
    int cores = 1;
    double start_time = now();
    double LbyR = -1;
    double *values, *EE;
    pthread_t **threads;

    while (1)
    {
        int c;
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },
            { "LbyR",      required_argument, 0, 'x' },
            { "ldim",      required_argument, 0, 'L' },
            { "lscale",    required_argument, 0, 'l' },
            { "precision", required_argument, 0, 'p' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "x:L:l:p:c:h", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch(c)
        {
            case 0:
              /* If this option set a flag, do nothing else now. */
              if (long_options[option_index].flag != 0)
                break;
            case 'x':
                LbyR = atof(optarg);
                break;
            case 'L':
                ldim = atoi(optarg);
                break;
            case 'l':
                lscale = atof(optarg);
                break;
            case 'c':
                cores = atoi(optarg);
                break;
            case 'p':
                precision = atof(optarg);
                break;
            case 'h':
                usage(stdout, argv[0]);
                exit(0);

            case '?':
              /* getopt_long already printed an error message. */
              break;

            default:
              abort();
        }
    }

    if(LbyR <= 0)
    {
        fprintf(stderr, "argument of -x must be nonnegative\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(precision <= 0)
    {
        fprintf(stderr, "argument of -p must be nonnegative\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(cores < 0)
    {
        fprintf(stderr, "argument of -t must be at least 1\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(ldim <= 0)
        ldim = lscale/LbyR;

    // disable buffering
    {
        fflush(stdin);
        fflush(stderr);
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);
    }

    fprintf(stderr, "# L/R=%g, precision=%g, ldim=%d, cores=%d\n", LbyR, precision, ldim, cores);

    values = (double *)xmalloc(ldim*sizeof(double));
    EE     = (double *)xmalloc(ldim*sizeof(double));

    for(int i = 0; i < ldim; i++)
    {
        values[i] = 0;
        EE[i]     = 0;
    }

    threads = (pthread_t **)xmalloc(cores*sizeof(pthread_t));
    for(int i = 0; i < cores; i++)
        threads[i] = NULL;

    int m = 0;
    while(m < ldim)
    {
        int bye = 0;
        param_t *r;

        // try to start threads
        for(int i = 0; i < cores; i++)
        {
            if(threads[i] == NULL)
                threads[i] = start_thread(LbyR, m++, ldim, precision);
            else if(pthread_tryjoin_np(*threads[i], (void *)&r) == 0)
            {
                values[r->m] = r->value;
                EE[r->m]     = r->logdet_EE;
                fprintf(stderr, "# m=%d, value=%.15g, logdet_EE=%.15g, logdet_MM=%.15g, time=%g\n", r->m, r->value, r->logdet_EE, r->logdet_MM, r->time);
                if(r->value/sumF(values, ldim) < precision)
                    bye = 1;
                xfree(r);
                threads[i] = NULL;
            }
        }

        if(bye)
            break;

        usleep(IDLE);
    }

    fprintf(stderr, "# waiting for last threads to finish\n");
    for(int i = 0; i < cores; i++)
    {
        param_t *r;
        if(threads[i] != NULL)
        {
            pthread_join(*threads[i], (void *)&r);
            values[r->m] = r->value;
            EE[r->m]     = r->logdet_EE;
            fprintf(stderr, "# m=%d, value=%.15g, logdet_EE=%.15g, logdet_MM=%.15g, time=%g\n", r->m, r->value, r->logdet_EE, r->logdet_MM, r->time);
            xfree(r);
            threads[i] = NULL;
        }
    }

    xfree(threads);

    printf("# L/R, value_perf, value_Drude, time\n");
    printf("%.15g, %.15g, %.15g, %.15g\n", LbyR, sumF(values, ldim), sumF(EE, ldim), now()-start_time);

    xfree(values);
    xfree(EE);

    return 0;
}

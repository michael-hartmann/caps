#define _GNU_SOURCE

#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "casimir_Fn.h"
#include "libcasimir.h"
#include "integration_perf.h"
#include "utils.h"

#define IDLE 1000 // in µs
#define LSCALE 7.0
#define PRECISION 5e-9
#define VERBOSE 1

double sumF(double *values, int lmax)
{
    int i;
    double F = 0;

    for(i = lmax-1; i >= 0; i--)
        F += values[i];

    return F;
}

void usage(FILE *stream, const char *name)
{
    fprintf(stream, "Usage: %s -x L/R [-l lscale -L lmax -p precision]\n\n", name);
    fprintf(stream, "\t-x L/R:    ratio of L and R, L/R > 0\n");
    fprintf(stream, "\t-n n:      Matsubara term\n");
    fprintf(stream, "\t-T T:      temperature\n");
    fprintf(stream, "\t-L lmax:   use lmax\n"); 
    fprintf(stream, "\t-l lscale: use lmax = lscale*R/L (ignored if -L is used), default: %g\n", LSCALE); 
    fprintf(stream, "\t-p prec:   use precision, default: %g\n", PRECISION); 
    fprintf(stream, "\t-c cores:  how many cores to use, default: 1\n");
}

pthread_t *start_thread(double LbyR, double T, int n, int m, int lmax, double precision)
{
    pthread_t *thread = xmalloc(sizeof(pthread_t));
    param_t *p        = xmalloc(sizeof(param_t));

    p->LbyR      = LbyR;
    p->T         = T;
    p->precision = precision;
    p->lmax      = lmax;
    p->n         = n;
    p->m         = m;
    p->logdet    = 0;
    p->time      = -1;

    pthread_create(thread, NULL, &logdetD, (void *)p);

    return thread;
}


void *logdetD(void *p)
{       
    casimir_t casimir;
    double start = now();
    double logdet;
    param_t *params = p;
    double LbyR      = params->LbyR;
    double T         = params->T;
    double precision = params->precision;
    int n            = params->n;
    int m            = params->m;
    int lmax         = params->lmax;
    integration_perf_t int_perf;

    casimir_integrate_perf_init(&int_perf, n*T, lmax);

    casimir_init(&casimir, 1/(1+LbyR), T);
    casimir_set_precision(&casimir, precision);
    casimir_set_verbose(&casimir, VERBOSE);
    casimir_set_lmax(&casimir, lmax);
    
    logdet = casimir_logdetD(&casimir, n, m, &int_perf);
    casimir_free(&casimir);
    
    params->logdet = logdet;
    params->time   = now()-start;

    return params;
}

int main(int argc, char *argv[])
{
    double precision = PRECISION;
    double lscale = LSCALE;
    int lmax = -1;
    int cores = 1;
    int i, n = -1;
    double start_time = now();
    double T = -1, LbyR = -1;
    double *values;
    pthread_t **threads;

    while (1)
    {
        int c;
        struct option long_options[] =
        {
          { "help",        no_argument,       0, 'h' },
          { "LbyR",        required_argument, 0, 'x' },
          { "temperature", required_argument, 0, 'T' },
          { "matsubara",   required_argument, 0, 'n' },
          { "lmax",        required_argument, 0, 'L' },
          { "lscale",      required_argument, 0, 'l' },
          { "precision",   required_argument, 0, 'p' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "x:L:l:p:n:T:c:h", long_options, &option_index);
      
        /* Detect the end of the options. */
        if (c == -1)
          break;
      
        switch (c)
        {
          case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
              break;
          case 'x':
              LbyR = atof(optarg);
              break;
          case 'L':
              lmax = atoi(optarg);
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
          case 'n':
              n = atoi(optarg);
              break;
          case 'T':
              T = atof(optarg);
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
    if(T <= 0)
    {
        fprintf(stderr, "argument of -T must be nonnegative\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(n < 0)
    {
        fprintf(stderr, "argument of -n must be positive\n\n");
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
    if(lmax <= 0)
        lmax = lscale/LbyR;

    // disable buffering
    {
        fflush(stdin);
        fflush(stderr);
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);
    }

    fprintf(stderr, "# L/R=%.15g, T=%.15g, n=%d, precision=%g, lmax=%d, cores=%d\n", LbyR, T, n, precision, lmax, cores);
    
    values = (double *)xmalloc(lmax*sizeof(double));

    for(i = 0; i < lmax; i++)
        values[i] = 0;

    threads = (pthread_t **)xmalloc(cores*sizeof(pthread_t));
    for(i = 0; i < cores; i++)
        threads[i] = NULL;


    {
        int m = 0;
        int finished = 0;
        int running = 0;

        while(1)
        {
            param_t *r;
            if(m >= lmax)
                finished = 1;

            // try to start threads
            for(i = 0; i < cores; i++)
            {
                if(threads[i] == NULL)
                {
                    if(!finished && m < lmax)
                    {
                        running++;
                        threads[i] = start_thread(LbyR, T, n, m++, lmax, precision);
                    }
                }
                else if(pthread_tryjoin_np(*threads[i], (void *)&r) == 0)
                {
                    running--;

                    values[r->m] = r->logdet;
                    if(n == 0)
                        values[r->m] /= 2;
                    if(r->m == 0)
                        values[r->m] /= 2;
                    fprintf(stderr, "# n=%d, m=%d, value=%.15g, time=%g\n", r->n, r->m, r->logdet, r->time);
                    if(r->logdet/sumF(values, lmax) < precision)
                        finished = 1;
                    xfree(r);
                    threads[i] = NULL;
                }
            }

            if(finished && running == 0)
                break;

            usleep(IDLE);
        }
    }

    xfree(threads);

    printf("# L/R, F_n, time\n");
    printf("%.15g, %.15g, %.15g\n", LbyR, T/M_PI*sumF(values, lmax), now()-start_time);

    xfree(values);

    return 0;
}

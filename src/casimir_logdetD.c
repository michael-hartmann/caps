#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "integration_perf.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

/* print usage */
static void usage(FILE *stream)
{
    fprintf(stream, "Usage: casimir_logdetD [OPTIONS]\n\
This program will calculate the free Casimir energy for the plane-sphere \n\
geometry for given n,m,T,L/R. \n\
\n\
Mandatory options:\n\
    -x L/R\n\
    -T Temperature\n\
    -n value of n\n\
    -m value of m\n\
\n\
Further options:\n\
    -l, --lscale\n\
        Specify parameter lscale. The vector space has to be truncated at some\n\
        value lmax. This program will use lmax=(R/L*lscale) (default: 3)\n\
\n\
    -L\n\
        Set lmax to given value. When -L is used, -l will be ignored\n\
\n\
    --buffering\n\
        Enable buffering. By default buffering for stderr and stdout is\n\
        disabled.\n\
\n\
    -v, --verbose\n\
        Be more verbose.\n\
\n\
    -h,--help\n\
        Show this help\n\
\n\
\n\
Compiled %s, %s\n", __DATE__, __TIME__);
}

int main(int argc, char *argv[])
{
    double T = -1;
    double lfac = 5;
    double LbyR = -1;
    int i, n = -1, m = -1;
    int lmax = 0;
    int buffering_flag = 0, verbose_flag = 0;

    printf("# %s", argv[0]);
    for(i = 1; i < argc; i++)
        printf(", %s", argv[i]);
    printf("\n");

    while (1)
    {
        int c;
        struct option long_options[] =
        {
          { "verbose",   no_argument,       &verbose_flag,   1 },
          { "buffering", no_argument,       &buffering_flag, 1 },
          { "help",      no_argument,       0, 'h' },
          { "lscale",    required_argument, 0, 'l' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "x:T:n:m:s:a:l:L:vqh", long_options, &option_index);
      
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
          case 'T':
              T = atof(optarg);
              break;
          case 'L':
              lmax = atoi(optarg);
          case 'v':
              verbose_flag = 1;
              break;
          case 'l':
              lfac = atof(optarg);
              break;
          case 'n':
              n = atoi(optarg);
              break;
          case 'm':
              m = atoi(optarg);
              break;
          case 'h':
              usage(stdout);
              exit(0);
      
          case '?':
            /* getopt_long already printed an error message. */
            break;
      
          default:
            abort();
        }
    }

    // disable buffering
    if(!buffering_flag)
    {
        fflush(stdin);
        fflush(stderr);
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);
    }

    if(lfac <= 0)
    {
        fprintf(stderr, "--lfac must be positive\n\n");
        usage(stderr);
        exit(1);
    }
    if(LbyR <= 0)
    {
        fprintf(stderr, "-x must be positive\n\n");
        usage(stderr);
        exit(1);
    }
    if(T <= 0)
    {
        fprintf(stderr, "positive value for -T required\n\n");
        usage(stderr);
        exit(1);
    }
    if(n < 0 || m < 0)
    {
        fprintf(stderr, "n,m >= 0\n\n");
        usage(stderr);
        exit(1);
    }

    if(lmax <= 0)
        lmax = MAX((int)ceil(lfac/LbyR), 5);

    printf("# L/R  = %g\n", LbyR);
    printf("# n    = %d\n", n);
    printf("# T    = %g\n", T);
    printf("# xi   = nT = %g\n", n*T);
    printf("# m    = %d\n", m);
    if(lmax)
        printf("# lmax = %d\n", lmax);
    else
        printf("# lfac = %g\n", lfac);
    printf("#\n");

    {
        const double Q = 1./(1.+LbyR);
        integration_perf_t int_perf;
        casimir_t casimir;
        double value, start_time = now();

        casimir_integrate_perf_init(&int_perf, n*T, lmax);

        casimir_init(&casimir, Q, T);
        casimir_set_verbose(&casimir, verbose_flag);
        casimir_set_lmax(&casimir, lmax);

        value = casimir_logdetD(&casimir, n, m, &int_perf);

        casimir_integrate_perf_free(&int_perf);
        casimir_free(&casimir);

        printf("# LbyR,T,n,m,value,lmax,time\n");
        printf("%g, %g, %d, %d, %.15g, %d, %g\n", LbyR, T, n, m, value, casimir.lmax, now()-start_time);
    }

    return 0;
}

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "integration_perf.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

/* print usage */
static void usage(FILE *stream)
{
    char info[1024] = { 0 };
    casimir_compile_info(info, sizeof(info));

    fprintf(stream,
"Usage: casimir_logdetD [OPTIONS]\n"
"This program will calculate the free Casimir energy for the plane-sphere\n"
"geometry for given n,m,T,L/R.\n"
"\n"
"Mandatory options:\n"
"    -x  L/R\n"
"    -nT imaginary frequency ξ\n"
"    -m  value of m\n"
"\n"
"Further options:\n"
"    -w, --omegap\n"
"       Set value of Plasma frequency omega_p of Drude metals in units of\n"
"       c/omegaP. If omitted, omegap = INFINITY.\n"
"\n"
"    -g, --gamma\n"
"       Set value of relaxation frequency gamma of Drude metals in units of\n"
"       c/(L+R). If omitted, gamma = 0.\n"
"\n"
"    -L\n"
"        Set lmax to given value. When -L is used, -l will be ignored.\n"
"\n"
"    -b, --buffering\n"
"        Enable buffering. By default buffering for stderr and stdout is\n"
"        disabled.\n"
"\n"
"    --trace THRESHOLD\n"
"        Try to calculate log det D(xi) using -Tr D(xi). If |trace|>THRESHOLD,\n"
"        fall back to log det D(xi).\n"
"\n"
"    --detalg DETALG\n"
"        Use DETALG to calculate determinant.\n"
"\n"
"    -D, --debug\n"
"        Enable debugging information.\n"
"\n"
"    -h,--help\n"
"        Show this help.\n"
"\n"
"%s\n", info);
}

int main(int argc, char *argv[])
{
    double start_time = now();

    /* geometry, Matsubara frequency */
    double LbyR = -1, nT = -1;
    int m = -1;

    /* material properties, be default: perfect reflectors */
    double gamma_ = 0, omegap = INFINITY;

    /* numerical parameters */
    int lmax = 0;
    double trace_threshold = -1;
    char detalg[64] = { 0 };

    /* flags */
    bool debug = false, buffering = false;

    while(1)
    {
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },
            { "buffering", no_argument,       0, 'b' },
            { "debug",     no_argument,       0, 'D' },

            { "nT",        required_argument, 0, 'T' },
            { "detalg",    required_argument, 0, 'd' },
            { "lscale",    required_argument, 0, 'l' },
            { "trace",     required_argument, 0, 't' },
            { "omegap",    required_argument, 0, 'w' },
            { "gamma",     required_argument, 0, 'g' },

            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long (argc, argv, "x:T:m:s:a:l:w:g:L:t:bqhDv", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch(c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if(long_options[option_index].flag != 0)
                    break;
            case 'x':
                LbyR = atof(optarg);
                break;
            case 'T':
                nT = atof(optarg);
                break;
            case 'w':
                omegap = atof(optarg);
                break;
            case 'g':
                gamma_ = atof(optarg);
                break;
            case 'L':
                lmax = atoi(optarg);
                break;
            case 'm':
                m = atoi(optarg);
                break;
            case 't':
                trace_threshold = atof(optarg);
                break;
            case 'd':
                strncpy(detalg, optarg, sizeof(detalg)/sizeof(char)-1);
                break;
            case 'b':
                buffering = true;
                break;
            case 'D':
                debug = true;
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

    /* disable buffering */
    if(!buffering)
    {
        fflush(stdin);
        fflush(stderr);
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);
    }

    /* check parameters */
    do {
        if(LbyR <= 0)
            fprintf(stderr, "-x must be positive.");
        else if(nT <= 0)
            fprintf(stderr, "--nT must be positive value.");
        else if(m < 0)
            fprintf(stderr, "m >= 0\n\n");
        else if(omegap < 0)
            fprintf(stderr, "--omegap, -w must be non negative.");
        else if(gamma_ < 0)
            fprintf(stderr, "--gamma, -g must be non negative.");
        else
            /* everything ok */
            break;

        /* error occured: print usage and exit */
        fprintf(stderr, "\n\n");
        usage(stderr);
        exit(1);
    } while(0);

    /* print command line options to stdout */
    printf("# %s", argv[0]);
    for(int i = 1; i < argc; i++)
        printf(", %s", argv[i]);
    printf("\n");

    casimir_t casimir;
    casimir_init(&casimir, LbyR, nT);

    if(lmax)
        casimir_set_lmax(&casimir, lmax);

    if(gamma_ > 0)
    {
        casimir_set_gamma_sphere(&casimir, gamma_);
        casimir_set_gamma_plane (&casimir, gamma_);
    }
    if(isfinite(omegap))
    {
        casimir_set_omegap_sphere(&casimir, omegap);
        casimir_set_omegap_plane (&casimir, omegap);
    }

    casimir_set_debug(&casimir, debug);

    if(strlen(detalg))
        casimir_set_detalg(&casimir, detalg);

    if(trace_threshold >= 0)
        casimir_set_trace_threshold(&casimir, trace_threshold);

    casimir_info(&casimir, stdout, "# ");
    printf("#\n");

    double logdet = casimir_logdetD(&casimir, 1, m);

    printf("# LbyR, nT, omegap, gamma, m, logdetD, lmax, time\n");
    printf("%g, %g, %g, %g, %d, %.15g, %d, %g\n", LbyR, nT, omegap, gamma_, m, logdet, casimir.lmax, now()-start_time);

    casimir_free(&casimir);

    return 0;
}

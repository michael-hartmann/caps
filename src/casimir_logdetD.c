#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

/* print usage */
static void usage(FILE *stream)
{
    char info[1024] = { 0 };
    casimir_compile_info(info, sizeof(info));

    fprintf(stream,
"Usage: casimir_logdetD [OPTIONS]\n\n"
"This program will calculate the free Casimir energy for the plane-sphere\n"
"geometry for given n,m,T,L/R.\n"
"\n"
"Mandatory options:\n"
"    -x, --LbyR  L/R\n"
"    --nT        imaginary frequency ξ in units of c/(L+R)\n"
"    -m          value of m\n"
"\n"
"Further options:\n"
"    -w, --omegap OMEGAP\n"
"       Set value of Plasma frequency omega_p of Drude metals in units of\n"
"       c/(L+R). If omitted, omegap = INFINITY.\n"
"\n"
"    -g, --gamma GAMMA\n"
"       Set value of relaxation frequency gamma of Drude metals in units of\n"
"       c/(L+R). If omitted, gamma = 0.\n"
"\n"
"    -L, --lmax LMAX\n"
"        Set lmax to LMAX. When -L is used, -l will be ignored.\n"
"\n"
"    -b, --buffering\n"
"        Enable buffering. By default buffering for stderr and stdout is\n"
"        disabled.\n"
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

    /* flags */
    bool debug = false, buffering = false;

    while(1)
    {
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },
            { "buffering", no_argument,       0, 'b' },
            { "debug",     no_argument,       0, 'D' },

            { "LbyR",      required_argument, 0, 'x' },
            { "nT",        required_argument, 0, 'T' },
            { "lmax",      required_argument, 0, 'L' },
            { "lscale",    required_argument, 0, 'l' },
            { "omegap",    required_argument, 0, 'w' },
            { "gamma",     required_argument, 0, 'g' },

            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long (argc, argv, "x:T:m:l:w:g:L:t:bhD", long_options, &option_index);

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
        else if(nT < 0)
            fprintf(stderr, "--nT must be non-negative value.");
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

    casimir_t *casimir;
    casimir = casimir_init(LbyR);

    if(lmax)
        casimir_set_lmax(casimir, lmax);

    /* XXX
    if(gamma_ >= 0 && isfinite(omegap))
        casimir_set_drude(casimir, omegap, gamma_, omegap, gamma_);
    */

    casimir_set_debug(casimir, debug);

    casimir_info(casimir, stdout, "# ");
    printf("#\n");

    if(nT > 0)
    {
        double logdet = casimir_logdetD(casimir, nT, m);

        printf("# L/R, ξ*(L+R)/c, ωp*(L+R)/c, γ*(L+R)/c, m, logdet(Id-M), lmax, time\n");
        printf("%g, %g, %g, %g, %d, %.16g, %d, %g\n", LbyR, nT, omegap, gamma_, m, logdet, casimir_get_lmax(casimir), now()-start_time);
    }
    else /* nT == 0 */
    {
        /* Calculate logdet(Id-M) for xi=0. In order to reduce the memory
         * footprint, we first calculate logdet(Id-EE) and then logdet(Id-MM).
         * Moreover, we calculate the traces of EE and MM.
         */
        matrix_t *EE, *MM;

        casimir_M0(casimir, m, &EE, NULL);

        double trace_EE = matrix_trace(EE);
        double logdet_EE = matrix_logdet(EE, -1, casimir->detalg);
        matrix_free(EE);

        casimir_M0(casimir, m, NULL, &MM);
        double trace_MM = matrix_trace(MM);
        double logdet_MM = matrix_logdet(MM, -1, casimir->detalg);
        matrix_free(MM);

        printf("# L/R, ξ*(L+R)/c, m, logdet(Id-EE), -trace(EE), logdet(Id-MM), -trace(MM), lmax, time\n");
        printf("%g, 0, %d, %.16g, %.16g, %.16g, %.16g, %d, %g\n", LbyR, m, logdet_EE, -trace_EE, logdet_MM, -trace_MM, casimir_get_lmax(casimir), now()-start_time);
    }

    casimir_free(casimir);

    return 0;
}

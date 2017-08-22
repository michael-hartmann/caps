#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "material.h"
#include "libcasimir.h"
#include "misc.h"
#include "utils.h"


/* print usage */
static void usage(FILE *stream)
{
    fprintf(stream,
"Usage: casimir_logdetD [OPTIONS]\n\n"
"This program will calculate the free Casimir energy for the plane-sphere\n"
"geometry for given n,m,T,L/R.\n"
"\n"
"Mandatory options:\n"
"    -L     separation L\n"
"    -R     radius R\n"
"    --nT   imaginary frequency ξ in units of c/(L+R)\n"
"    -m     value of m\n"
"\n"
"Further options:\n"
"    -l, --ldim LDIM\n"
"        Set ldim to LDIM. When -L is used, -l will be ignored.\n"
"\n"
"    -f, --material FILENAME\n"
"        Use material described by FILENAME.\n"
"\n"
"    -b, --buffering\n"
"        Enable buffering. By default buffering for stderr and stdout is\n"
"        disabled.\n"
"\n"
"    -d, --dense\n"
"        Compute the dense matrix and use LU decomposition to calculate the determinant\n"
"\n"
"    -h,--help\n"
"        Show this help.\n"
"\n");
}

int main(int argc, char *argv[])
{
    double start_time = now();

    char filename[512] = { 0 };

    /* geometry, Matsubara frequency */
    double L = 0, R = 0, nT = -1;
    int m = -1;

    /* numerical parameters */
    int ldim = 0;

    /* flags */
    bool buffering = false, dense = false;

    while(1)
    {
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },
            { "buffering", no_argument,       0, 'b' },
            { "dense",     no_argument,       0, 'd' },

            { "material",  required_argument, 0, 'f' },
            { "nT",        required_argument, 0, 'T' },
            { "ldim",      required_argument, 0, 'l' },

            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long (argc, argv, "L:R:T:m:l:f:t:bdh", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch(c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if(long_options[option_index].flag != 0)
                    break;
            case 'L':
                L = atof(optarg);
                break;
            case 'R':
                R = atof(optarg);
                break;
            case 'T':
                nT = atof(optarg);
                break;
            case 'l':
                ldim = atoi(optarg);
                break;
            case 'm':
                m = atoi(optarg);
                break;
            case 'f':
                strncpy(filename, optarg, sizeof(filename)-sizeof(char));
                break;
            case 'b':
                buffering = true;
                break;
            case 'd':
                dense = true;
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
        if(L <= 0)
            fprintf(stderr, "-L must be positive.");
        if(R <= 0)
            fprintf(stderr, "-R must be positive.");
        else if(nT < 0)
            fprintf(stderr, "--nT must be non-negative value.");
        else if(m < 0)
            fprintf(stderr, "m >= 0\n\n");
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
    casimir = casimir_init(L/R);

    material_t *material = NULL;
    if(strlen(filename) > 0)
    {
        material = material_init(filename, L+R);
        if(material == NULL)
        {
            fprintf(stderr, "Can't read %s or invalid format\n", filename);
            usage(stderr);
            exit(1);
        }

        casimir_set_epsilonm1(casimir, material_epsilonm1, material);
    }

    /* set dimension of vector space */
    if(ldim)
        casimir_set_ldim(casimir, ldim);

    if(dense)
        casimir_set_detalg(casimir, DETALG_LU);

    casimir_info(casimir, stdout, "# ");
    printf("#\n");

    if(nT > 0)
    {
        double logdet = casimir_logdetD(casimir, nT, m);

        printf("# L, R, ξ*(L+R)/c, m, logdet(Id-M), ldim, time\n");
        printf("%g, %g, %g, %d, %.16g, %d, %g\n", L, R, nT, m, logdet, casimir_get_ldim(casimir), now()-start_time);
    }
    else /* nT == 0 */
    {
        double logdet_EE, logdet_MM;
        casimir_logdetD0(casimir, m, 0, &logdet_EE, &logdet_MM, NULL);

        printf("# L, R, ξ*(L+R)/c, m, logdet(Id-EE), logdet(Id-MM), lmax, time\n");
        printf("%g, %g, 0, %d, %.16g, %.16g, %d, %g\n", L, R, m, logdet_EE, logdet_MM, casimir_get_ldim(casimir), now()-start_time);
    }

    casimir_free(casimir);

    if(material != NULL)
        material_free(material);

    return 0;
}

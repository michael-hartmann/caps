#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "material.h"
#include "libcaps.h"
#include "misc.h"
#include "utils.h"


/* print usage */
static void usage(FILE *stream)
{
    fprintf(stream,
"Usage: capc_logdetD [OPTIONS]\n\n"
"This program will calculate the free Casimir energy for the plane-sphere\n"
"geometry for given n,m,T,L/R.\n"
"\n"
"Mandatory options:\n"
"    -L     separation L\n"
"    -R     radius R\n"
"    --xi   imaginary frequency ξ in units of (L+R)/c\n"
"    -m     value of m\n"
"\n"
"Further options:\n"
"    -l, --ldim LDIM\n"
"        Set ldim to LDIM.\n"
"\n"
"    -f, --material FILENAME\n"
"        Use material described by FILENAME.\n"
"\n"
"    -d, --detalg DETALG\n"
"        Compute the matrix using DETALG (LU, QR, CHOLESKY, HODLR)\n"
"\n"
"    -i, --iepsrel IEPSREL\n"
"        Relative accuracy to evaluate integrals\n"
"\n"
"    -h,--help\n"
"        Show this help.\n"
"\n"
"\n"
"Environment variables:\n"
"   CAPS_DUMP:\n"
"        If this variable is set, the round-trip matrix will be dumped in numpy\n"
"        format to the filename contained in CAPS_DUMP. Please note that the\n"
"        round-trip matrix will only be dumped if detalg is QR, LU or CHOLESKY.\n"
"\n"
"   CAPS_CACHE_ELEMS:\n"
"        Determines the size of the cache for the integrals I."
"\n");
}

int main(int argc, char *argv[])
{
    double iepsrel = 0;
    double start_time = now();
    detalg_t detalg = DETALG_HODLR;

    char filename[512] = { 0 };

    /* geometry, Matsubara frequency */
    double L = 0, R = 0, xi_ = -1;
    int m = -1;

    /* numerical parameters */
    int ldim = 0;

    while(1)
    {
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },

            { "iepsrel",   required_argument, 0, 'i' },
            { "detalg",    required_argument, 0, 'd' },
            { "material",  required_argument, 0, 'f' },
            { "xi",        required_argument, 0, 'x' },
            { "ldim",      required_argument, 0, 'l' },

            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "L:R:T:m:l:f:d:i:bh", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch(c)
        {
            case 0:
                break;
            case 'L':
                L = atof(optarg);
                break;
            case 'R':
                R = atof(optarg);
                break;
            case 'x':
                xi_ = atof(optarg);
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
            case 'd':
                if(strcaseequal(optarg, "HODLR"))
                    detalg = DETALG_HODLR;
                else if(strcaseequal(optarg, "LU"))
                    detalg = DETALG_LU;
                else if(strcaseequal(optarg, "QR"))
                    detalg = DETALG_QR;
                else if(strcaseequal(optarg, "CHOLESKY"))
                    detalg = DETALG_CHOLESKY;
                else
                {
                    fprintf(stderr, "Unknown algorithm: %s\n\n", optarg);
                    usage(stderr);
                    exit(1);
                }
                break;
            case 'i':
                iepsrel = atof(optarg);
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
    disable_buffering();

    /* check parameters */
    do {
        if(L <= 0)
            fprintf(stderr, "-L must be positive.");
        if(R <= 0)
            fprintf(stderr, "-R must be positive.");
        else if(xi_ < 0)
            fprintf(stderr, "--xi must be non-negative value.");
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

    caps_t *caps;
    caps = caps_init(R,L);

    if(iepsrel > 0)
        caps_set_epsrel(caps, iepsrel);

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

        caps_set_epsilonm1(caps, material_epsilonm1, material);
    }

    /* set dimension of vector space */
    if(ldim)
        caps_set_ldim(caps, ldim);

    caps_set_detalg(caps, detalg);

    caps_info(caps, stdout, "# ");
    printf("#\n");

    if(xi_ > 0)
    {
        double logdet = caps_logdetD(caps, xi_, m);

        printf("# L, R, ξ*(L+R)/c, m, logdet(Id-M), ldim, time\n");
        printf("%g, %g, %g, %d, %.16g, %d, %g\n", L, R, xi_, m, logdet, caps_get_ldim(caps), now()-start_time);
    }
    else /* xi_ == 0 */
    {
        double logdet_EE, logdet_MM;
        caps_logdetD0(caps, m, 0, &logdet_EE, &logdet_MM, NULL);

        printf("# L, R, ξ*(L+R)/c, m, logdet(Id-EE), logdet(Id-MM), lmax, time\n");
        printf("%g, %g, 0, %d, %.16g, %.16g, %d, %g\n", L, R, m, logdet_EE, logdet_MM, caps_get_ldim(caps), now()-start_time);
    }

    caps_free(caps);

    if(material != NULL)
        material_free(material);

    return 0;
}

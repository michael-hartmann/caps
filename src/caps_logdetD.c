#include <stdio.h>

#include "argparse.h"
#include "material.h"
#include "libcaps.h"

static const char *usage_epilog =
"\n\n"
"Environment variables:\n"
"   CAPS_DUMP:\n"
"        If this variable is set, the round-trip matrix will be dumped in numpy\n"
"        format to the filename contained in CAPS_DUMP. Please note that the\n"
"        round-trip matrix will only be dumped if detalg is QR, LU or CHOLESKY.\n"
"\n"
"   CAPS_CACHE_ELEMS:\n"
"        Determines the size of the cache for the integrals I.";

int main(int argc, char *argv[])
{
    double start_time = now();

    /* geometry, Matsubara frequency, filename */
    double L = 0, R = 0, xi_ = -1;
    int m = -1;
    char *filename = NULL;

    /* numerical parameters */
    int ldim = 0;
    double iepsrel = 0;
    detalg_t detalg = DETALG_HODLR;
    char *detalg_str = NULL;

    const char *const usage[] = {
        "caps_logetD -L separation -R radius --xi frequency -m M [further options]",
        NULL,
    };

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Mandatory options"),
        OPT_DOUBLE('L', "separation", &L, "smallest seperation between plane and spheres in m", NULL, 0, 0),
        OPT_DOUBLE('R', "radius", &R, "radius of the sphere in m", NULL, 0, 0),
        OPT_DOUBLE('x', "xi", &xi_, "imaginary frequency ξ in units of (L+R)/c", NULL, 0, 0),
        OPT_INTEGER('m', NULL, &m, "value of m", NULL, 0, 0),
        OPT_GROUP("Further options"),
        OPT_INTEGER('l', "ldim", &ldim, "dimension of vector space", NULL, 0, 0),
        OPT_STRING('f', "material", &filename, "use material described by file", NULL, 0, 0),
        OPT_STRING('d', "detalg", &detalg_str, "relative error for internal integration", NULL, 0, 0),
        OPT_DOUBLE('i', "iepsrel", &iepsrel, "relative accuracy for evaluation of integrals", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse, "\nCompute log det(Id-M^m(xi)) in the plane-sphere geometry", usage_epilog);
    argc = argparse_parse(&argparse, argc, argv);

    /* disable buffering */
    disable_buffering();

    /* check parameters */
    do {
        if(L <= 0)
            fprintf(stderr, "separation L must be positive");
        else if(R <= 0)
            fprintf(stderr, "radius R must be positive");
        else if(xi_ < 0)
            fprintf(stderr, "xi must be non-negative value");
        else if(m < 0)
            fprintf(stderr, "m must be non-negative integer");
        else if(ldim < 0)
            fprintf(stderr, "ldim must be positive integer");
        else if(detalg_str && (caps_detalg_from_string(detalg_str, &detalg) != 1))
            fprintf(stderr, "invalid value for detalg");
        else
            /* everything ok */
            break;

        /* error occured: print usage and exit */
        fprintf(stderr, "\n\n");
        argparse_usage_stream(&argparse, stderr);
        return 1;
    } while(0);

    caps_t *caps = caps_init(R,L);

    if(iepsrel > 0)
        caps_set_epsrel(caps, iepsrel);

    caps_set_detalg(caps, detalg);

    /* set dimension of vector space */
    if(ldim > 0)
        caps_set_ldim(caps, ldim);


    material_t *material = NULL;
    if(filename != NULL)
    {
        if((material = material_init(filename, L+R)) == NULL)
        {
            fprintf(stderr, "Can't read %s or invalid format\n", filename);
            argparse_usage(&argparse);
            return 1;
        }

        caps_set_epsilonm1(caps, material_epsilonm1, material);
    }

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

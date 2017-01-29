#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unistd.h>

#include "casimir.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

/* default values for precision and lfac */
#define DEFAULT_PRECISION 1e-8
#define DEFAULT_LFAC 5
#define MIN_LMAX 20

/** @brief Find character in string
 *
 * This function counts how many times the character c in string str appears.
 *
 * @param str string
 * @param c character
 * retval how many times c is in str
 */
static int cinstr(const char *str, char c)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c)
            i++;

    return i;
}

/** @brief Find n-th occurence of character in string
 *
 * This function returns a pointer to the n-th occurrence of the character c in
 * the string s. If the character c occures less than n times, NULL is
 * returned.
 *
 * @param str string
 * @param c character
 * @param n occurence
 * @retval NULL if c occures less than n times in str
 * @retval ptr pointer to n-th occurence of c
 */
static const char *indexn(const char *str, char c, int n)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c && ++i == n)
            return str;

    return NULL;
}

/* print usage */
void usage(FILE *stream)
{
    char msg[4096];

    casimir_compile_info(msg, sizeof(msg));

    fprintf(stream, "Usage: casimir [OPTIONS]\n"
"This program will calculate the free Casimir energy F(T,L/R) for the\n"
"plane-sphere geometry for aspect ratio L/R and temperature T. The output is in\n"
"units of ħc/(L+R).\n"
"\n"
"Mandatory options:\n"
"    -x, --LbyR L/R\n"
"        Separation L between surface of sphere and plane divided by radius\n"
"        of sphere, where L/R > 0.\n"
"        If you want to calculate several points, you may pass a start and stop\n"
"        value, and the number of points to be calculated.\n"
"        Examples:\n"
"            $ ./casimir -T 1 -x 0.5,0.9,5\n"
"            This will calculate the free Casimir energy for L/R=0.5,0.6,0.7,0.8\n"
"            and 0,9 (linear scale)\n"
"            $ ./casimir -T 1 -x 0.5,0.9,5,log\n"
"            This will calculate five free energies for L/R=0.5,...,0,9, but using\n"
"            a logarithmic scale.\n"
"\n"
"    -T TEMPERATURE\n"
"        Temperature in units of ħc/(2π*kB*(L+R)). You may use the same\n"
"        syntax like for -x to calculate several points.\n"
"\n"
"Further options:\n"
"    -g, --gamma\n"
"        Set value of relaxation frequency γ of Drude metals in units of\n"
"        c/(L+R). If omitted, γ = 0.\n"
"\n"
"    -w, --omegap\n"
"        Set value of Plasma frequency ωp of Drude metals in units of\n"
"        c/(L+R). If omitted, ωp = INFINITY.\n"
"\n"
"    -l, --lscale\n"
"        Specify parameter lscale. The vector space has to be truncated for\n"
"        some maximum value lmax. This program will use\n"
"        lmax=MAX(R/L*lscale,%d) (default: %d)\n"
"\n"
"    -L, --lmax LMAX\n"
"        Set lmax to LMAX. When -L is specified and positive -l will be ignored.\n"
"\n"
"    -p, --precision\n"
"        Set precision to given value. The value determines when the sum over\n"
"        the Matsubara frequencies and the sum over m is truncated. The sum is\n"
"        truncated when |F_n(m)/F_n(0)| < precision. (default: %g)\n"
"\n"
"    --buffering\n"
"        Enable buffering. By default buffering for stderr and stdout is\n"
"        disabled.\n"
"\n"
"    --verbose\n"
"        Be more verbose.\n"
"\n"
"    -q, --quiet\n"
"        The progress is printed to stderr unless this flag is set.\n"
"\n"
"    -h,--help\n"
"        Show this help\n"
"\n"
"\n"
"%s\n", MIN_LMAX, DEFAULT_LFAC, DEFAULT_PRECISION, msg);
}

/* parse a range given for LbyR or T from the command line.
 * Examples:
 * 1) "value"            => list = { value, value, 1, SCALE_LIN }
 * 2) "start,stop,N"     => list = { start, stop, N, SCALE_LIN }
 * 3) "start,stop,N,log" => list = { start, stop, N, SCALE_LOG }
 */
void parse_range(const char param, const char *_optarg, double list[])
{
    int elems = cinstr(_optarg, ','); /* commas in _optarg */
    list[3] = SCALE_LIN;

    switch(elems)
    {
        case 0:
            /* no comma => example 1) */
            list[0] = list[1] = atof(_optarg);
            list[2] = 1;
            break;
        case 3:
            /* 3 commas => example 3) */
            if(strncasecmp(indexn(_optarg, ',', 3)+1, "log", 3) == 0)
                list[3] = SCALE_LOG;
            /* here no break! */
        case 2:
            /* 2 commas => example 2) */
            list[0] = atof(_optarg);
            list[1] = atof(indexn(_optarg, ',', 1)+1);
            list[2] = atoi(indexn(_optarg, ',', 2)+1);

            /* N must be positive */
            if(list[2] <= 0)
            {
                fprintf(stderr, "error parsing parameter -%c\n\n", param);
                usage(stderr);
                exit(1);
            }

            if(list[0] == list[1])
            {
                fprintf(stderr, "start value must not be stop value");
                usage(stderr);
                exit(1);
            }

            /* ensure that start < stop */
            if(list[0] > list[1])
            {
                /* swap list[0] and list[1] */
                double temp = list[0];
                list[0] = list[1];
                list[1] = temp;
            }
            break;

        default:
            fprintf(stderr, "Can't parse range %s.\n\n", _optarg);
            usage(stderr);
            exit(1);
    }
}


static double linspace(double start, double stop, int N, int i)
{
    if(N == 1)
        return start;
    else
        return start+(stop-start)*i/(N-1);
}

static double logspace(double start, double stop, int N, int i)
{
    if(N == 1)
        return start;
    else
        return start*pow(pow(stop/start, 1./(N-1)), i);
}

int main(int argc, char *argv[])
{
    double gamma_ = 0, omegap = INFINITY;
    double precision = DEFAULT_PRECISION;
    double lfac      = DEFAULT_LFAC;
    double lT[4]     = { 0,0,0,SCALE_LIN }; /* start, stop, N, lin/log */
    double lLbyR[4]  = { 0,0,0,SCALE_LIN }; /* start, stop, N, lin/log */
    int lmax = 0;
    int buffering_flag = 0, quiet_flag = 0, verbose_flag = 0;

    /* parse command line options */
    while (1)
    {
        struct option long_options[] =
        {
            { "quiet",     no_argument,       &quiet_flag,     1 },
            { "verbose",   no_argument,       &verbose_flag,   1 },
            { "buffering", no_argument,       &buffering_flag, 1 },
            { "help",      no_argument,       0, 'h' },
            { "LbyR",      required_argument, 0, 'x' },
            { "lscale",    required_argument, 0, 'l' },
            { "lmax",      required_argument, 0, 'L' },
            { "precision", required_argument, 0, 'p' },
            { "gamma",     required_argument, 0, 'g' },
            { "omegap",    required_argument, 0, 'w' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long (argc, argv, "x:T:l:L:p:g:w:qh", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option sets a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
            case 'x':
                parse_range('x', optarg, lLbyR);
                break;
            case 'T':
                parse_range('T', optarg, lT);
                break;
            case 'L':
                lmax = atoi(optarg);
                break;
            case 'q':
                quiet_flag = 1;
                break;
            case 'l':
                lfac = atof(optarg);
                break;
            case 'p':
                precision = atof(optarg);
                break;
            case 'g':
                gamma_ = atof(optarg);
                break;
            case 'w':
                omegap = atof(optarg);
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

    /* check command line arguments */
    do
    {
        if(lfac <= 0)
            fprintf(stderr, "wrong argument for -l, --lscale: lscale must be positive\n\n");
        else if(precision <= 0)
            fprintf(stderr, "wrong argument for -p, --precision: precision must be positive\n\n");
        else if(lLbyR[0] <= 0 || lLbyR[1] <= 0)
            fprintf(stderr, "wrong argument for -x: x=L/R must be positive\n\n");
        else if(lT[0] <= 0 || lT[1] <= 0)
            fprintf(stderr, "wrong argument for -T: temperature must be positive\n\n");
        else if(gamma_ < 0)
            fprintf(stderr, "wrong argument for --gamma: gamma must be nonnegative\n\n");
        else if(omegap <= 0)
            fprintf(stderr, "wrong argument for --omegap: omegap must be positive\n\n");
        else
            /* everythink ok */
            break;

        /* print usage and exit */
        usage(stderr);
        exit(1);
    } while(0);

    if(!quiet_flag)
    {
        char msg[4096];
        casimir_compile_info(msg, sizeof(msg));
        printf("# %s\n#\n", msg);
    }

    int i = 0;
    for(int iLbyR = 0; iLbyR < lLbyR[2]; iLbyR++)
        for(int iT = 0; iT < lT[2]; iT++)
        {
            double start_time = now();
            int nmax;
            double F,LbyR,T;

            if(lLbyR[3] == SCALE_LIN)
                LbyR = linspace(lLbyR[0], lLbyR[1], lLbyR[2], iLbyR);
            else
                LbyR = logspace(lLbyR[0], lLbyR[1], lLbyR[2], iLbyR);

            if(lT[3] == SCALE_LIN)
                T = linspace(lT[0], lT[1], lT[2], iT);
            else
                T = logspace(lT[0], lT[1], lT[2], iT);

            casimir_t *casimir = casimir_init(LbyR, T);

            casimir_set_verbose(casimir, verbose_flag);
            casimir_set_precision(casimir, precision);

            /* XXX
            if(isfinite(omegap))
                casimir_set_drude(&casimir, omegap, gamma_, omegap, gamma_);
            */

            if(lmax > 0)
                casimir_set_lmax(casimir, lmax);
            else
                casimir_set_lmax(casimir, MAX((int)ceil(lfac/LbyR), MIN_LMAX));


            if(!quiet_flag)
                casimir_info(casimir, stdout, "# ");

            F = casimir_F(casimir, &nmax);
            casimir_free(casimir);

            if(!quiet_flag)
                printf("#\n");

            /* if quiet, print this line just once */
            if(i == 0 || !quiet_flag)
                printf("# L/R, 2π*kB*T*(L+R)/(ħc), ωp*(L+R)/c, γ*(L+R)/c, F*(L+R)/(ħc), lmax, nmax, time\n");


            printf("%g, %g, %g, %g, %.12g, %d, %d, %g\n",
                LbyR, T,
                omegap, gamma_,
                F, casimir_get_lmax(casimir), nmax, now()-start_time
            );

            if(!quiet_flag)
            {
                double progress = ++i*100/(lLbyR[2]*lT[2]);
                fprintf(stderr, "# %6.2f%%, L/R=%g, T=%g\n", progress, LbyR, T);
                if(progress != 100)
                    printf("#\n#\n");
            }
        }

    return 0;
}

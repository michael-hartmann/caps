#define _DEFAULT_SOURCE /* make usleep work */

#include <ctype.h>
#include <getopt.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "casimir.h"

#include "quadpack.h"
#include "fcqs.h"
#include "libcasimir.h"
#include "material.h"
#include "misc.h"
#include "utils.h"

#define EPSREL 1e-6
#define CUTOFF 1e-9
#define ETA 7.
#define IDLE 25

#define STATE_RUNNING 1
#define STATE_IDLE    0

/* @brief Create casimir_mpi object
 *
 * @param [in] L separation between sphere and plate in meter
 * @param [in] R radius of sphere in meter
 * @param [in] T temperature in Kelvin
 * @param [in] filename filename of material description or NULL
 * @param [in] omegap plasma frequency of the Drude model in eV
 * @param [in] gamma_ relaxation frequency of the Drude model eV
 * @param [in] ldim dimension of vector space
 * @param [in] cutoff cutoff for summation over m
 * @param [in] iepsrel relative accuracy for integration of k for matrix elements
 * @param [in] cores number of cores to use
 * @param [in] verbose flag if verbose
 * @retval object casimir_mpi_t object
 */
casimir_mpi_t *casimir_mpi_init(double L, double R, double T, char *filename, double omegap, double gamma_, int ldim, double cutoff, double iepsrel, int cores, bool verbose)
{
    casimir_mpi_t *self = xmalloc(sizeof(casimir_mpi_t));

    self->L       = L;
    self->R       = R;
    self->T       = T;
    self->omegap  = omegap;
    self->gamma   = gamma_;
    self->ldim    = ldim;
    self->cutoff  = cutoff;
    self->iepsrel = iepsrel;
    self->cores   = cores;
    self->verbose = verbose;
    self->tasks   = xmalloc(cores*sizeof(casimir_task_t *));
    self->alpha   = 2*L/(L+R);

    /* number of determinants we have computed */
    self->determinants = 0;

    TERMINATE(strlen(filename) > 511, "filename too long: %s", filename);

    memset(self->filename, '\0', sizeof(self->filename));
    strncpy(self->filename, filename, sizeof(self->filename)-sizeof(char));

    self->tasks[0] = NULL;
    for(int i = 1; i < cores; i++)
    {
        casimir_task_t *task = xmalloc(sizeof(casimir_task_t));
        task->index    = -1;
        task->state    = STATE_IDLE;
        self->tasks[i] = task;
    }

    return self;
}

static void _mpi_stop(int cores)
{
    double buf[] = { -1, -1, -1, -1, -1, -1, -1 };

    /* stop all remaining slaves */
    for(int i = 1; i < cores; i++)
        MPI_Send(buf, 7, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
}

void casimir_mpi_free(casimir_mpi_t *self)
{
    _mpi_stop(self->cores);

    for(int i = 1; i < self->cores; i++)
    {
        xfree(self->tasks[i]);
        self->tasks[i] = NULL;
    }

    xfree(self->tasks);
    self->tasks = NULL;

    xfree(self);
}


int casimir_mpi_get_running(casimir_mpi_t *self)
{
    int running = 0;

    for(int i = 1; i < self->cores; i++)
        if(self->tasks[i]->state == STATE_RUNNING)
            running++;

    return running;
}

int casimir_mpi_submit(casimir_mpi_t *self, int index, double xi, int m)
{
    for(int i = 1; i < self->cores; i++)
    {
        casimir_task_t *task = self->tasks[i];

        if(task->state == STATE_IDLE)
        {
            double buf[] = { xi, self->L, self->R, self->omegap, self->gamma, m, self->iepsrel, self->ldim };

            task->index = index;
            task->xi    = xi;
            task->m     = m;
            task->state = STATE_RUNNING;

            MPI_Send (buf,              8, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send (self->filename, 512, MPI_CHAR,   i, 0, MPI_COMM_WORLD);
            MPI_Irecv(&task->recv,     1,  MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &task->request);

            return 1;
        }
    }

    return 0;
}

int casimir_get_determinants(casimir_mpi_t *self)
{
    return self->determinants;
}

int casimir_mpi_retrieve(casimir_mpi_t *self, casimir_task_t **task_out)
{
    *task_out = NULL;

    for(int i = 1; i < self->cores; i++)
    {
        casimir_task_t *task = self->tasks[i];

        if(task->state == STATE_RUNNING)
        {
            int flag = 0;
            MPI_Status status;

            MPI_Test(&task->request, &flag, &status);

            if(flag)
            {
                MPI_Wait(&task->request, &status);

                task->value = task->recv;
                task->state = STATE_IDLE;
                self->determinants += 1;

                *task_out = self->tasks[i];

                return 1;
            }
        }
    }

    return 0;
}

static double wrapper_integrand(double xi_, void *args)
{
    const double t0 = now();
    casimir_mpi_t *casimir_mpi = (casimir_mpi_t *)args;
    const double xi = xi_/casimir_mpi->alpha;
    const double v = integrand(xi, casimir_mpi);
    printf("# xi=%.12g, logdetD=%.15g, t=%g\n", xi, v, now()-t0);
    return v;
}

double integrand(double xi, casimir_mpi_t *casimir_mpi)
{
    int m;
    double terms[1024] = { NAN };
    bool verbose = casimir_mpi->verbose;
    const double mmax = sizeof(terms)/sizeof(double);
    const double cutoff = casimir_mpi->cutoff;

    /* gather all data */
    for(m = 0; m < mmax; m++)
    {
        while(1)
        {
            casimir_task_t *task = NULL;

            /* send job */
            if(casimir_mpi_submit(casimir_mpi, m, xi, m))
                break;

            /* retrieve jobs */
            while(casimir_mpi_retrieve(casimir_mpi, &task))
            {
                double v = terms[task->m] = task->value;

                if(verbose)
                    fprintf(stderr, "# m=%d, xi=%.12g, logdetD=%.12g\n", task->m, xi, task->value);

                if(v == 0 || v/terms[0] < cutoff)
                    goto done;
            }

            usleep(IDLE);
        }
    }

    TERMINATE(true, "sum did not converge, sorry. :(");

    done:

    /* retrieve all remaining running jobs */
    while(casimir_mpi_get_running(casimir_mpi) > 0)
    {
        casimir_task_t *task = NULL;

        while(casimir_mpi_retrieve(casimir_mpi, &task))
        {
            terms[task->m] = task->value;
            if(verbose)
                fprintf(stderr, "# m=%d, xi=%.12g, logdetD=%.12g\n", task->m, xi, task->value);
        }

        usleep(IDLE);
    }

    terms[0] /= 2; /* m = 0 */
    return kahan_sum(terms, m);
}

int main(int argc, char *argv[])
{
    int cores, rank;
    MPI_Comm new_comm;

    /* initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, rank == 0, 0, &new_comm);
    MPI_Comm_size(MPI_COMM_WORLD, &cores);

    if(rank == 0)
        master(argc, argv, cores);
    else
        slave(MPI_COMM_WORLD, rank);

    MPI_Finalize();

    return 0;
}

void master(int argc, char *argv[], const int cores)
{
    bool verbose = false, fcqs = false;
    char filename[512] = { 0 };
    int ldim = 0;
    double L = 0, R = 0, T = 0, omegap = INFINITY, gamma_ = 0;
    double cutoff = CUTOFF, epsrel = EPSREL, eta = ETA;
    double iepsrel = CASIMIR_EPSREL;
    material_t *material = NULL;
    char time_str[128];
    time_t rawtime;
    struct tm *info;

    #define EXIT() do { _mpi_stop(cores); return; } while(0)

    /* parse command line options */
    while (1)
    {
        struct option long_options[] = {
            { "help",        no_argument,       0, 'h' },
            { "verbose",     no_argument,       0, 'v' },
            { "fcqs",        no_argument,       0, 'F' },
            { "temperature", required_argument, 0, 'T' },
            { "eta",         required_argument, 0, 'E' },
            { "ldim",        required_argument, 0, 'l' },
            { "cutoff",      required_argument, 0, 'c' },
            { "epsrel",      required_argument, 0, 'e' },
            { "iepsrel",     required_argument, 0, 'i' },
            { "material",    required_argument, 0, 'f' },
            { "omegap",      required_argument, 0, 'w' },
            { "gamma",       required_argument, 0, 'g' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long(argc, argv, "R:L:T:l:c:e:E:f:i:w:g:FvVh", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option sets a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
            case 'L':
                L = atof(optarg);
                break;
            case 'R':
                R = atof(optarg);
                break;
            case 'T':
                T = atof(optarg);
                break;
            case 'l':
                ldim = atoi(optarg);
                break;
            case 'E':
                eta = atof(optarg);
                break;
            case 'c':
                cutoff = atof(optarg);
                break;
            case 'i':
                iepsrel = atof(optarg);
                break;
            case 'e':
                epsrel = atof(optarg);
                break;
            case 'w':
                omegap = atof(optarg);
                break;
            case 'g':
                gamma_ = atof(optarg);
                break;
            case 'F':
                fcqs = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'f':
                strncpy(filename, optarg, sizeof(filename)-sizeof(char));
                break;
            case 'V':
                casimir_build(stdout, NULL);
                exit(0);
            case 'h':
                usage(stdout);
                EXIT();

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    if(L <= 0)
    {
        fprintf(stderr, "separation must be positive.\n\n");
        usage(stderr);
        EXIT();
    }
    if(R <= 0)
    {
        fprintf(stderr, "radius of sphere must be positive.\n\n");
        usage(stderr);
        EXIT();
    }
    if(T < 0)
    {
        fprintf(stderr, "temperature must be non-negative.\n\n");
        usage(stderr);
        EXIT();
    }
    if(cutoff <= 0)
    {
        fprintf(stderr, "cutoff must be positive.\n\n");
        usage(stderr);
        EXIT();
    }
    if(epsrel <= 0)
    {
        fprintf(stderr, "epsrel must be positive.\n\n");
        usage(stderr);
        EXIT();
    }
    if(iepsrel < 0)
    {
        fprintf(stderr, "iepsrel must be non-negative.\n\n");
        usage(stderr);
        EXIT();
    }
    if(eta <= 0)
    {
        fprintf(stderr, "eta must be positive.\n\n");
        usage(stderr);
        EXIT();
    }
    if(strlen(filename))
    {
        if(access(filename, R_OK) != 0)
        {
            fprintf(stderr, "file %s does not exist or is not readable.\n\n", filename);
            usage(stderr);
            EXIT();
        }
        material = material_init(filename, L+R);
        if(material == NULL)
        {
            fprintf(stderr, "file %s has wrong format.\n\n", filename);
            usage(stderr);
            EXIT();
        }
    }
    if(!isinf(omegap))
    {
        if(omegap <= 0)
        {
            fprintf(stderr, "omegap must be positive\n\n");
            usage(stderr);
            EXIT();
        }
        if(gamma_ < 0)
        {
            fprintf(stderr, "gamma must be non-negative\n\n");
            usage(stderr);
            EXIT();
        }
    }

    if(cores < 2)
    {
        fprintf(stderr, "need at least 2 cores\n");
        EXIT();
    }

    const double LbyR = L/R;

    /* if ldim was not set */
    if(ldim <= 0)
        ldim = ceil(eta/LbyR);

    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(time_str, sizeof(time_str),"%c", info);

    casimir_build(stdout, "# ");
    printf("# pid: %d\n", (int)getpid());
    printf("# start time: %s\n", time_str);
    printf("#\n");

    casimir_mpi_t *casimir_mpi = casimir_mpi_init(L, R, T, filename, omegap, gamma_, ldim, cutoff, iepsrel, cores, verbose);

    /* estimate, cf. eq. (6.33) */
    const double alpha = 2*LbyR/(1+LbyR);

    printf("# LbyR = %.15g\n", LbyR);
    printf("# RbyL = %.15g\n", 1/LbyR);
    printf("# L = %.15g\n", L);
    printf("# R = %.15g\n", R);
    printf("# T = %.15g\n", T);
    printf("# cutoff = %g\n", cutoff);
    printf("# epsrel = %g\n", epsrel);
    printf("# iepsrel = %g\n", iepsrel);
    printf("# ldim = %d\n", ldim);
    printf("# cores = %d\n", cores);
    printf("# alpha = %.15g\n", alpha);
    if(strlen(filename))
        printf("# filename = %s\n", filename);
    else if(!isinf(omegap))
    {
        printf("# omegap = %.15g\n", omegap);
        printf("# gamma = %.15g\n", gamma_);
    }


    double F = NAN;
    if(T == 0)
    {
        double integral = 0;
        double abserr = 0;
        int ier, neval;

        if(fcqs)
            printf("# quad = Fourier-Chebsheb quadrature scheme\n");
        else
            printf("# quad = adaptive Gauss-Kronrod\n");
        printf("#\n");

        if(fcqs)
            integral = fcqs_semiinf(wrapper_integrand, casimir_mpi, &epsrel, &neval, 1, &ier);
        else
        {
            integral = dqagi(wrapper_integrand, 0, 1, 0, epsrel, &abserr, &neval, &ier, casimir_mpi);
            epsrel = fabs(abserr/integral);
        }

        printf("#\n");
        printf("# ier=%d, integral=%.15g, neval=%d, epsrel=%g\n", ier, integral, neval, epsrel);

        WARN(ier != 0, "ier=%d", ier);

        /* free energy for T=0 */
        F = integral/alpha/M_PI;
    }
    else
    {
        const double T_scaled = 2*M_PI*CASIMIR_kB*(R+L)*T/(CASIMIR_hbar*CASIMIR_c);
        /* finite temperature */
        double v[2048] = { 0 };

        /* xi = 0 */
        {
            const double t0 = now();

            casimir_t *casimir = casimir_init(LbyR);
            casimir_set_ldim(casimir, ldim);

            if(material == NULL && isinf(omegap) && gamma_ == 0)
            {
                printf("# model = perfect conductors\n");
                v[0] = casimir_logdetD0_perf(casimir, cutoff);
            }
            else if(material == NULL && gamma_ == 0)
            {
                printf("# model = plasma\n");
                v[0] = casimir_logdetD0_plasma(casimir, omegap, cutoff);
            }
            else if(material == NULL)
            {
                printf("# model = drude\n");
                v[0] = casimir_logdetD0_drude(casimir);
            }
            else
            {
                double omegap_low, gamma_low;
                material_get_extrapolation(material, &omegap_low, &gamma_low, NULL, NULL);

                double omegap_scaled = omegap_low*(L+R)/CASIMIR_c;

                if(gamma_low == 0)
                {
                    printf("# model = optical data (xi=0: Plasma)\n");
                    v[0] = casimir_logdetD0_plasma(casimir, omegap_scaled, cutoff);
                }
                else
                {
                    printf("# model = optical data (xi=0: Drude)\n");
                    printf("# plasma = %.15g (logdetD(xi=0) for plasma model with omegap=%geV)\n", casimir_logdetD0_plasma(casimir, omegap_scaled, cutoff), omegap_low);

                    v[0] = casimir_logdetD0_drude(casimir);
                }
            }

            casimir_free(casimir);

            printf("#\n");
            printf("# xi=0, logdetD=%.15g, t=%g\n", v[0], now()-t0);
        }


        for(size_t n = 1; n < sizeof(v)/sizeof(double); n++)
        {
            const double t0 = now();
            const double xi = n*T_scaled;
            v[n] = integrand(xi, casimir_mpi);
            printf("# xi=%.15g, logdetD=%.15g, t=%g\n", xi, v[n], now()-t0);

            if(fabs(v[n]/v[0]) < epsrel)
            {
                v[0] /= 2;
                F = T_scaled/M_PI*kahan_sum(v, n+1);
                break;
            }
        }
    }

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(time_str, sizeof(time_str), "%c", info);

    printf("#\n");
    printf("# %d determinants computed\n", casimir_get_determinants(casimir_mpi));
    printf("# stop time: %s\n", time_str);
    printf("#\n");
    printf("# L/R, L, R, T, ldim, F*(L+R)/(ħc)\n");
    printf("%.16g, %.16g, %.16g, %.16g, %d, %.16g\n", LbyR, L, R, T, ldim, F);

    if(material != NULL)
        material_free(material);

    casimir_mpi_free(casimir_mpi);
}

void slave(MPI_Comm master_comm, int rank)
{
    char filename[512] = { 0 };
    double userdata[2] = { 0 };
    double buf[8] = { 0 };

    MPI_Status status;
    MPI_Request request;

    while(1)
    {
        casimir_t  *casimir  = NULL;
        material_t *material = NULL;

        memset(buf,      0, sizeof(buf));
        memset(filename, 0, sizeof(filename));
        memset(userdata, 0, sizeof(userdata));

        MPI_Recv(buf, 8, MPI_DOUBLE, 0, 0, master_comm, &status);

        /* Matsubara frequency */
        const double xi = buf[0];

        /* signal to quit */
        if(xi < 0)
            break;

        /* geometry */
        const double L = buf[1]; /* in m */;
        const double R = buf[2]; /* in m */
        const double LbyR = L/R;

        /* material properties (scaled) */
        const double omegap = buf[3]/(CASIMIR_hbar_eV*CASIMIR_c)*(L+R);
        const double gamma_ = buf[4]/(CASIMIR_hbar_eV*CASIMIR_c)*(L+R);

        const int m          = (int)buf[5];
        const double iepsrel = buf[6];
        const int ldim       = (int)buf[7];

        /* get filename */
        MPI_Recv(filename, 512, MPI_CHAR, 0, 0, master_comm, &status);

        casimir = casimir_init(LbyR);
        TERMINATE(casimir == NULL, "casimir object is null");
        casimir_set_ldim(casimir, ldim);

        if(iepsrel > 0)
            casimir_set_epsrel(casimir, iepsrel);

        /* set material properties */
        if(strlen(filename))
        {
            material = material_init(filename, L+R);
            TERMINATE(material == NULL, "material_init failed");
            casimir_set_epsilonm1(casimir, material_epsilonm1, material);
        }
        else if(!isinf(omegap))
        {
            userdata[0] = omegap;
            userdata[1] = gamma_;
            casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);
        }

        double logdet = casimir_logdetD(casimir, xi, m);
        TERMINATE(isnan(logdet), "L/R=%.10g, xi=%.10g, m=%d, ldim=%d", LbyR, xi, m, ldim);
        casimir_free(casimir);

        if(material != NULL)
            material_free(material);

        MPI_Isend(&logdet, 1, MPI_DOUBLE, 0, 0, master_comm, &request);
        MPI_Wait(&request, &status);
    }
}

void usage(FILE *stream)
{
    fprintf(stream,
"Usage: casimir [OPTIONS]\n\n"
"This program computes the free Casimir energy F(T,L/R) for the plane-sphere\n"
"geometry for L, R and T. The output is in units of ħc/(L+R).\n"
"\n"
"This program uses MPI for parallization and needs at leas to cores to run.\n"
"\n"
"The free eenergy at T=0 is calculated using integration:\n"
"   F(L/R) = \\int_0^\\infty dξ logdet D(ξ)\n"
"The integration is performed either using an adaptive Gauß-Kronrod\n"
"quadrature or a Fourier-Chebshev quadrature scheme. The integrand decays\n"
"exponentially, c.f. Eq. (6.33) of [1].\n"
"\n"
"References:\n"
"  [1] Hartmann, Negative Casimir entropies in the plane-sphere geometry, 2014\n"
"\n"
"Mandatory options:\n"
"    -L L\n"
"        Separation L between sphere and plane, L>0 (in m).\n"
"\n"
"    -R R\n"
"        Radius R of the sphere, R>0 (in m).\n"
"\n"
"Further options:\n"
"    -T, --temperature TEMPERATURE\n"
"        Set temperature to TEMPERATURE (default: 0; in K)\n"
"\n"
"    -l, --ldim LDIM\n"
"        Set ldim to the value LDIM. (default: ldim=ceil(eta*R/L))\n"
"\n"
"    --eta ETA\n"
"        Set eta to the value ETA. eta is used to determine ldim if not set\n"
"        by --ldim. (default: eta=%g)\n"
"\n"
"    -c, --cutoff CUTOFF\n"
"        Stop summation over m for a given value of ξ if\n"
"            logdet(D^m(ξ))/logdet(D^0(ξ) < CUTOFF\n"
"        (default: %g)\n"
"\n"
"    -e, --epsrel EPSREL\n"
"       Request relative accuracy of EPSREL for integration over xi if T=0,\n"
"       or stop criterion logdetD(n)/logdetD(n=0) < EPSREL for T>0\n"
"       (default: %g)\n"
"\n"
"    -i, --iepsrel IEPSREL\n"
"       Set relative accuracy of integration over k for the matrix elements\n"
"       to IEPSREL (default: %g)\n"
"\n"
"    -F, --fcqs\n"
"      Use Fourier-Chebshev quadrature scheme to compute integral over xi. This\n"
"      is usually faster than using Gauss-Kronrod. (only used when T=0)\n"
"\n"
"    -f, --material FILENAME\n"
"        Filename of the material description file. If set, --omegap and\n"
"        --gamma will be ignored.\n"
"\n"
"    --omegap OMEGAP\n"
"        Model the metals using the Drude/Plasma model and set Plasma\n"
"        frequency to OMEGAP. (DEFAULT: perfect conductors; in eV)\n"
"\n"
"    --gamma GAMMA\n"
"        Set dissipation of Drude model to GAMMA. Ignored if no value for\n"
"        --omegap is given. (DEFAULT: perfect conductors; in eV)\n"
"\n"
"    -v, --verbose\n"
"        Also print results for each m\n"
"\n"
"    -V\n"
"        Print information about build to stdout and exit\n"
"\n"
"    -h, --help\n"
"        Show this help\n",
    ETA, CUTOFF, EPSREL, CASIMIR_EPSREL);
}

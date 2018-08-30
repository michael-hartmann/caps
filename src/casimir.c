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
#include "psd.h"
#include "utils.h"

#define EPSREL 1e-6
#define CUTOFF 1e-9
#define LDIM_MIN 20
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
    self->alpha   = 2*L/(L+R); /* used to scale integration if T=0 */

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

static double integrand(double xi_, void *args)
{
    const double t0 = now();
    casimir_mpi_t *casimir_mpi = (casimir_mpi_t *)args;
    const double xi = xi_/casimir_mpi->alpha;
    const double v = F_xi(xi, casimir_mpi);
    printf("# xi=%.15g, logdetD=%.15g, t=%g\n", xi, v, now()-t0);
    return v;
}

void F_HT(casimir_mpi_t *casimir_mpi, double omegap, double *drude, double *pr, double *plasma)
{
    const double omegap_orig = casimir_mpi->omegap;
    const double gamma_orig  = casimir_mpi->gamma;

    /* Drude
     * The actual value of omegap and gamma for the high-temperature limit are
     * irrelevant. What is relevant is that gamma is positive and non-zero.
     */
    if(drude != NULL)
    {
        casimir_mpi->omegap = 1;
        casimir_mpi->gamma  = 1;
        *drude = F_xi(0, casimir_mpi);
    }

    /* PR */
    if(pr != NULL)
    {
        casimir_mpi->omegap = INFINITY;
        casimir_mpi->gamma  = 0;
        *pr = F_xi(0, casimir_mpi);
    }

    /* plasma */
    if(plasma != NULL)
    {
        casimir_mpi->omegap = omegap;
        casimir_mpi->gamma  = 0;
        *plasma = F_xi(0, casimir_mpi);
    }

    casimir_mpi->omegap = omegap_orig;
    casimir_mpi->gamma  = gamma_orig;
}

double F_xi(double xi, casimir_mpi_t *casimir_mpi)
{
    int m;
    double drude_HT = NAN;
    double terms[4096] = { NAN };
    bool verbose = casimir_mpi->verbose;
    const double mmax = sizeof(terms)/sizeof(double);
    const double cutoff = casimir_mpi->cutoff;

    if(xi == 0)
    {
        /* compute Drude contribution */
        casimir_t *casimir = casimir_init(casimir_mpi->L/casimir_mpi->R);
        casimir_set_ldim(casimir, casimir_mpi->ldim);
        if(casimir_mpi->iepsrel > 0)
            casimir_set_epsrel(casimir, casimir_mpi->iepsrel);
        drude_HT = casimir_logdetD0_drude(casimir);
        casimir_free(casimir);

        if(!isinf(casimir_mpi->omegap) && casimir_mpi->gamma > 0)
            /* omegap finite => Drude model */
            return drude_HT;
    }

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

    if(xi == 0)
        return drude_HT + kahan_sum(terms, m);
    else
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
    bool verbose = false, fcqs = false, ht = false;
    char filename[512] = { 0 };
    int ldim = 0;
    double L = 0, R = 0, T = 0, omegap = INFINITY, gamma_ = 0;
    double cutoff = CUTOFF, epsrel = EPSREL, eta = ETA;
    double iepsrel = CASIMIR_EPSREL;
    material_t *material = NULL;
    char time_str[128];
    time_t rawtime;
    struct tm *info;
    int psd_order = 0;

    #define EXIT() do { _mpi_stop(cores); return; } while(0)

    /* parse command line options */
    while (1)
    {
        struct option long_options[] = {
            { "help",        no_argument,       0, 'h' },
            { "verbose",     no_argument,       0, 'v' },
            { "fcqs",        no_argument,       0, 'F' },
            { "version",     no_argument,       0, 'V' },
            { "ht",          no_argument,       0, 'H' },
            { "temperature", required_argument, 0, 'T' },
            { "eta",         required_argument, 0, 'E' },
            { "ldim",        required_argument, 0, 'l' },
            { "cutoff",      required_argument, 0, 'c' },
            { "epsrel",      required_argument, 0, 'e' },
            { "iepsrel",     required_argument, 0, 'i' },
            { "material",    required_argument, 0, 'f' },
            { "omegap",      required_argument, 0, 'w' },
            { "gamma",       required_argument, 0, 'g' },
            { "psd",         required_argument, 0, 'P' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long(argc, argv, "R:L:T:l:c:e:E:f:i:w:g:P:FvVHh", long_options, &option_index);

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
            case 'H':
                ht = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'f':
                strncpy(filename, optarg, sizeof(filename)-sizeof(char));
                break;
            case 'P':
                psd_order = atoi(optarg);
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
    if(psd_order < 0)
    {
        fprintf(stderr, "order for Pade spectrum decomposition must be positive.\n\n");
        usage(stderr);
        EXIT();
    }
    if(strlen(filename))
    {
        if(strlen(filename) > 511)
        {
            fprintf(stderr, "filename %s must not be longer than 511 characters.\n\n", filename);
            usage(stderr);
            EXIT();
        }
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
        ldim = MAX(LDIM_MIN, ceil(eta/LbyR));

    /* disable buffering */
    disable_buffering();

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(time_str, sizeof(time_str),"%c", info);

    casimir_build(stdout, "# ");
    printf("# pid: %d\n", (int)getpid());
    printf("# start time: %s\n", time_str);
    printf("#\n");

    printf("# LbyR = %.15g\n", LbyR);
    printf("# RbyL = %.15g\n", 1/LbyR);
    printf("# L = %.15g\n", L);
    printf("# R = %.15g\n", R);
    if(ht)
        printf("# high-temperature limit\n");
    else
        printf("# T = %.15g\n", T);
    printf("# cutoff = %g\n", cutoff);
    printf("# epsrel = %g\n", epsrel);
    printf("# iepsrel = %g\n", iepsrel);
    printf("# ldim = %d\n", ldim);
    printf("# cores = %d\n", cores);
    if(strlen(filename))
        printf("# filename = %s\n", filename);
    else if(!isinf(omegap))
    {
        printf("# omegap = %.15g\n", omegap);
        printf("# gamma = %.15g\n", gamma_);
    }
    if(psd_order)
        printf("# use Pade spectrum decomposition of order %d\n", psd_order);

    casimir_mpi_t *casimir_mpi = casimir_mpi_init(L, R, T, filename, omegap, gamma_, ldim, cutoff, iepsrel, cores, verbose);

    /* high-temperature limit */
    if(ht)
    {
        double drude, pr, plasma;

        if(!isinf(omegap))
        {
            F_HT(casimir_mpi, omegap, &drude, &pr, &plasma);
            printf("#\n");
            printf("# L/R, L, R, ldim, omegap, F_Drude/(kB*T), F_PR/(kB*T), F_Plasma/(kB*T)\n");
            printf("%.16g, %.16g, %.16g, %d, %g, %.16g, %.16g, %.16g\n", LbyR, L, R, ldim, omegap, drude, pr, plasma);
        }
        else
        {
            F_HT(casimir_mpi, 0, &drude, &pr, NULL);
            printf("#\n");
            printf("# L/R, L, R, ldim, F_Drude/(kB*T), F_PR/(kB*T)\n");
            printf("%.16g, %.16g, %.16g, %d, %.16g, %.16g\n", LbyR, L, R, ldim, drude, pr);
        }

        casimir_mpi_free(casimir_mpi);
        return;
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
            integral = fcqs_semiinf(integrand, casimir_mpi, &epsrel, &neval, 1, &ier);
        else
        {
            integral = dqagi(integrand, 0, 1, 0, epsrel, &abserr, &neval, &ier, casimir_mpi);
            epsrel = fabs(abserr/integral);
        }

        printf("#\n");
        printf("# ier=%d, integral=%.15g, neval=%d, epsrel=%g\n", ier, integral, neval, epsrel);

        WARN(ier != 0, "ier=%d", ier);

        /* free energy for T=0 */
        F = integral/casimir_mpi->alpha/M_PI;
    }
    else
    {
        /* finite temperature */

        double drude_HT = NAN, plasma_HT = NAN, pr_HT = NAN;
        const double T_scaled = 2*M_PI*CASIMIR_kB*(R+L)*T/(CASIMIR_hbar*CASIMIR_c);
        double v[4096] = { 0 };

        /* xi = 0 */
        {
            const double t0 = now();

            casimir_t *casimir = casimir_init(LbyR);
            casimir_set_ldim(casimir, ldim);

            if(material == NULL && isinf(omegap))
            {
                F_HT(casimir_mpi, 0, NULL, &pr_HT, NULL);
                printf("# model = perfect reflectors\n");
                v[0] = pr_HT;
            }
            else if(material == NULL && gamma_ == 0)
            {
                F_HT(casimir_mpi, omegap, NULL, NULL, &plasma_HT);
                printf("# model = plasma\n");
                v[0] = plasma_HT;
            }
            else if(material == NULL)
            {
                F_HT(casimir_mpi, 0, &drude_HT, NULL, NULL);
                printf("# model = drude\n");
                v[0] = drude_HT;
            }
            else
            {
                double omegap_low, gamma_low;
                material_get_extrapolation(material, &omegap_low, &gamma_low, NULL, NULL);

                if(gamma_low == 0)
                {
                    F_HT(casimir_mpi, omegap_low*CASIMIR_hbar_eV, NULL, NULL, &plasma_HT);
                    printf("# model = optical data (xi=0: Plasma)\n");
                    v[0] = plasma_HT;
                }
                else
                {
                    F_HT(casimir_mpi, omegap_low*CASIMIR_hbar_eV, &drude_HT, NULL, &plasma_HT);
                    printf("# model = optical data (xi=0: Drude)\n");
                    printf("# plasma = %.15g (logdetD(xi=0) for plasma model with omegap=%geV)\n", plasma_HT, omegap_low*CASIMIR_hbar_eV);

                    v[0] = drude_HT;
                }
            }

            casimir_free(casimir);

            printf("#\n");
            printf("# xi=0, logdetD=%.15g, t=%g\n", v[0], now()-t0);
        }


        if(psd_order)
        {
            /* Pade spectrum decomposition (PSD), see psd.c and [1]
             * References:
             *      [1] Hu, Xu, Yan, J. Chem. Phys. 133, 101106 (2010)
             */
            double *psd_xi  = xcalloc(psd_order, sizeof(double));
            double *psd_eta = xcalloc(psd_order, sizeof(double));

            psd(psd_order, psd_xi, psd_eta);

            for(int n = 0; n < psd_order; n++)
            {
                const double xi = psd_xi[n]*T_scaled/(2*M_PI);
                const double t0 = now();
                v[n+1] = psd_eta[n]*F_xi(xi, casimir_mpi);
                printf("# xi=%.15g, logdetD=%.15g, t=%g\n", xi, v[n+1], now()-t0);
            }

            xfree(psd_xi);
            xfree(psd_eta);
        }
        else
        {
            /* Matsubara spectrum decomposition (MSD) */
            for(size_t n = 1; n < sizeof(v)/sizeof(double); n++)
            {
                const double t0 = now();
                const double xi = n*T_scaled;
                v[n] = F_xi(xi, casimir_mpi);
                printf("# xi=%.15g, logdetD=%.15g, t=%g\n", xi, v[n], now()-t0);

                if(fabs(v[n]/v[0]) < epsrel)
                    break;
            }
        }

        v[0] /= 2; /* half weight */
        F = T_scaled/M_PI*kahan_sum(v, sizeof(v)/sizeof(v[0]));
    }

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(time_str, sizeof(time_str), "%c", info);

    printf("#\n");
    printf("# %d determinants computed\n", casimir_get_determinants(casimir_mpi));
    printf("# stop time: %s\n", time_str);
    printf("#\n");
    printf("# L/R, L, R, T, ldim, E*(L+R)/(ħc)\n");
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
        double logdet = NAN;

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

        /* high-temperature case */
        if(xi == 0)
        {
            if(isinf(omegap))
                /* MM mode of PR */
                casimir_logdetD0(casimir, m, 0, NULL, &logdet, NULL);
            else
                /* plasma */
                casimir_logdetD0(casimir, m, omegap, NULL, NULL, &logdet);
        }
        else
        {
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

            logdet = casimir_logdetD(casimir, xi, m);
            TERMINATE(isnan(logdet), "L/R=%.10g, xi=%.15g, m=%d, ldim=%d", LbyR, xi, m, ldim);
        }

        MPI_Isend(&logdet, 1, MPI_DOUBLE, 0, 0, master_comm, &request);
        casimir_free(casimir);

        if(material != NULL)
            material_free(material);
        MPI_Wait(&request, &status);
    }
}

void usage(FILE *stream)
{
    fprintf(stream,
"Usage: casimir [OPTIONS]\n\n"
"This program computes the free Casimir energy E(T,L,R) for the plane-sphere\n"
"geometry. L denotes the smallest separation between sphere and plane, R is\n"
"the radius of the sphere, and T is the temperature.\n"
"\n"
"This program uses MPI for parallization and needs at least two cores to run.\n"
"\n"
"The free energy at T=0 is calculated using integration:\n"
"   E(L,R,T=0) = ∫dξ log det(1-M(ξ)),  ξ=0...∞,\n"
"where M denotes the round-trip operator. The integration is performed either\n"
"using an adaptive Gauß-Kronrod quadrature or a Fourier-Chebshev quadrature\n"
"scheme.\n"
"\n"
"References:\n"
"  [1] Hartmann, The Casimir effect in the plane-sphere geometry and the\n"
"      Proximity Force Approximation, phd thesis, 2018\n"
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
"        Set temperature to TEMPERATURE. (default: 0; in K)\n"
"\n"
"    -l, --ldim LDIM\n"
"        Set ldim to the value LDIM. (default: ldim=max(%d, ceil(eta*R/L)))\n"
"\n"
"    --eta ETA\n"
"        Set eta to the value ETA. eta is used to determine ldim if not set\n"
"        by --ldim. (default: eta=%g)\n"
"\n"
"    -c, --cutoff CUTOFF\n"
"        Stop summation over m for a given value of ξ if\n"
"            logdet(1-M^m(ξ))/logdet(1-M^0(ξ) < CUTOFF.\n"
"        (default: %g)\n"
"\n"
"    -e, --epsrel EPSREL\n"
"       Request relative accuracy of EPSREL for integration over xi if T=0,\n"
"       or stop criterion logdetD(n)/logdetD(n=0) < EPSREL for T>0.\n"
"       (default: %g)\n"
"\n"
"    -i, --iepsrel IEPSREL\n"
"       Set relative accuracy of integration over k for the matrix elements\n"
"       to IEPSREL. (default: %g)\n"
"\n"
"    -F, --fcqs\n"
"      Use Fourier-Chebshev quadrature scheme to compute integral over xi. This\n"
"      is usually faster than using Gauss-Kronrod. (only for T=0; experimental)\n"
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
"    -H, --ht\n"
"        Compute the Casimir free energy in the high-temperature limit for\n"
"        perfect reflectors, Drude and plasma model. The value for the plasma\n"
"        model is only computed if a plasma frequency is given by --omegap.\n"
"\n"
"    -P, --psd N\n"
"        Instead of the Matsubara spectrum decomposition, use the Pade spectrum\n"
"        decomposition (PSD) of order N. The PSD might converge faster than the\n"
"        but there are no checks wether the sum over the PSD frequencies\n"
"        converged. (experimental)\n"
"\n"
"    -v, --verbose\n"
"        Also print results for each m.\n"
"\n"
"    -V, --version\n"
"        Print information about build to stdout and exit.\n"
"\n"
"    -h, --help\n"
"        Show this help.\n",
    LDIM_MIN, ETA, CUTOFF, EPSREL, CASIMIR_EPSREL);
}

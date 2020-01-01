#define _DEFAULT_SOURCE /* make usleep work */

#include <getopt.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

#include "buf.h"
#include "fcqs.h"
#include "libcaps.h"
#include "material.h"
#include "misc.h"
#include "psd.h"
#include "quadpack.h"
#include "utils.h"

/* preprocessor define */
#define EPSREL 1e-6 /**< default value for --epsrel */
#define CUTOFF 1e-9 /**< default value for --cutoff */
#define LDIM_MIN 20 /**< minimum value for --ldim */
#define ETA 7.      /**< default value for --eta */
#define IDLE 1      /**< idle time in ms */

#define STATE_RUNNING 1
#define STATE_IDLE    0

/* local typedefs */

typedef struct {
    int index, m;
    double xi_;
    double recv;
    double value;
    MPI_Request request;
    int state;
} caps_task_t;

typedef struct {
    double L, R, T, omegap, gamma, cutoff, iepsrel, alpha;
    int ldim, cores;
    bool verbose;
    caps_task_t **tasks;
    int determinants;
    char filename[512];
    double cache[4096][2];
    int cache_elems;
} caps_mpi_t;

/* prototypes */
static caps_mpi_t *caps_mpi_init(double L, double R, double T, char *filename, char *resume, double omegap, double gamma_, int ldim, double cutoff, double iepsrel, int cores, bool verbose);
static void caps_mpi_free(caps_mpi_t *self);
static int caps_mpi_submit(caps_mpi_t *self, int index, double xi, int m);
static int caps_mpi_retrieve(caps_mpi_t *self, caps_task_t **task_out);
static int caps_mpi_get_running(caps_mpi_t *self);
static int caps_get_determinants(caps_mpi_t *self);

static void usage(FILE *stream);

static void F_HT(caps_mpi_t *caps_mpi, double omegap, double *drude, double *pr, double *plasma);
static double F_xi(double xi, caps_mpi_t *caps_mpi);
static void master(int argc, char *argv[], const int cores);
static void slave(MPI_Comm master_comm, const int rank);


/** @brief Write time into string
 *
 * Write current time in a human readable format into string s. The output is
 * similar to "Aug 30 2018 14:37:35".
 *
 * @param s string
 * @param len maximum length of array s
 */
static void time_as_string(char *s, size_t len)
{
    time_t rawtime;
    struct tm *info;

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(s, len, "%c", info);
}


/* sleep for ms miliseconds */
static void msleep(unsigned int ms)
{
#ifdef _WIN32
    Sleep(ms);
#else
    usleep(ms*1000);
#endif
}

static bool is_readable(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if(f == NULL)
        return false;

    fclose(f);
    return true;
}

/* @brief Create caps_mpi object
 *
 * @param [in] L separation between sphere and plate in meter
 * @param [in] R radius of sphere in meter
 * @param [in] T temperature in Kelvin
 * @param [in] filename filename of material description or NULL
 * @param [in] resume filename of partial output to be resumed
 * @param [in] omegap plasma frequency of the Drude model in eV
 * @param [in] gamma_ relaxation frequency of the Drude model eV
 * @param [in] ldim dimension of vector space
 * @param [in] cutoff cutoff for summation over m
 * @param [in] iepsrel relative accuracy for integration of k for matrix elements
 * @param [in] cores number of cores to use
 * @param [in] verbose flag if verbose
 * @retval object caps_mpi_t object
 */
static caps_mpi_t *caps_mpi_init(double L, double R, double T, char *filename, char *resume, double omegap, double gamma_, int ldim, double cutoff, double iepsrel, int cores, bool verbose)
{
    caps_mpi_t *self = xmalloc(sizeof(caps_mpi_t));

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
    self->tasks   = xmalloc(cores*sizeof(caps_task_t *));
    self->alpha   = 2*L/(L+R); /* used to scale integration if T=0 */

    /* number of determinants we have computed */
    self->determinants = 0;

    TERMINATE(strlen(filename) > 511, "filename too long: %s", filename);
    memset(self->filename, '\0', sizeof(self->filename));
    strncpy(self->filename, filename, sizeof(self->filename)-sizeof(char));

    /* cache for resume */
    self->cache_elems = 0;
    if(resume && strlen(resume) > 0)
    {
        char line[512];
        FILE *fh = fopen(resume, "r");
        TERMINATE(fh == NULL, "cannot open %s for reading", resume);

        while(fgets(line, sizeof(line)/sizeof(char), fh) != NULL)
        {
            char *p1 = strstr(line, "# xi*(L+R)/c=");
            char *p2 = strstr(line, "logdetD=");
            if(p1 && p2)
            {
                char *p3;

                p1 += 13;
                p3 = strchr(p1, ',');
                TERMINATE(p3 == NULL, "%s has wrong format", resume);
                *p3 = '\0';
                self->cache[self->cache_elems][0] = strtodouble(p1); /* xi */

                p2 += 8;
                p3 = strchr(p2, ',');
                TERMINATE(p3 == NULL, "%s has wrong format", resume);
                *p3 = '\0';
                self->cache[self->cache_elems][1] = strtodouble(p2); /* logdetD */

                self->cache_elems++;
            }
        }

        fclose(fh);
    }

    self->tasks[0] = NULL;
    for(int i = 1; i < cores; i++)
    {
        caps_task_t *task = xmalloc(sizeof(caps_task_t));
        task->index    = -1;
        task->state    = STATE_IDLE;
        self->tasks[i] = task;
    }

    return self;
}

/* stop all remaining slaves */
static void _mpi_stop(int cores)
{
    double buf[] = { -1, -1, -1, -1, -1, -1, -1 };

    for(int i = 1; i < cores; i++)
        MPI_Send(buf, 7, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
}

/** @brief Free mpi object
 *
 * Stop all running jobs and free allocated memory.
 *
 * @param [in] self caps_mpit_t object
 */
static void caps_mpi_free(caps_mpi_t *self)
{
    _mpi_stop(self->cores);

    for(int i = 1; i < self->cores; i++)
        xfree(self->tasks[i]);

    xfree(self->tasks);
    xfree(self);
}


/** @brief Get number of running jobs
 *
 * @param [in] self caps_mpit_t object
 * @retval running number of processes that are running
 */
static int caps_mpi_get_running(caps_mpi_t *self)
{
    int running = 0;

    for(int i = 1; i < self->cores; i++)
        if(self->tasks[i]->state == STATE_RUNNING)
            running++;

    return running;
}

/* xi_ = ξ(L+R)/c */
static int caps_mpi_submit(caps_mpi_t *self, int index, double xi_, int m)
{
    for(int i = 1; i < self->cores; i++)
    {
        caps_task_t *task = self->tasks[i];

        if(task->state == STATE_IDLE)
        {
            double buf[] = { xi_, self->L, self->R, self->omegap, self->gamma, m, self->iepsrel, self->ldim };

            task->index = index;
            task->xi_   = xi_;
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

/** @brief Get number of computed determinants
 *
 * Get the number of determinants that have been computed.
 *
 * @param [in] self caps_mpit_t object
 * @retval determinants number of computed determinants
 */
static int caps_get_determinants(caps_mpi_t *self)
{
    return self->determinants;
}

static int caps_mpi_retrieve(caps_mpi_t *self, caps_task_t **task_out)
{
    *task_out = NULL;

    for(int i = 1; i < self->cores; i++)
    {
        caps_task_t *task = self->tasks[i];

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

/* x = 2ξL/c */
static double integrand(double x, void *args)
{
    double logdetD;
    const double t0 = now();
    caps_mpi_t *caps_mpi = (caps_mpi_t *)args;
    const double xi_ = x/caps_mpi->alpha; /* xi_=ξ(L+R)/c; α=2*L/(L+R) */

    /* For large values of ξ the integrand logdet(Id-M(ξ)) is almost 0, but the
     * actual computation of the matrix elements might yield warnings and
     * errors. To save computation time and prevent warnings and errors, we
     * estimate the integrand using the PFA assuming perfect reflectors. If the
     * aspect ratio is sufficiently high to use the PFA estimate (we choose
     * R/L>10), we estimate the value of the Matsubara frequency ξ where
     *      logdet(Id-M(ξ_cutoff))=logdet_cutoff=1e-100.
     * If ξ>ξ_cutoff and R/L>10, we return 0.
     *
     * Assuming perfect reflectors and that the PFA is valid, one finds
     *      logdet(Id-M(ξ)) = -1/2 R/L Li_3(exp(-ξL/c)).
     * As ξ is assumed to be large, the argument of the polylog becomes small,
     * and we can use Li_3(x)=~x. With this, one finds the estimate above.
     */
    const double LbyR = caps_mpi->L/caps_mpi->R;
    const double logdet_cutoff = 1e-100;
    const double xi_cutoff = -(1+1/LbyR)*log(2*LbyR*logdet_cutoff)/2;

    if(LbyR < 0.1 && xi_ > xi_cutoff)
        logdetD = 0;
    else
        logdetD = F_xi(xi_, caps_mpi);

    printf("# xi*(L+R)/c=%.16g, logdetD=%.16g, t=%g\n", xi_, logdetD, now()-t0);
    return logdetD;
}

/* omegap in eV; drude, pr and plasma in units of kB*T */
static void F_HT(caps_mpi_t *caps_mpi, double omegap, double *drude, double *pr, double *plasma)
{
    const double omegap_orig = caps_mpi->omegap;
    const double gamma_orig  = caps_mpi->gamma;

    /* Drude
     * The actual value of omegap and gamma for the high-temperature limit are
     * irrelevant. What is relevant is that gamma is positive and non-zero.
     */
    if(drude != NULL)
    {
        caps_mpi->omegap = 1;
        caps_mpi->gamma  = 1;
        *drude = F_xi(0, caps_mpi);
    }

    /* PR */
    if(pr != NULL)
    {
        caps_mpi->omegap = INFINITY;
        caps_mpi->gamma  = 0;
        *pr = F_xi(0, caps_mpi);
    }

    /* plasma */
    if(plasma != NULL)
    {
        caps_mpi->omegap = omegap;
        caps_mpi->gamma  = 0;
        *plasma = F_xi(0, caps_mpi);
    }

    caps_mpi->omegap = omegap_orig;
    caps_mpi->gamma  = gamma_orig;
}

/* xi_ = ξ(L+R)/c */
static double F_xi(double xi_, caps_mpi_t *caps_mpi)
{
    int m;
    double drude_HT = NAN;
    double terms[4096] = { NAN };
    bool verbose = caps_mpi->verbose;
    const double mmax = sizeof(terms)/sizeof(double);
    const double cutoff = caps_mpi->cutoff;

    /* look in cache if resumed */
    for(int j = 0; j < caps_mpi->cache_elems; j++)
    {
        const double xi_cache = caps_mpi->cache[j][0];
        //printf("xi_cache=%g vs xi=%g\n", xi_cache, xi);

        if(xi_ == xi_cache || fabs(1-xi_cache/xi_) < 1e-11)
        {
            const double logdetD = caps_mpi->cache[j][1];
            return logdetD;
        }
    }

    if(xi_ == 0)
    {
        /* compute Drude contribution */
        caps_t *caps = caps_init(caps_mpi->R, caps_mpi->L);
        caps_set_ldim(caps, caps_mpi->ldim);
        if(caps_mpi->iepsrel > 0)
            caps_set_epsrel(caps, caps_mpi->iepsrel);
        drude_HT = caps_ht_drude(caps);
        caps_free(caps);

        if(!isinf(caps_mpi->omegap) && caps_mpi->gamma > 0)
            /* omegap finite => Drude model */
            return drude_HT;
    }

    /* gather all data */
    for(m = 0; m < mmax; m++)
    {
        while(1)
        {
            caps_task_t *task = NULL;

            /* send job */
            if(caps_mpi_submit(caps_mpi, m, xi_, m))
                break;

            /* retrieve jobs */
            while(caps_mpi_retrieve(caps_mpi, &task))
            {
                double v = terms[task->m] = task->value;

                if(verbose)
                    fprintf(stderr, "# m=%d, xi_=%.16g, logdetD=%.16g\n", task->m, xi_, task->value);

                if(v == 0 || v/terms[0] < cutoff)
                    goto done;
            }

            msleep(IDLE);
        }
    }

    TERMINATE(true, "sum did not converge, sorry. :(");

    done:

    /* retrieve all remaining running jobs */
    while(caps_mpi_get_running(caps_mpi) > 0)
    {
        caps_task_t *task = NULL;

        while(caps_mpi_retrieve(caps_mpi, &task))
        {
            terms[task->m] = task->value;
            if(verbose)
                fprintf(stderr, "# m=%d, xi_=%.16g, logdetD=%.16g\n", task->m, xi_, task->value);
        }

        msleep(IDLE);
    }

    terms[0] /= 2; /* m = 0 */

    if(xi_ == 0)
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

static void master(int argc, char *argv[], const int cores)
{
    bool verbose = false, fcqs = false, ht = false;
    char filename[512] = { 0 };
	char resume[512] = { 0 };
    int ldim = 0;
    double L = 0, R = 0, T = 0, omegap = INFINITY, gamma_ = 0;
    double cutoff = CUTOFF, epsrel = EPSREL, eta = ETA;
    double iepsrel = CAPS_EPSREL;
    material_t *material = NULL;
    char time_str[128];
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
            { "psd",         no_argument,       0, 'p' },
            { "temperature", required_argument, 0, 'T' },
            { "eta",         required_argument, 0, 'E' },
            { "ldim",        required_argument, 0, 'l' },
            { "cutoff",      required_argument, 0, 'c' },
            { "epsrel",      required_argument, 0, 'e' },
            { "iepsrel",     required_argument, 0, 'i' },
            { "material",    required_argument, 0, 'f' },
			{ "resume",      required_argument, 0, 'r' },
            { "omegap",      required_argument, 0, 'w' },
            { "gamma",       required_argument, 0, 'g' },
            { "psd-order",   required_argument, 0, 'P' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long(argc, argv, "R:L:T:l:c:e:E:f:r:i:w:g:P:pFvVHh", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'L':
                L = strtodouble(optarg);
                break;
            case 'R':
                R = strtodouble(optarg);
                break;
            case 'T':
                T = strtodouble(optarg);
                break;
            case 'l':
                ldim = atoi(optarg);
                break;
            case 'E':
                eta = strtodouble(optarg);
                break;
            case 'c':
                cutoff = strtodouble(optarg);
                break;
            case 'i':
                iepsrel = strtodouble(optarg);
                break;
            case 'e':
                epsrel = strtodouble(optarg);
                break;
            case 'w':
                omegap = strtodouble(optarg);
                break;
            case 'g':
                gamma_ = strtodouble(optarg);
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
            case 'r':
                strncpy(resume, optarg, sizeof(resume)-sizeof(char));
                break;
            case 'p':
                psd_order = -1; /* auto */
                break;
            case 'P':
                psd_order = atoi(optarg);
                break;
            case 'V':
                caps_build(stdout, NULL);
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
        if(strlen(filename) > 511)
        {
            fprintf(stderr, "filename %s must not be longer than 511 characters.\n\n", filename);
            usage(stderr);
            EXIT();
        }
        if(!is_readable(filename))
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
    if(fcqs && (ht || T > 0))
    {
        fprintf(stderr, "flag --fcqs can only be used when T=0\n\n");
        usage(stderr);
        EXIT();
    }

    if(cores < 2)
    {
        fprintf(stderr, "This program needs at least 2 cores to run.\n");
        fprintf(stderr, "Have you started the program using mpirun?\n");
        fprintf(stderr, "Example using five cores:\n");
        fprintf(stderr, "  $ mpirun -n 5 %s [OPTIONS]\n", argv[0]);
        EXIT();
    }

    const double LbyR = L/R;

    if(psd_order < 0)
    {
        const double Teff = 4*CAPS_PI*CAPS_KB/CAPS_HBAR/CAPS_C*T*L;
        psd_order = ceil( (1-1.5*log10(epsrel))/sqrt(Teff) );
    }

    /* if ldim was not set */
    if(ldim <= 0)
        ldim = MAX(LDIM_MIN, ceil(eta/LbyR));

    /* disable buffering */
    disable_buffering();

    time_as_string(time_str, sizeof(time_str)/sizeof(time_str[0]));

    caps_build(stdout, "# ");
#ifndef _WIN32
    printf("# pid: %d\n", getpid());
#endif
    printf("# start time: %s\n", time_str);
    printf("#\n");

    printf("# LbyR = %.16g\n", LbyR);
    printf("# RbyL = %.16g\n", 1/LbyR);
    printf("# L = %.16g\n", L);
    printf("# R = %.16g\n", R);
    if(ht)
        printf("# high-temperature limit\n");
    else
    {
        printf("# T = %.16g\n", T);
        if(T > 0)
        {
            if(psd_order)
                printf("# using Pade spectrum decomposition (PSD) of order %d\n", psd_order);
            else
                printf("# using Matsubara spectrum decomposition (MSD)\n");
        }
    }
    printf("# cutoff = %g\n", cutoff);
    printf("# epsrel = %g\n", epsrel);
    printf("# iepsrel = %g\n", iepsrel);
    printf("# ldim = %d\n", ldim);
    printf("# cores = %d\n", cores);
    if(strlen(filename))
        printf("# filename = %s\n", filename);
    else if(!isinf(omegap))
    {
        printf("# omegap = %.16g\n", omegap);
        printf("# gamma = %.16g\n", gamma_);
    }
	if(strlen(resume))
        printf("# resume = %s\n", resume);

    caps_mpi_t *caps_mpi = caps_mpi_init(L, R, T, filename, resume, omegap, gamma_, ldim, cutoff, iepsrel, cores, verbose);

    /* high-temperature limit */
    if(ht)
    {
        double drude, pr, plasma;

        if(!isinf(omegap))
        {
            F_HT(caps_mpi, omegap, &drude, &pr, &plasma);

            time_as_string(time_str, sizeof(time_str)/sizeof(time_str[0]));
            printf("#\n");
            printf("# stop time: %s\n", time_str);
            printf("#\n");
            printf("# L/R, L, R, ldim, omegap, E_Drude/(kB*T), E_PR/(kB*T), E_Plasma/(kB*T)\n");
            printf("%.16g, %.16g, %.16g, %d, %g, %.16g, %.16g, %.16g\n", LbyR, L, R, ldim, omegap, drude, pr, plasma);
        }
        else
        {
            F_HT(caps_mpi, 0, &drude, &pr, NULL);

            time_as_string(time_str, sizeof(time_str)/sizeof(time_str[0]));
            printf("#\n");
            printf("# stop time: %s\n", time_str);
            printf("#\n");
            printf("# L/R, L, R, ldim, E_Drude/(kB*T), E_PR/(kB*T)\n");
            printf("%.16g, %.16g, %.16g, %d, %.16g, %.16g\n", LbyR, L, R, ldim, drude, pr);
        }

        caps_mpi_free(caps_mpi);
        return;
    }

    double F = NAN;
    if(T == 0)
    {
        double integral = 0, abserr = 0;
        int ier, neval;

        if(fcqs)
            printf("# quad = Fourier-Chebshev quadrature scheme\n");
        else
            printf("# quad = adaptive Gauss-Kronrod\n");
        printf("#\n");

        if(fcqs)
            integral = fcqs_semiinf(integrand, caps_mpi, &epsrel, &neval, 1, &ier);
        else
        {
            integral = dqagi(integrand, 0, 1, 0, epsrel, &abserr, &neval, &ier, caps_mpi);
            epsrel = fabs(abserr/integral);
        }

        printf("#\n");
        printf("# ier=%d, integral=%.16g, neval=%d, epsrel=%g\n", ier, integral, neval, epsrel);

        WARN(ier != 0, "ier=%d", ier);

        /* free energy for T=0 */
        F = integral/caps_mpi->alpha/CAPS_PI;
    }
    else
    {
        /* finite temperature */
        double drude_HT = NAN, plasma_HT = NAN, pr_HT = NAN;
        const double T_scaled = 2*CAPS_PI*CAPS_KB*(R+L)*T/(CAPS_HBAR*CAPS_C);
        double *v = NULL;

        /* xi = 0 */
        {
            const double t0 = now();

            caps_t *caps = caps_init(R, L);
            caps_set_ldim(caps, ldim);

            if(material == NULL && isinf(omegap))
            {
                F_HT(caps_mpi, 0, NULL, &pr_HT, NULL);
                printf("# model = perfect reflectors\n");
                buf_push(v, pr_HT);
            }
            else if(material == NULL && gamma_ == 0)
            {
                F_HT(caps_mpi, omegap, NULL, NULL, &plasma_HT);
                printf("# model = plasma\n");
                buf_push(v, plasma_HT);
            }
            else if(material == NULL)
            {
                F_HT(caps_mpi, 0, &drude_HT, NULL, NULL);
                printf("# model = drude\n");
                buf_push(v, drude_HT);
            }
            else
            {
                double omegap_low, gamma_low;
                material_get_extrapolation(material, &omegap_low, &gamma_low, NULL, NULL);
                omegap_low *= CAPS_HBAR_EV; /* convert from rad/s to eV */
                gamma_low  *= CAPS_HBAR_EV; /* convert from rad/s to eV */

                if(gamma_low == 0)
                {
                    F_HT(caps_mpi, omegap_low, NULL, NULL, &plasma_HT);
                    printf("# model = optical data (xi=0: Plasma)\n");
                    buf_push(v, plasma_HT);
                }
                else
                {
                    F_HT(caps_mpi, omegap_low, &drude_HT, NULL, &plasma_HT);
                    printf("# model = optical data (xi=0: Drude)\n");
                    printf("# plasma = %.16g (logdetD(xi=0) for plasma model with omegap=%geV)\n", plasma_HT, omegap_low);

                    buf_push(v, drude_HT);
                }
            }

            caps_free(caps);

            printf("#\n");
            printf("# xi*(L+R)/c=0, logdetD=%.16g, t=%g\n", v[0], now()-t0);
        }


        if(psd_order > 0)
        {
            /* Pade spectrum decomposition (PSD), see psd.c.
             * Reference: Hu, Xu, Yan, J. Chem. Phys. 133, 101106 (2010),
             *            https://doi.org/10.1063/1.3602466
             */
            double *psd_xi  = xcalloc(psd_order, sizeof(double));
            double *psd_eta = xcalloc(psd_order, sizeof(double));

            psd(psd_order, psd_xi, psd_eta);

            for(int n = 0; n < psd_order; n++)
            {
                const double xi = psd_xi[n]*T_scaled/(2*CAPS_PI);
                const double t0 = now();
                buf_push(v, psd_eta[n]*F_xi(xi, caps_mpi));
                printf("# xi*(L+R)/c=%.16g, logdetD=%.16g, t=%g\n", xi, v[n+1], now()-t0);
            }

            xfree(psd_xi);
            xfree(psd_eta);
        }
        else
        {
            /* Matsubara spectrum decomposition (MSD) */
            for(size_t n = 1; true; n++)
            {
                const double t0 = now();
                const double xi = n*T_scaled;
                buf_push(v, F_xi(xi, caps_mpi));
                printf("# xi*(L+R)/c=%.16g, logdetD=%.16g, t=%g\n", xi, v[n], now()-t0);

                if(fabs(v[n]/v[0]) < epsrel)
                    break;
            }
        }

        v[0] /= 2; /* half weight */
        F = T_scaled/CAPS_PI*kahan_sum(v, buf_size(v));

        buf_free(v);
    }

    time_as_string(time_str, sizeof(time_str)/sizeof(time_str[0]));

    printf("#\n");
    printf("# %d determinants computed\n", caps_get_determinants(caps_mpi));
    printf("# stop time: %s\n", time_str);
    printf("#\n");
    printf("# L/R, L, R, T, ldim, E*(L+R)/(hbar*c)\n");
    printf("%.16g, %.16g, %.16g, %.16g, %d, %.16g\n", LbyR, L, R, T, ldim, F);

    if(material != NULL)
        material_free(material);

    caps_mpi_free(caps_mpi);
}

static void slave(MPI_Comm master_comm, const int rank)
{
    char filename[512] = { 0 };
    double userdata[2] = { 0 };
    double buf[8] = { 0 };

    while(1)
    {
        double logdet = NAN;

        caps_t  *caps  = NULL;
        material_t *material = NULL;

        memset(buf,      0, sizeof(buf));
        memset(filename, 0, sizeof(filename));
        memset(userdata, 0, sizeof(userdata));

        MPI_Status status;
        MPI_Recv(buf, 8, MPI_DOUBLE, 0, 0, master_comm, &status);

        /* Matsubara frequency; xi_ = ξ(L+R)/c */
        const double xi_ = buf[0];

        /* signal to quit */
        if(xi_ < 0)
            break;

        /* geometry */
        const double L = buf[1]; /* in m */;
        const double R = buf[2]; /* in m */
        const double LbyR = L/R;

        const double omegap = buf[3]/CAPS_HBAR_EV; /* plasma frequency in rad/s */
        const double gamma_ = buf[4]/CAPS_HBAR_EV; /* relaxation frequency in rad/s */

        const int m          = (int)buf[5];
        const double iepsrel = buf[6];
        const int ldim       = (int)buf[7];

        /* get filename */
        MPI_Recv(filename, 512, MPI_CHAR, 0, 0, master_comm, &status);

        caps = caps_init(R,L);
        TERMINATE(caps == NULL, "caps object is null");
        caps_set_ldim(caps, ldim);

        if(iepsrel > 0)
            caps_set_epsrel(caps, iepsrel);

        /* high-temperature case */
        if(xi_ == 0)
        {
            if(isinf(omegap))
                /* MM mode of PR */
                caps_logdetD0(caps, m, 0, NULL, &logdet, NULL);
            else
                /* plasma */
                caps_logdetD0(caps, m, omegap, NULL, NULL, &logdet);
        }
        else
        {
            /* set material properties */
            if(strlen(filename))
            {
                material = material_init(filename, L+R);
                TERMINATE(material == NULL, "material_init failed");
                caps_set_epsilonm1(caps, material_epsilonm1, material);
            }
            else if(!isinf(omegap))
            {
                userdata[0] = omegap;
                userdata[1] = gamma_;
                caps_set_epsilonm1(caps, caps_epsilonm1_drude, userdata);
            }

            logdet = caps_logdetD(caps, xi_, m);
            TERMINATE(isnan(logdet), "L/R=%.16g, xi_=%.16g, m=%d, ldim=%d", LbyR, xi_, m, ldim);
        }

        caps_free(caps);
        if(material != NULL)
            material_free(material);

        MPI_Request request;
        MPI_Isend(&logdet, 1, MPI_DOUBLE, 0, 0, master_comm, &request);
        MPI_Wait(&request, &status);
    }
}

static void usage(FILE *stream)
{
    fprintf(stream,
"Usage: caps [OPTIONS]\n\n"
"This program computes the free Casimir energy E(T,L,R) for the plane-sphere\n"
"geometry. L denotes the smallest separation between sphere and plane, R is the\n"
"radius of the sphere, and T is the temperature.\n"
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
"  [2] Hartmann, Ingold, Maia Neto, Advancing numerics for the Casimir effect\n"
"      to experimentally relevant aspect ratios, Phys. Scr. 93, 114003 (2018)\n"
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
"        Set eta to the value ETA. eta is used to determine ldim if not set by\n"
"        --ldim. (default: eta=%g)\n"
"\n"
"    -c, --cutoff CUTOFF\n"
"        Stop summation over m for a given value of ξ if\n"
"            logdet(1-M^m(ξ))/logdet(1-M^0(ξ) < CUTOFF.\n"
"        (default: %g)\n"
"\n"
"    -e, --epsrel EPSREL\n"
"       Request relative accuracy of EPSREL for integration over xi if T=0, or\n"
"       stop criterion logdetD(n)/logdetD(n=0) < EPSREL for T>0.\n"
"       (default: %g)\n"
"\n"
"    -i, --iepsrel IEPSREL\n"
"       Set relative accuracy of integration over k for the matrix elements to\n"
"       IEPSREL. (default: %g)\n"
"\n"
"    -F, --fcqs\n"
"      Use Fourier-Chebshev quadrature scheme to compute integral over xi. This\n"
"      is usually faster than using Gauss-Kronrod. (only for T=0; experimental)\n"
"\n"
"    -f, --material FILENAME\n"
"        Filename of the material description file. If set, --omegap and\n"
"        --gamma will be ignored.\n"
"\n"
"    -r, --resume FILENAME\n"
"        Resume the computation from a partial output file. (experimental)\n"
"\n"
"    --omegap OMEGAP\n"
"        Model the metals using the Drude/Plasma model and set plasma\n"
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
"    -p, --psd\n"
"        Instead of the Matsubara spectrum decomposition, use the Pade spectrum\n"
"        decomposition (PSD) of order N. N is chosen automatically. The PSD\n"
"        converges faster than the MSD. (experimental)\n"
"\n"
"    -P, --psd-order N\n"
"        Use Pade spectrum decomposition (PSD) of order N. In contrast to the\n"
"        option --psd, you can chose the order N of the PSD. (experimental)\n"
"\n"
"    -v, --verbose\n"
"        Also print results for each m.\n"
"\n"
"    -V, --version\n"
"        Print information about build to stdout and exit.\n"
"\n"
"    -h, --help\n"
"        Show this help.\n",
    LDIM_MIN, ETA, CUTOFF, EPSREL, CAPS_EPSREL);
}

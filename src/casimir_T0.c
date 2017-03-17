#define _BSD_SOURCE /* make usleep work */

#include <ctype.h>
#include <getopt.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "casimir_T0.h"

#include "gausslaguerre.h"
#include "quadpack.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

#define EPSREL 1e-6
#define CUTOFF 1e-9
#define ETA 7.
#define IDLE 25

#define STATE_RUNNING 1
#define STATE_IDLE    0

static double _pfa(double LbyR, double T, double omegap, double gamma_)
{
    double userdata[2];
    casimir_t *casimir = casimir_init(LbyR);
    if(!isinf(omegap))
    {
        userdata[0] = omegap;
        userdata[1] = gamma_;
        casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);
    }
    double pfa = casimir_pfa(casimir, T);

    casimir_free(casimir);

    return pfa;
}

void casimir_mpi_init(casimir_mpi_t *self, double LbyR, double omegap, double gamma_, int ldim, double cutoff, int cores, bool verbose)
{
    self->LbyR    = LbyR;
    self->omegap  = omegap;
    self->gamma   = gamma_;
    self->ldim    = ldim;
    self->cutoff  = cutoff;
    self->cores   = cores;
    self->verbose = verbose;
    self->tasks   = xmalloc(cores*sizeof(casimir_task_t *));
    self->alpha   = 2*LbyR/(1+LbyR);
    self->k       = 1;

    for(int i = 0; i < cores; i++)
    {
        casimir_task_t *task = xmalloc(sizeof(casimir_task_t));
        task->index    = -1;
        task->state    = STATE_IDLE;
        self->tasks[i] = task;
    }
}

void casimir_mpi_free(casimir_mpi_t *self)
{
    double buf[] = { -1, -1, -1, -1 };

    /* stop all remaining slaves */
    for(int i = 1; i < self->cores; i++)
        MPI_Send(buf, 4, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

    for(int i = 0; i < self->cores; i++)
    {
        xfree(self->tasks[i]);
        self->tasks[i] = NULL;
    }

    xfree(self->tasks);
    self->tasks = NULL;
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
            double buf[] = { xi, self->LbyR, self->omegap, self->gamma, m, self->ldim };

            task->index = index;
            task->xi    = xi;
            task->m     = m;
            task->state = STATE_RUNNING;

            MPI_Send (buf,        6, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Irecv(task->recv, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &task->request);

            return 1;
        }
    }

    return 0;
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

                task->value = task->recv[0];
                task->state = STATE_IDLE;

                *task_out = task;

                return 1;
            }
        }
    }

    return 0;
}

static double wrapper_integrand(double xi, void *args)
{
    casimir_mpi_t *casimir_mpi = (casimir_mpi_t *)args;
    return integrand(xi/casimir_mpi->alpha, casimir_mpi);
}

double integrand(double xi, casimir_mpi_t *casimir_mpi)
{
    int m;
    double terms[2048] = { NAN };
    bool verbose = casimir_mpi->verbose;
    const double nmax = sizeof(terms)/sizeof(double);
    const double cutoff = casimir_mpi->cutoff;
    const double t0 = now();

    /* gather all data */
    for(m = 0; m < nmax; m++)
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
    double logdetD = kahan_sum(terms, m);

    printf("# k=%d, xi=%.12g, logdetD=%.15g, t=%g\n", casimir_mpi->k++, xi, logdetD, now()-t0);

    return logdetD;
}

int main(int argc, char *argv[])
{
    int cores, rank, ret;
    MPI_Comm new_comm;

    /* initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, rank == 0, 0, &new_comm);
    MPI_Comm_size(MPI_COMM_WORLD, &cores);

    if(rank == 0)
        ret = master(argc, argv, cores);
    else
        ret = slave(MPI_COMM_WORLD, rank);

    MPI_Finalize();

    return ret;
}

int master(int argc, char *argv[], int cores)
{
    bool verbose = false, zero = 0;
    int order = -1, ldim = 0, ret = 0;
    double LbyR = -1, cutoff = CUTOFF, epsrel = EPSREL, omegap = INFINITY, gamma_ = 0, T = 0;
    casimir_mpi_t casimir_mpi;

    #define EXIT(n) do { ret = n; goto out; } while(0)

    /* parse command line options */
    while (1)
    {
        int c;
        struct option long_options[] = {
            { "help",        no_argument,       0, 'h' },
            { "verbose",     no_argument,       0, 'v' },
            { "zero",        no_argument,       0, 'z' },
            { "zero",        no_argument,       0, 'z' },
            { "LbyR",        required_argument, 0, 'x' },
            { "temperature", required_argument, 0, 'T' },
            { "ldim",        required_argument, 0, 'L' },
            { "order",       required_argument, 0, 'N' },
            { "cutoff",      required_argument, 0, 'c' },
            { "epsrel",      required_argument, 0, 'e' },
            { "omegap",      required_argument, 0, 'w' },
            { "gamma",       required_argument, 0, 'g' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "x:T:L:N:c:e:p:w:g:zdvh", long_options, &option_index);

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
                LbyR = atof(optarg);
                break;
            case 'T':
                T = atof(optarg);
                break;
            case 'L':
                ldim = atoi(optarg);
                break;
            case 'c':
                cutoff = atof(optarg);
                break;
            case 'e':
                epsrel = atof(optarg);
                break;
            case 'N':
                order = atoi(optarg);
                break;
            case 'w':
                omegap = atof(optarg);
                break;
            case 'g':
                gamma_ = atof(optarg);
                break;
            case 'v':
                verbose = true;
                break;
            case 'z':
                zero = true;
                break;
            case 'h':
                usage(stdout);
                EXIT(0);

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    if(LbyR <= 0)
    {
        fprintf(stderr, "LbyR must be positive.\n\n");
        usage(stderr);
        EXIT(1);
    }
    if(T < 0)
    {
        fprintf(stderr, "Temperature must be non-negative.\n\n");
        usage(stderr);
        EXIT(1);
    }
    if(cutoff <= 0)
    {
        fprintf(stderr, "cutoff must be positive.\n\n");
        usage(stderr);
        EXIT(1);
    }
    if(epsrel <= 0)
    {
        fprintf(stderr, "epsrel must be positive.\n\n");
        usage(stderr);
        EXIT(1);
    }
    if(!isinf(omegap))
    {
        if(omegap <= 0)
        {
            fprintf(stderr, "omegap must be positive\n\n");
            usage(stderr);
            EXIT(1);
        }
        if(gamma_ < 0)
        {
            fprintf(stderr, "gamma must be non-negative\n\n");
            usage(stderr);
            EXIT(1);
        }
    }

    if(cores < 2)
    {
        fprintf(stderr, "need at least 2 cores\n");
        EXIT(1);
    }

    /* if ldim was not set */
    if(ldim <= 0)
        ldim = ceil(ETA/LbyR);

    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    /* estimate, cf. eq. (6.33) */
    const double alpha = 2*LbyR/(1+LbyR);

    /* compute pfa */
    double pfa = _pfa(LbyR, T, omegap, gamma_);

    printf("# LbyR   = %.15g\n", LbyR);
    printf("# T      = %.15g\n", T);
    printf("# cutoff = %g\n", cutoff);
    printf("# ldim   = %d\n", ldim);
    printf("# cores  = %d\n", cores);
    printf("# alpha  = %.15g\n", alpha);
    printf("# F_PFA  = %.15g\n", pfa);
    if(!isinf(omegap))
    {
        printf("# omegap = %.15g\n", omegap);
        printf("# gamma  = %.15g\n", gamma_);
    }

    casimir_mpi_init(&casimir_mpi, LbyR, omegap, gamma_, ldim, cutoff, cores, verbose);

    double F;
    if(T == 0)
    {
        double integral = 0;

        if(order > 0)
        {
            double *xk, *ln_wk;
            order = gausslaguerre_nodes_weights(order, &xk, &ln_wk);

            printf("# quadrature: Gauss-Laguerre (order = %d)\n", order);
            printf("#\n");

            if(zero)
                wrapper_integrand(0, &casimir_mpi);

            for(int k = 0; k < order; k++)
            {
                double xi = xk[k]/alpha;
                double v = integrand(xi, &casimir_mpi);

                integral += exp(ln_wk[k]+xk[k])*v;
            }
        }
        else
        {
            double abserr;
            int ier, neval;

            printf("# quadrature: adaptive Gauss-Kronrod (epsrel = %g)\n", epsrel);
            printf("#\n");

            if(zero)
                wrapper_integrand(0, &casimir_mpi);

            integral = dqagi(wrapper_integrand, 0, 1, 0, epsrel, &abserr, &neval, &ier, &casimir_mpi);

            printf("#\n");
            printf("# ier=%d, integral=%.15g, neval=%d, abserr=%g, absrel=%g\n", ier, integral, neval, abserr, fabs(abserr/integral));

            WARN(ier != 0, "ier=%d", ier);
        }

        /* free energy for T=0 */
        F = integral/alpha/M_PI;
    }
    else
    {
        /* fnite temperature */
        double v[4096] = { 0 };

        F = NAN;

        for(size_t n = 0; n < sizeof(v)/sizeof(double); n++)
        {
            double xi = n*T;
            v[n] = integrand(xi, &casimir_mpi);

            if(fabs(v[n]/v[0]) < epsrel)
            {
                v[0] /= 2;
                F = T/M_PI*kahan_sum(v, n+1);
                break;
            }
        }

        TERMINATE(isnan(F), "something went wrong: free energy is not a number");
    }

    printf("#\n");
    printf("# L/R, T, ldim, F_PFA*(L+R)/(ħc), F*(L+R)/(ħc), F/F_pfa\n");
    printf("%g, %.15g, %d, %.12g, %.12g, %.12g\n", LbyR, T, ldim, pfa, F, F/pfa);

    casimir_mpi_free(&casimir_mpi);
out:

    return ret;
}

int slave(MPI_Comm master_comm, int rank)
{
    double userdata[2], buf[6];
    MPI_Status status;
    MPI_Request request;

    while(1)
    {
        MPI_Recv(buf, 6, MPI_DOUBLE, 0, 0, master_comm, &status);

        const double xi     = buf[0];
        const double LbyR   = buf[1];
        const double omegap = buf[2];
        const double gamma_ = buf[3];
        const int m         = buf[4];
        const int ldim      = buf[5];

        if(xi < 0)
            break;

        casimir_t *casimir = casimir_init(LbyR);
        if(!isinf(omegap))
        {
            userdata[0] = omegap;
            userdata[1] = gamma_;
            casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);
        }

        casimir_set_ldim(casimir, ldim);
        double logdet = casimir_logdetD(casimir, xi, m);
        TERMINATE(isnan(logdet), "LbyR=%.10g, xi=%.10g, m=%d, ldim=%d", LbyR, xi, m, ldim);
        casimir_free(casimir);

        MPI_Isend(&logdet, 1, MPI_DOUBLE, 0, 0, master_comm, &request);
        MPI_Wait(&request, &status);
    }

    return 0;
}

void usage(FILE *stream)
{
    fprintf(stream,
"Usage: casimir_T0 [OPTIONS]\n\n"
"This program will calculate the free Casimir energy F(T=0,L/R) for the\n"
"plane-sphere geometry for given L/R and temperature T. The output is in\n"
"units of ħc/(L+R).\n"
"\n"
"This program uses MPI for parallization. The program needs at least two\n"
"cores to run.\n"
"\n"
"The free eenergy at T=0 is calculated using integration:\n"
"   F(L/R) = \\int_0^\\infty dξ logdet D(ξ)\n"
"The integration is done using Gauss-Laguerre integration or adaptive\n"
"Gauß-Kronrod quadrature. The integrand decays exponentially, c.f.\n"
"Eq. (6.33) of [1].\n"
"\n"
"References:\n"
"  [1] Hartmann, Negative Casimir entropies in the plane-sphere geometry, 2014\n"
"\n"
"Mandatory options:\n"
"    -x, --LbyR L/R\n"
"        Separation L between sphere and plane divided by radius of sphere,\n"
"        where L/R > 0.\n"
"\n"
"Further options:\n"
"    -T, --temperature TEMPERATURE\n"
"        Set temperature to TEMPERATURE (default: 0)\n"
"\n"
"    -L LDIM\n"
"        Set ldim to the value LDIM. (default: ldim=ceil(%g*R/L))\n"
"\n"
"    -c, --cutoff CUTOFF\n"
"        Stop summation over m for a given value of ξ if\n"
"            logdet(D^m(ξ))/logdet(D^0(ξ) < CUTOFF\n"
"        (default: %g)\n"
"\n"
"    -e, --epsrel\n"
"       Request relative accuracy of EPSREL for Gauß-Kronrod integration\n"
"       (default: %g)\n"
"\n"
"    -N, --order\n"
"        Order of Gauss-Laguerre integration. If not specified, use adaptive\n"
"        Gauß-Kronrod quadrature.\n"
"\n"
"    --omegap OMEGAP\n"
"        Model the metals using the Drude/Plasma model and set Plasma\n"
"        frequency to OMEGAP. (DEFAULT: perfect conductors)\n"
"\n"
"    --gamma GAMMA\n"
"        Set dissipation of Drude model to GAMMA. Ignored if no value for\n"
"        --omegap is given. (DEFAULT: perfect conductors)\n"
"\n"
"    -z, --zero\n"
"        Also compute the term for xi=0; does only work for PC or Drude\n"
"\n"
"    -v, --verbose\n"
"        Also print results for each m\n"
"\n"
"    -h, --help\n"
"        Show this help\n",
    ETA, EPSREL, CUTOFF);
}

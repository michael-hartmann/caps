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
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

#define PRECISION 1e-12
#define ORDER 50
#define LFAC 6.
#define IDLE 25

#define STATE_RUNNING 1
#define STATE_IDLE    0

void casimir_mpi_init(casimir_mpi_t *self, double LbyR, int ldim, double precision, int cores)
{
    self->LbyR      = LbyR;
    self->ldim      = ldim;
    self->precision = precision;
    self->cores     = cores;
    self->tasks     = xmalloc(cores*sizeof(casimir_task_t *));

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
            double buf[] = { xi, self->LbyR, m, self->ldim };

            task->index = index;
            task->xi    = xi;
            task->m     = m;
            task->state = STATE_RUNNING;

            MPI_Send (buf,        4, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
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

double integrand(double xi, casimir_mpi_t *casimir_mpi)
{
    int m;
    double terms[10000] = { NAN };
    bool debug = false;
    const double nmax = sizeof(terms)/sizeof(double);
    const double precision = casimir_mpi->precision;

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

                if(debug)
                    fprintf(stderr, "# m=%d, xi=%.12g, logdetD=%.12g\n", task->m, xi, task->value);

                if(!isnan(terms[0]) && v/terms[0] < precision)
                    goto done;
            }

            usleep(IDLE);
        }
    }

    done:

    /* retrieve all remaining running jobs */
    while(casimir_mpi_get_running(casimir_mpi) > 0)
    {
        casimir_task_t *task = NULL;

        while(casimir_mpi_retrieve(casimir_mpi, &task))
        {
            terms[task->m] = task->value;
            if(debug)
                fprintf(stderr, "# m=%d, xi=%.12g, logdetD=%.12g\n", task->m, xi, task->value);
        }

        usleep(IDLE);
    }

    terms[0] /= 2; /* m = 0 */
    return kahan_sum(terms, m);
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
    int order = ORDER, ldim = 0, ret = 0;
    double LbyR = -1, lfac = LFAC, precision = PRECISION;
    casimir_mpi_t casimir_mpi;

    #define EXIT(n) do { ret = n; goto out; } while(0)

    /* parse command line options */
    while (1)
    {
        int c;
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },
            { "LbyR",      required_argument, 0, 'x' },
            { "ldim",      required_argument, 0, 'L' },
            { "order",     required_argument, 0, 'N' },
            { "lscale",    required_argument, 0, 'l' },
            { "precision", required_argument, 0, 'p' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "x:L:N:l:c:p:dh", long_options, &option_index);

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
            case 'L':
                ldim = atoi(optarg);
                break;
            case 'l':
                lfac = atof(optarg);
                break;
            case 'p':
                precision = atof(optarg);
                break;
            case 'N':
                order = atoi(optarg);
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
    if(order < 1)
    {
        fprintf(stderr, "order must be positive.\n\n");
        usage(stderr);
        EXIT(1);
    }
    if(precision <= 0)
    {
        fprintf(stderr, "precision must be positive.\n\n");
        usage(stderr);
        EXIT(1);
    }
    if(ldim <= 0)
    {
        if(lfac > 0)
        {
            ldim = lfac/LbyR;
        }
        else if(lfac < 0)
        {
            fprintf(stderr, "lfac must be positive\n\n");
            usage(stderr);
            EXIT(1);
        }
        else
        {
            fprintf(stderr, "ldim must be positive\n\n");
            usage(stderr);
            EXIT(1);
        }
    }

    if(cores < 2)
    {
        fprintf(stderr, "need at least 2 cores\n");
        EXIT(1);
    }

    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    printf("# LbyR  = %.15g\n", LbyR);
    printf("# order = %d\n", order);
    printf("# prec  = %g\n", precision);
    printf("# ldim  = %d\n", ldim);
    printf("# cores = %d\n", cores);
    printf("#\n");

    casimir_mpi_init(&casimir_mpi, LbyR, ldim, precision, cores);

    if(true)
    {
        double *xk, *ln_wk;
        /* estimate, cf. eq. (6.33) */
        double alpha = 2*LbyR/(1+LbyR);
        double integral = 0;
        order = gausslaguerre_nodes_weights(order, &xk, &ln_wk);

        for(int k = 0; k < order; k++)
        {
            double t0 = now();
            double xi = xk[k]/alpha;
            double v = integrand(xi, &casimir_mpi);
            printf("# k=%d, xi=%g, logdetD=%.15g, t=%g\n", k, xi, v, now()-t0);

            integral += exp(ln_wk[k]+xk[k])*v;
        }

        /* free energy for T=0 */
        double F0 = integral/alpha/M_PI;

        printf("#\n");
        printf("# L/R, ldim, order, alpha, F(T=0)*(L+R)/(ħc)\n");
        printf("%g, %d, %d, %.15g, %.12g\n", LbyR, ldim, order, alpha, F0);
    }
    else
    {

    }


    casimir_mpi_free(&casimir_mpi);
out:

    return ret;
}

int slave(MPI_Comm master_comm, int rank)
{
    double buf[5];
    MPI_Status status;
    MPI_Request request;

    while(1)
    {
        MPI_Recv(buf, 5, MPI_DOUBLE, 0, 0, master_comm, &status);

        const double xi   = buf[0];
        const double LbyR = buf[1];
        const int m    = (int)buf[2];
        const int ldim = buf[3];

        if(xi < 0)
            break;

        casimir_t *casimir = casimir_init(LbyR);
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
    fprintf(stderr,
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
"The integration is done using Gauss-Laguerre integration. The integrand\n"
"decays exponentially, c.f. eq. (6.33) of [1].\n"
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
"    -l, --lscale\n"
"        Specify parameter lscale. The vector space has to be truncated for\n"
"        some value ldim. This program will use ldim=(R/L*lscale) (default: %g)\n"
"\n"
"    -L LDIM\n"
"        Set ldim to the value LMAX. When -L is specified, -l will be ignored\n"
"\n"
"    -p, --precision\n"
"        Set precision to given value (default: %g)\n"
"\n"
"    -N, --order\n"
"        Order of Gauss-Laguerre integration (default: %d)\n"
"\n"
"    -h,--help\n"
"        Show this help\n",
    LFAC, PRECISION, ORDER);
}

#define _BSD_SOURCE /* make usleep work */

#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#include "casimir_T0.h"

#include "integration_perf.h"
#include "gausslaguerre.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

#define PRECISION 1e-12
#define ORDER 50
#define LFAC 6.
#define IDLE 10

#define STATE_RUNNING 1
#define STATE_IDLE    0

void casimir_mpi_init(casimir_mpi_t *self, double LbyR, int lmax, double precision, int cores)
{
    self->LbyR      = LbyR;
    self->lmax      = lmax;
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
            double buf[] = { xi, self->LbyR, m, self->lmax };

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
    int order = ORDER, lmax = 0, ret = 0;
    double F0;
    double alpha, LbyR = -1, lfac = LFAC, precision = PRECISION;
    double integral = 0, *xk, *ln_wk;
    double **values = NULL;
    casimir_mpi_t casimir_mpi;

    #define EXIT(n) do { ret = n; goto out; } while(0)

    /* parse command line options */
    while (1)
    {
        int c;
        struct option long_options[] = {
            { "help",      no_argument,       0, 'h' },
            { "LbyR",      required_argument, 0, 'x' },
            { "lmax",      required_argument, 0, 'L' },
            { "order",     required_argument, 0, 'N' },
            { "lscale",    required_argument, 0, 'l' },
            { "precision", required_argument, 0, 'p' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "x:L:N:l:c:p:h", long_options, &option_index);

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
                lmax = atoi(optarg);
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
    if(lmax <= 0)
    {
        if(lfac > 0)
        {
            lmax = lfac/LbyR;
        }
        else if(lfac < 0)
        {
            fprintf(stderr, "lfac must be positive\n\n");
            usage(stderr);
            EXIT(1);
        }
        else
        {
            fprintf(stderr, "lmax must be positive\n\n");
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

    /* estimate, cf. eq. (6.33) */
    alpha = 2*LbyR/(1+LbyR);
    order = gausslaguerre_nodes_weights(order, &xk, &ln_wk);

    printf("# LbyR  = %.15g\n", LbyR);
    printf("# alpha = %.15g\n", alpha);
    printf("# order = %d\n", order);
    printf("# prec  = %g\n", precision);
    printf("# lmax  = %d\n", lmax);
    printf("# cores = %d\n", cores);
    printf("#\n");

    casimir_mpi_init(&casimir_mpi, LbyR, lmax, precision, cores);

    /* allocate memory */
    values = xmalloc(order*sizeof(double *));
    for(int i = 0; i < order; i++)
    {
        values[i] = xmalloc(lmax*sizeof(double));

        for(int m = 0; m < lmax; m++)
            values[i][m] = NAN;
    }

    /* gather all data */
    for(int m = 0; m < lmax; m++)
    {
        for(int i = 0; i < order; i++)
        {
            int k = 0;

            /* k is either 0 or the last value in the list that is non NAN */
            while(k < (m-1) && !isnan(values[i][k+1]))
                k++;

            /* check if we still have to calculate this */
            if(k <= 0 || values[i][k]/values[i][0] >= precision)
                while(1)
                {
                    casimir_task_t *task = NULL;

                    /* send job */
                    if(casimir_mpi_submit(&casimir_mpi, i, xk[i]/alpha, m))
                        break;

                    /* retrieve jobs */
                    while(casimir_mpi_retrieve(&casimir_mpi, &task))
                        values[task->index][task->m] = task->value;

                    usleep(IDLE);
                }
        }

        printf("# m = %d\n", m);
    }

    /* retrieve all running jobs */
    while(casimir_mpi_get_running(&casimir_mpi) > 0)
    {
        casimir_task_t *task = NULL;

        while(casimir_mpi_retrieve(&casimir_mpi, &task))
            values[task->index][task->m] = task->value;

        usleep(IDLE);
    }

    casimir_mpi_free(&casimir_mpi);

    /* do Gauss-Legendre quadrature */
    for(int i = 0; i < order; i++)
    {
        double value = 0;
        for(int m = 0; m < lmax; m++)
        {
            if(isnan(values[i][m]))
                break;

            if(m == 0)
                value = values[i][0]/2;
            else
                value += values[i][m];
        }
        printf("# k=%d, x=%.10Lg, logdetD(xi = x/alpha)=%.16Lg\n", i, xk[i], value);
        integral += exp(ln_wk[i]+xk[i])*value;
    }

    /* free energy for T=0 */
    F0 = integral/alpha/M_PI;

    printf("#\n");
    printf("# L/R, lmax, order, alpha, F(T=0)*(L+R)/(ħc)\n");
    printf("%g, %d, %d, %.15g, %.12Lg\n", LbyR, lmax, order, alpha, F0);

    /* free memory */
    for(int i = 0; i < order; i++)
    {
        xfree(values[i]);
        values[i] = NULL;
    }
    xfree(values);
    values = NULL;

out:
    return ret;
}

int slave(MPI_Comm master_comm, int rank)
{
    int m, lmax;
    double buf[5];
    double logdet,xi,LbyR;
    MPI_Status status;
    MPI_Request request;
    casimir_t casimir;

    while(1)
    {
        MPI_Recv(buf, 5, MPI_DOUBLE, 0, 0, master_comm, &status);

        xi   = buf[0];
        LbyR = buf[1];
        m    = (int)buf[2];
        lmax = buf[3];

        if(xi < 0)
            break;

        casimir_init(&casimir, LbyR, xi);
        casimir_set_lmax(&casimir, lmax);
        logdet = casimir_logdetD(&casimir, 1, m);
        casimir_free(&casimir);

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
"        some value lmax. This program will use lmax=(R/L*lscale) (default: %g)\n"
"\n"
"    -L LMAX\n"
"        Set lmax to the value LMAX. When -L is specified, -l will be ignored\n"
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

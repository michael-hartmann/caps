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
#include "edouble.h"
#include "utils.h"

#define PRECISION 1e-12
#define ORDER 50
#define LFAC 6.
#define IDLE_MASTER 5000
#define IDLE_SLAVE 5000

#define STATE_RUNNING 1
#define STATE_IDLE    0

/* calculate integrand logdetD(xi) = \sum_{m=0}^{m=lmax} logdetD^m(xi) */
double integrand(double xi, double LbyR, int lmax, double precision)
{
    casimir_t casimir;
    double v = 0, v0 = 0;
    int m = 0, use_trace = 0;

    /* initialize Casimir object and set lmax */
    casimir_init(&casimir, LbyR, xi);
    casimir_set_lmax(&casimir, lmax);

    /* sum up all contributions from m=0 to m=lmax */
    while(m <= lmax)
    {
        double v_m;

        if(use_trace)
            v_m = -casimir_trM(&casimir, 1, m);
        else
        {
            v_m = casimir_logdetD(&casimir, 1, m);

            /* For large arguments of xi the calculation of logdetD^m becomes
             * inaccurate. This corresponds to large distances (large xi <->
             * large distances), so we can use the approximation
             *      logdetD^m \approx -Tr M^m(xi).
             *
             * The factor 1e-8 was determined by experience.
             */
            if(fabs(v_m) < 1e-8)
            {
                v_m = -casimir_trM(&casimir, 1, m);
                use_trace = 1;
            }
        }

        /* divide contribution from m=0 by 2 */
        if(m == 0)
        {
            v0 = v_m;
            v_m /= 2;
        }

        v += v_m;

        /* if contribution to the sum is smaller than v_m/v0, we're done */
        if(fabs(v_m/v0) < precision)
            break;

        m++;
    }

    casimir_free(&casimir);

    return v;
}


void casimir_mpi_init(casimir_mpi_t *self, double LbyR, int lmax, double precision, int cores)
{
    int i;

    self->LbyR      = LbyR;
    self->lmax      = lmax;
    self->precision = precision;
    self->cores     = cores;
    self->tasks     = xmalloc(cores*sizeof(casimir_task_t *));

    for(i = 0; i < cores; i++)
    {
        casimir_task_t *task = xmalloc(sizeof(casimir_task_t));
        task->index = -1;
        task->state = STATE_IDLE;
        self->tasks[i] = task;
    }
}

void casimir_mpi_free(casimir_mpi_t *self)
{
    int i;
    double buf[] = { -1, -1, -1, -1, -1 };

    for(i = 1; i < self->cores; i++)
        MPI_Send(buf, 5, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

    for(i = 0; i < self->cores; i++)
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
    int i;

    for(i = 1; i < self->cores; i++)
        if(self->tasks[i]->state == STATE_RUNNING)
            running++;

    return running;
}

int casimir_mpi_submit(casimir_mpi_t *self, int index, double xi)
{
    int i;

    for(i = 1; i < self->cores; i++)
    {
        casimir_task_t *task = self->tasks[i];

        if(task->state == STATE_IDLE)
        {
            double buf[] = { index, xi, self->LbyR, self->lmax, self->precision };

            task->index = index;
            task->xi    = xi;
            task->state = STATE_RUNNING;

            MPI_Send (buf,        5, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Irecv(task->recv, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &task->request);

            return 1;
        }
    }

    return 0;
}


int casimir_mpi_retrieve(casimir_mpi_t *self, casimir_task_t **task_out)
{
    int i;

    *task_out = NULL;

    for(i = 1; i < self->cores; i++)
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
    int i, order = ORDER, lmax = 0, ret = 0;
    edouble F0;
    double alpha, LbyR = -1, lfac = LFAC, precision = PRECISION;
    edouble integral = 0, *xk, *ln_wk;
    double values[order];
    casimir_mpi_t casimir_mpi;

    #define EXIT(n) do { ret = n; goto out; } while(0)

    /* parse command line options */
    while (1)
    {
        int c;
        struct option long_options[] =
        {
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

    /* do Gauss-Legendre quadrature */
    i = 0;
    while(i < order)
    {
        casimir_task_t *task = NULL;

        /* send jobs as long we have work to do */
        while(i < order && casimir_mpi_submit(&casimir_mpi, i, xk[i]/alpha))
            i++;

        while(1)
        {
            int flag = casimir_mpi_retrieve(&casimir_mpi, &task);

            if(flag)
            {
                values[task->index] = task->value;
                printf("# k=%d, x=%.15Lg, logdetD(xi = x/alpha)=%.15g\n", task->index, xk[task->index], values[task->index]);
            }
            else if(i != order || casimir_mpi_get_running(&casimir_mpi) == 0)
                break;
        }

        usleep(IDLE_MASTER);
    }

    casimir_mpi_free(&casimir_mpi);

    for(i = 0; i < order; i++)
        integral += expe(ln_wk[i]+xk[i])*values[i];

    /* free energy or T=0 */
    F0 = integral/alpha/M_PI;

    printf("#\n");
    printf("# L/R, lmax, order, alpha, F(T=0)\n");
    printf("%.15g, %d, %d, %.15g, %.15Lg\n", LbyR, lmax, order, alpha, F0);

out:
    return ret;
}

int slave(MPI_Comm master_comm, int rank)
{
    int k, lmax;
    double buf[5];
    double v,xi,precision,LbyR;
    MPI_Status status;
    MPI_Request request;

    while(1)
    {
        MPI_Recv(buf, 5, MPI_DOUBLE, 0, 0, master_comm, &status);

        k         = buf[0];
        xi        = buf[1];
        LbyR      = buf[2];
        lmax      = buf[3];
        precision = buf[4];

        if(k < 0)
            break;

        v = integrand(xi, LbyR, lmax, precision);

        buf[0] = v;

        MPI_Isend(buf, 1, MPI_DOUBLE, 0, 0, master_comm, &request);
        MPI_Wait(&request, &status);

        usleep(IDLE_SLAVE);
    }

    return 0;
}

void usage(FILE *stream)
{
    fprintf(stderr,
"Usage: casimir_T0 [OPTIONS]\n"
"This program will calculate the free Casimir energy F(T=0,L/R) for the\n"
"plane-sphere geometry for given L/R and temperature T. The output is in scaled\n"
"units.\n"
"\n"
"This program uses MPI for parallization. The program needs at least two cores\n"
"to run.\n"
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
"        Order of Gauss-Laguerre integrateion (default: %d)\n"
"\n"
"    -h,--help\n"
"        Show this help\n",
    LFAC, PRECISION, ORDER);
}

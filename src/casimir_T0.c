#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"

#include "casimir_T0.h"

#include "integration_perf.h"
#include "gausslaguerre.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "edouble.h"


#define PRECISION 1e-12
#define ORDER 50
#define LFAC 6.

int master(int argc, char *argv[], int cores);

double integrand(double xi, double LbyR, int lmax, double precision)
{
    casimir_t casimir;
    integration_perf_t int_perf;
    double v = 0, v0 = 0;
    int m = 0, use_trace = 0;
    const double Q = 1/(1+LbyR);

    casimir_init(&casimir, Q, xi);
    casimir_set_lmax(&casimir, lmax);

    casimir_integrate_perf_init(&int_perf, 1*casimir.T, casimir.lmax);

    while(m <= lmax)
    {
        double v_m;

        if(use_trace)
            v_m = -casimir_trM(&casimir, 1, m, &int_perf);
        else
        {
            v_m = casimir_logdetD(&casimir, 1, m, &int_perf);
            if(fabs(v_m) < 1e-8)
            {
                v_m = -casimir_trM(&casimir, 1, m, &int_perf);
                use_trace = 1;
            }
        }

        if(m == 0)
        {
            v0 = v_m;
            v_m /= 2;
        }

        v += v_m;

        if(fabs(v_m/v0) < precision)
            break;

        m++;
    }


    casimir_integrate_perf_free(&int_perf);
    casimir_free(&casimir);

    return v;
}


int submit_job(int destination, MPI_Request *request, double *recv, int k, double xi, double LbyR, int lmax, double precision)
{
    double buf[] = { k, xi, LbyR, lmax, precision };

    MPI_Send(buf, 5, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD);
    MPI_Irecv(recv, 2, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD, request);

    return 0;
}

int retrieve_job(MPI_Request *request, double *buf, int *index, double *value)
{
    MPI_Status status;
    MPI_Wait(request, &status);

    *index = buf[0];
    *value = buf[1];

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
    double F0, alpha, LbyR = -1, lfac = LFAC, precision = PRECISION;
    edouble integral = 0, *xk, *ln_wk;
    int k = 0;
    int i;
    double recv[cores][2];
    MPI_Request requests[cores];

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
    
    if(LbyR < 0)
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

    for(i = 1; i < MIN(order,cores); i++)
    {
        submit_job(i, &requests[i], recv[i], k, xk[k]/alpha, LbyR, lmax, precision);
        k++;
    }

    /* do Gauss-Legendre quadrature */
    while(1)
    {
        int index;
        double value;

        for(i = 1; i < cores && k < order; i++)
        {
            int flag = 0;
            MPI_Status status;

            MPI_Test(&requests[i], &flag, &status);
            if(flag)
            {
                retrieve_job(&requests[i], recv[i], &index, &value);

                printf("# k=%d, x=%.15g, logdetD(xi = x/alpha)=%.15g\n", index, (double)xk[index], value);
                integral += expe(ln_wk[index]+xk[index])*value;

                submit_job(i, &requests[i], recv[i], k, xk[k]/alpha, LbyR, lmax, precision);
                k++;
            }

            usleep(5*1000);
        }

        if(k >= order)
        {

            for(i = 1; i < cores; i++)
            {
                retrieve_job(&requests[i], recv[i], &index, &value);

                printf("# k=%d, x=%.15g, logdetD(xi = x/alpha)=%.15g\n", index, (double)xk[index], value);
                integral += expe(ln_wk[index]+xk[index])*value;
            }

            break;
        }
    }

    /* free energy or T=0 */
    F0 = (double)(integral/alpha/M_PI);
    
    printf("#\n");
    printf("# L/R, order, alpha, F(T=0)\n");
    printf("%.15g, %d, %.15g, %.15g\n", LbyR, order, alpha, F0);

out:
    for(i = 1; i < cores; i++)
    {
        double buf[] = { -1, -1, -1, -1, -1 };
        MPI_Send(buf, 5, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }

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

        buf[0] = k;
        buf[1] = v;

        MPI_Isend(buf, 2, MPI_DOUBLE, 0, 0, master_comm, &request);
        MPI_Wait(&request, &status);
    }

    return 0;
}

void usage(FILE *stream) {
    fprintf(stderr,
"Usage: casimir_T0 [OPTIONS]\n"
"This program will calculate the free Casimir energy F(T=0,L/R) for the\n"
"plane-sphere geometry for given L/R and temperature T. The output is in scaled\n"
"units.\n"
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

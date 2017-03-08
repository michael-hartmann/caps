#ifndef CASIMIR_T0
#define CASIMIR_T0

#include <stdlib.h>

typedef struct {
    int index, m;
    double xi;
    double recv[1];
    double value;
    MPI_Request request;
    int state;
} casimir_task_t;

typedef struct {
    double LbyR, omegap, gamma, precision, alpha;
    int k, ldim, cores;
    casimir_task_t **tasks;
} casimir_mpi_t;

void casimir_mpi_init(casimir_mpi_t *self, double LbyR, double omegap, double gamma_, int ldim, double precision, int cores);
void casimir_mpi_free(casimir_mpi_t *self);
int casimir_mpi_submit(casimir_mpi_t *self, int index, double xi, int m);
int casimir_mpi_retrieve(casimir_mpi_t *self, casimir_task_t **task_out);
int casimir_mpi_get_running(casimir_mpi_t *self);

void usage(FILE *stream);

double integrand(double xi, casimir_mpi_t *casimir_mpi);
int master(int argc, char *argv[], int cores);
int slave(MPI_Comm master_comm, int rank);

void stop_process(int task);
int submit_job(int process, MPI_Request *request, double *recv, int k, double xi, double LbyR, int ldim, double precision);
int retrieve_job(MPI_Request *request, double *buf, int *index, double *value);

#endif

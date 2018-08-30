#ifndef CASIMIR_T0
#define CASIMIR_T0

#include <stdbool.h>
#include <stdlib.h>

typedef struct {
    int index, m;
    double xi;
    double recv;
    double value;
    MPI_Request request;
    int state;
} casimir_task_t;

typedef struct {
    double L, R, T, omegap, gamma, cutoff, iepsrel, alpha;
    int ldim, cores;
    bool verbose;
    casimir_task_t **tasks;
    int determinants;
    char filename[512];
} casimir_mpi_t;

casimir_mpi_t *casimir_mpi_init(double L, double R, double T, char *filename, double omegap, double gamma_, int ldim, double cutoff, double iepsrel, int cores, bool verbose);
void casimir_mpi_free(casimir_mpi_t *self);
int casimir_mpi_submit(casimir_mpi_t *self, int index, double xi, int m);
int casimir_mpi_retrieve(casimir_mpi_t *self, casimir_task_t **task_out);
int casimir_mpi_get_running(casimir_mpi_t *self);
int casimir_get_determinants(casimir_mpi_t *self);

void usage(FILE *stream);

void F_HT(casimir_mpi_t *casimir_mpi, double omegap, double *drude, double *pr, double *plasma);
double F_xi(double xi, casimir_mpi_t *casimir_mpi);
void master(int argc, char *argv[], int cores);
void slave(MPI_Comm master_comm, int rank);

void stop_process(int task);
int submit_job(int process, MPI_Request *request, double *recv, int k, double xi, double LbyR, int ldim, double cutoff);
int retrieve_job(MPI_Request *request, double *buf, int *index, double *value);

#endif

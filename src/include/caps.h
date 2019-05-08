#ifndef CAPS_H
#define CAPS_H

#include <stdbool.h>

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

caps_mpi_t *caps_mpi_init(double L, double R, double T, char *filename, char *resume, double omegap, double gamma_, int ldim, double cutoff, double iepsrel, int cores, bool verbose);
void caps_mpi_free(caps_mpi_t *self);
int caps_mpi_submit(caps_mpi_t *self, int index, double xi, int m);
int caps_mpi_retrieve(caps_mpi_t *self, caps_task_t **task_out);
int caps_mpi_get_running(caps_mpi_t *self);
int caps_get_determinants(caps_mpi_t *self);

void usage(FILE *stream);

void F_HT(caps_mpi_t *caps_mpi, double omegap, double *drude, double *pr, double *plasma);
double F_xi(double xi, caps_mpi_t *caps_mpi);
void master(int argc, char *argv[], int cores);
void slave(MPI_Comm master_comm, int rank);

void stop_process(int task);
int submit_job(int process, MPI_Request *request, double *recv, int k, double xi, double LbyR, int ldim, double cutoff);
int retrieve_job(MPI_Request *request, double *buf, int *index, double *value);

#endif

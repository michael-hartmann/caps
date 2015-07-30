#ifndef CASIMIR_T0
#define CASIMIR_T0

#include <stdlib.h>

void usage(FILE *stream);
double integrand(double xi, double LbyR, int lmax, double precision);

int master(int argc, char *argv[], int cores);
int slave(MPI_Comm master_comm, int rank);

int submit_job(int destination, MPI_Request *request, double *recv, int k, double xi, double LbyR, int lmax, double precision);
int retrieve_job(MPI_Request *request, double *buf, int *index, double *value);

#endif

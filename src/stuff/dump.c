#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        fprintf(stderr, "Usage: %s filename\n", argv[0]);
        return 1;
    }

    matrix_float80 *M = matrix_float80_load(argv[1]);
    if(M == NULL)
    {
        fprintf(stderr, "Couldn't load matrix\n");
        return 1;
    }

    int dim = M->dim;

    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            printf("%d, %d, %Lg\n", i,j, matrix_get(M,i,j));

    matrix_float80_free(M);

    return 0;
}

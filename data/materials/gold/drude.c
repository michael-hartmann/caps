#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double omegap1, gamma1;
    double omegap2, gamma2;
    double xi_min, xi_max;
    size_t points;
    double *xi, *epsm1;
} metal_t;


void metal_free(metal_t *metal)
{
    free(metal->xi);
    free(metal->epsm1);
    free(metal);
}

metal_t *metal_init(const char *filename, double omegap1, double gamma1, double omegap2, double gamma2)
{
    size_t size = 128, points = 0;
    char line[2048];
    FILE *f = fopen(filename, "r");

    if(f == NULL)
        return NULL;

    metal_t *metal = malloc(sizeof(metal_t));
    metal->omegap1 = omegap1;
    metal->gamma1  = gamma1;
    metal->omegap2 = omegap2;
    metal->gamma2  = gamma2;

    metal->xi      = malloc(size*sizeof(double));
    metal->epsm1   = malloc(size*sizeof(double));

    while(fgets(line, sizeof(line)/sizeof(char), f) != NULL)
    {
        if(line[0] == '#')
            continue;

        char *p = strchr(line, ' ');
        if(p != NULL)
        {
            double xi    = atof(line);
            double epsm1 = atof(p)-1;

            /* increase buffer */
            if(points >= size)
            {
                size *= 2;
                metal->xi    = realloc(metal->xi,    size*sizeof(double));
                metal->epsm1 = realloc(metal->epsm1, size*sizeof(double));
            }

            metal->xi[points]    = xi;
            metal->epsm1[points] = epsm1;

            if(points > 0 && xi <= metal->xi[points-1])
            {
                metal_free(metal);
                return NULL;
            }


            points++;
        }
    }

    metal->xi_min = metal->xi[0];
    metal->xi_max = metal->xi[points-1];
    metal->points = points;

    return metal;
}

double epsilonm1(double xi, void *args)
{
    metal_t *metal = (metal_t *)args;
    const double xi_min = metal->xi_min, xi_max = metal->xi_max;

    if(xi < xi_min)
    {
        /* use drude */
        return pow(metal->omegap1,2)/(xi*(xi+metal->gamma1));
    }
    else if(xi > xi_max)
    {
        /* use foo */
        return pow(metal->omegap2,2)/(xi*(xi+metal->gamma2));
    }

    /* do binary search */
    int lower = 0, upper = metal->points-1;

    /* find indices of upper and lower value */
    while(upper-lower != 1)
    {
        int middle = (upper+lower)/2;
        double xi_middle = metal->xi[middle];
        if(xi_middle > xi)
            upper = middle;
        else
            lower = middle;
    }

    const double xi_lower = metal->xi[lower], xi_upper = metal->xi[upper];
    const double epsm1_lower = metal->epsm1[lower], epsm1_upper = metal->epsm1[upper];

    /* linear interpolation */
    const double epsm1 = epsm1_lower + (xi-xi_lower)*(epsm1_upper-epsm1_lower)/(xi_upper-xi_lower);

    //printf("xi_min = %g, xi_max = %g, points = %zu, lower=%g, upper=%g => %g\n", xi_min, xi_max, points, metal->xi[lower], metal->xi[upper], epsm1);
    return epsm1;
}


int main(int argc, char *argv[])
{
    double xi = atof(argv[1]);

    metal_t *metal = metal_init("GoldEpsIm.dat", 0, 0, 0, 0);
    double epsm1 = epsilonm1(xi,metal);
    metal_free(metal);

    printf("%g %g\n", xi, epsm1);

    return 0;
}

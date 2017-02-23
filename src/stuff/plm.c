#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define M_PI 3.141592653589793
#define R0 25

double Clm(int l, int m)
{
    return exp(2*lgamma(m+0.5)+lgamma(l+1) - log(M_PI) - m*log(2) -lgamma(l+m+1.5) -lgamma(m+1));
}

/* Evaluation of Pl(cosh(xi)) using an asymptotic expansion provided that
 * (l+1)*sinh(xi) >= 25. */
double Pl(int l, double x)
{
    const double xi = acosh(x);
    const double expxi = exp(xi);
    const double sinhxi = sqrt((x+1)*(x-1));

    printf("(l+1)*sinh(xi) >= R0 = %s\n", ((l+1)*sinhxi >= R0) ? "true" : "false");

    double sum = 0;
    double sinhxi_m = 1; /* sinh(xi)**m */
    double expxi_m  = 1; /* exp(xi)**m */
    for(int m = 0; m < 18; m++)
    {
        //printf("m=%d, Clm(l=%d,m=%d)=%g, cosh=%g, sinhxi=%g\n", m, l,m,Clm(l,m), (exp(m*xi)+exp(-(m+2*l+1)*xi)), sinhxi_m);

        sum += Clm(l,m) * (expxi_m+exp(-(m+2*l+1)*xi)) /sinhxi_m;
        sinhxi_m *= sinhxi;
        expxi_m  *= expxi;
    }

    return log(2/(M_PI*sinhxi))/2 + log(sum) + (l+0.5)*xi - log(2);
}

double dPl(int l, double x)
{
    double lnPl  = Pl(l,   x);
    double lnPlm = Pl(l-1, x);

    return log((l*x)/((x+1)*(x-1))) + lnPl + log1p( -exp(lnPlm-lnPl)/x );
}

double Plm(int l, int m, double x)
{
    double prefactor, v0, v1;
    double root = 1/sqrt((x+1)*(x-1));

    prefactor = Pl(l,x);
    v0 = 1;
    v1 = -exp(dPl(l,x)-prefactor)/root;

    for(int mm = 1; mm < m; mm++)
    {
        //printf("v0=%g, v1=%g\n", v0, v1);
        double v = (l+mm)*(l-mm+1)*v0 + 2*mm*x*v1*root;
        v0 = v1;
        v1 = v;

        if(fabs(v) > 1e150)
        {
            prefactor += log(1e150);
            v0 /= 1e150;
            v1 /= 1e150;
        }
    }

    return prefactor+log(fabs(v1));
}

int main(int argc, char *argv[])
{
    int l,m;
    double x;

    if(argc < 4)
    {
        fprintf(stderr, "%s l m x\n", argv[0]);
        return 1;
    }

    l = atoi(argv[1]);
    m = atoi(argv[2]);
    x = atof(argv[3]);

    double  v =  Pl(l,x);
    double dv = dPl(l,x);

    printf("l=%d, x=%g, Pl(x)=%.15g, Pl'(x)=%.15g\n", l, x, v, dv);
    printf("m=%d, Plm(x)=%.15g)\n", m, Plm(l,m,x));

    return 0;
}


/**
 * @name Proximity Force Approximation
 */
/*@{*/

/* integrand */
static double _integrand1(double q, void *args)
{
    double calLbyd, xi, k, rTE, rTM;

    casimir_pfa_t *pfa = (casimir_pfa_t *)args;
    calLbyd = pfa->calLbyd;
    xi = pfa->xi;
    casimir_t *casimir = pfa->casimir;

    k = sqrt(pow_2(calLbyd*q/2)-pow_2(xi)); /* k scaled */
    casimir->rp(casimir, xi, k, &rTE, &rTM);

    return q*(log1p(-pow_2(rTE)*exp(-q)) + log1p(-pow_2(rTM)*exp(-q)));
}

static double _integrand_xi(double z, void *args)
{
    int ier, neval;
    double integral, abserr;

    casimir_pfa_t *pfa = (casimir_pfa_t *)args;
    pfa->xi = z*pfa->calLbyd; /* xi scaled */

    integral = dqagi(_integrand1, 2*z, 1, 0, 1e-12, &abserr, &neval, &ier, args);

    WARN(ier != 0, "ier=%d, integral=%g, abserr=%g, abrel=%g, neval=%d", ier, integral, abserr, fabs(abserr/integral), neval);

    return integral;
}

/* integrand of proximity force approximation */
static double _casimir_pp(double t, void *args)
{
    casimir_pfa_t *pfa = (casimir_pfa_t *)args;
    const double calLbyd = (1+1/pfa->casimir->LbyR)/t;
    const double T = pfa->T;

    pfa->calLbyd = calLbyd;

    if(T == 0)
    {
        double abserr, integral;
        int ier, neval;

        integral = dqagi(_integrand_xi, 0, 1, 0, 1e-10, &abserr, &neval, &ier, args);

        WARN(ier != 0, "ier=%d, integral=%g, abserr=%g, abrel=%g, neval=%d", ier, integral, abserr, fabs(abserr/integral), neval);
        return integral/(16*pow_2(M_PI)*pow_2(t))*calLbyd;
    }
    else
    {
        /* do summation over Matsubara frequencies */
        const double v0 = _integrand_xi(0, args);
        double sum = v0/2;

        for(int n = 1;; n++)
        {
            const double v = _integrand_xi(n*T/calLbyd, args);
            sum += v;

            if(v/v0)
                return T/pow_2(4*M_PI)*sum/pow_2(t);
        }
    }
}

/* @brief Compute free energy using PFA
 *
 * Compute free energy using the proximity force approximation using
 *      F_PFA = 2*pi*R*L*Int_1^(1+R/L) dt F_pp(Lt)/A
 * where F_pp/A is the free energy per area for the plane-plane geometry.
 *
 * This function returns the free energy F in units of (L+R)/(hbar*c), i.e.,
 * it returns (L+R)/(hbar*c)*F.
 *
 * @param [in] casimir object
 * @param [in] T temperature
 * @param retval F free energy
 */
double casimir_pfa(casimir_t *casimir, double T)
{
    int ier, neval;
    double abserr, integral, LbyR = casimir->LbyR;

    casimir_pfa_t pfa =
    {
        .casimir = casimir,
        .T       = T
    };

    integral = dqags(_casimir_pp, 1, 1+1/LbyR, 0, 1e-10, &abserr, &neval, &ier, &pfa);

    return 2*M_PI/LbyR*integral;
}

/*@}*/

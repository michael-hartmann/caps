#include <math.h>

#include "libcasimir.h"
#include "integration_perf.h"
#include "unittest.h"
#include "sfunc.h"

#include "test_integration_perf.h"

static void _integrals(int l1, int l2, int m, double nT, casimir_integrals_t *cint);
static void _integrals(int l1, int l2, int m, double nT, casimir_integrals_t *cint)
{
    integration_perf_t int_perf;
    int lmax = MAX(l1,l2);
    casimir_integrate_perf_init(&int_perf, nT, lmax);
    casimir_integrate_perf(&int_perf, l1, l2, m, cint, NULL);
    casimir_integrate_perf_free(&int_perf);
}

int test_integration_perf(void)
{
    casimir_integrals_t cint;
    integration_perf_t int_perf;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    _integrals(1,1,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -2.980829253011726);
    AssertAlmostEqual(&test, cint.lnB_TM, -2.0645385211375711716721);
    AssertAlmostEqual(&test, cint.lnC_TM, -2.5753641449035618548786);
    AssertAlmostEqual(&test, cint.lnD_TM, -2.5753641449035618548786);

    _integrals(241,73,1,30,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 406.63047665158294437064);
    AssertAlmostEqual(&test, cint.lnB_TM, 419.71230683599700819362);
    AssertAlmostEqual(&test, cint.lnC_TM, 412.57255550309976814896);
    AssertAlmostEqual(&test, cint.lnD_TM, 413.766977852356385781);

    _integrals(241,1,1,30,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 249.75276347175786475423);
    AssertAlmostEqual(&test, cint.lnB_TM, 258.05248402595679167552);
    AssertAlmostEqual(&test, cint.lnC_TM, 251.17334248392289626321);
    AssertAlmostEqual(&test, cint.lnD_TM, 256.62788585958419530558);

    _integrals(241,241,1,30,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 838.84852861683729524124);
    AssertAlmostEqual(&test, cint.lnB_TM, 853.98316452183914507246);
    AssertAlmostEqual(&test, cint.lnC_TM, 846.41479992430049881808);

    _integrals(3,2,1,2,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -4.094372316589062);
    AssertAlmostEqual(&test, cint.lnB_TM, -1.970116759119433);
    AssertAlmostEqual(&test, cint.lnC_TM, -3.298725852652321);

    casimir_integrate_perf_init(&int_perf, 0.005, 4);
    casimir_integrate_perf(&int_perf, 4, 4, 0, &cint, NULL);
    AssertAlmostEqual(&test, cint.lnB_TM, 56.977025325953406);
    casimir_integrate_perf_free(&int_perf);

    _integrals(4,4,1,0.005,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +2.4806179125126554e17*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -2.2226323455151368e24*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.9457269656680333e20*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +6.9457269656680333e20*-2);

    _integrals(40,40,1,0.25,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +1.5754477603435539e159*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -6.3723632215476122e166*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -9.9568222699306801e162*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.9568222699306801e162*-2);

    _integrals(40,40,40,1,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +6.4140686579381969e91*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.0147301906459434e95*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -2.5352219594503741e93*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +2.5352219594503736e93*-2);

    _integrals(7,4,3,8.5,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +4.8180365200137397e-9*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.3731640166794149e-8*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.7659079909128738e-9*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.44463292099617e-9*-2);

    _integrals(40,40,0,2.5,&cint);
    casimir_integrate_perf(&int_perf, 40, 40, 0, &cint, NULL);

    _integrals(100,41,0,2.5,&cint);
    casimir_integrate_perf(&int_perf, 100, 41, 0, &cint, NULL);

    return test_results(&test, stderr);
}

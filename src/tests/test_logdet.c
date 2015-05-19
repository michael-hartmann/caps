#include "libcasimir.h"
#include "unittest.h"
#include "integration_drude.h"
#include "integration_perf.h"

#include "test_logdet.h"

int test_logdet(int cores)
{
    unittest_t test;
    casimir_t casimir;
    integration_perf_t int_perf;
    const double RbyScriptL = 0.97;
    const double T = 0.1;
    const int lmax = 200;
    double logdet;

    unittest_init(&test, "logdet M", "calculate logdet");

    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_lmax(&casimir, lmax);
    casimir_set_cores(&casimir, cores);


    logdet = casimir_logdetD(&casimir, 0, 0, NULL);
    AssertAlmostEqual(&test, logdet, -3.45236396285874);

    logdet = casimir_logdetD(&casimir, 0, 1, NULL);
    AssertAlmostEqual(&test, logdet, -2.63586999367158);

    logdet = casimir_logdetD(&casimir, 0, 10, NULL);
    AssertAlmostEqual(&test, logdet, -0.0276563864490425);

    casimir_integrate_perf_init(&int_perf, T, lmax);
    logdet = casimir_logdetD(&casimir, 1, 1, &int_perf);
    AssertAlmostEqual(&test, logdet, -2.63900987016801);
    casimir_integrate_perf_free(&int_perf);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}

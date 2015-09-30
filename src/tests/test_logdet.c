#include "libcasimir.h"
#include "unittest.h"
#include "integration_drude.h"
#include "integration_perf.h"

#include "test_logdet.h"

int test_logdet(void)
{
    unittest_t test;
    casimir_t casimir;
    integration_perf_t int_perf;
    double logdet;

    unittest_init(&test, "logdet M", "calculate logdet");

    casimir_init(&casimir, 0.006, 1.006);
    casimir_set_lmax(&casimir, 500);
    casimir_integrate_perf_init(&int_perf, 1.006, 500);

    logdet = casimir_logdetD(&casimir, 1, 1, &int_perf);
    AssertAlmostEqual(&test, logdet, -9.27424711542347);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}

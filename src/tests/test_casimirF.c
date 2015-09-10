#include "libcasimir.h"
#include "unittest.h"

#include "test_casimirF.h"

int test_casimirF(int cores)
{
    unittest_t test;
    casimir_t casimir;
    double F;

    unittest_init(&test, "casimirF", "Compare free energies");

    casimir_init(&casimir, 1/0.85-1, 2.7);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, cores);
    casimir_set_lmax(&casimir, 30);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -1.34361893570375);
    casimir_free(&casimir);

    casimir_init(&casimir, 1/0.7-1, 1);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, cores);
    casimir_set_lmax(&casimir, 15);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -0.220709222562969);
    casimir_free(&casimir);

    casimir_init(&casimir, 1/0.85-1, 2.7);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, cores);
    casimir_set_lmax(&casimir, 30);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -1.34361893570375);
    casimir_free(&casimir);

    return test_results(&test, stderr);
}

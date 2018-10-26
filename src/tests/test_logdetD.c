#include "libcasimir.h"
#include "misc.h"
#include "unittest.h"

#include "test_logdetD.h"

int test_logdetD()
{
    unittest_t test;
    casimir_t *casimir;

    /* xi=1 */
    double v1[] = {
        -7.284253322517682,   /* m=0 */
        -6.463971987842913,   /* m=1 */
        -4.590158706563189,   /* m=2 */
        -3.328160339577356,   /* m=3 */
        -2.440910567353524,   /* m=4 */
        -1.804043640922947,   /* m=5 */
        -1.340576287829648,   /* m=6 */
        -1.000032667657758,   /* m=7 */
        -0.7480828914339441,  /* m=8 */
        -0.5607500921873892,  /* m=9 */
        -0.4209565757132674,  /* m=10 */
        -0.3163606713158237,  /* m=11 */
        -0.2379466893689776,  /* m=12 */
        -0.1790756135473637,  /* m=13 */
        -0.1348293030091449,  /* m=14 */
        -0.1015481646873064,  /* m=15 */
        -0.07650001369497129, /* m=16 */
        -0.0576399686181587,  /* m=17 */
        -0.04343470369372667, /* m=18 */
        -0.03273287499638459, /* m=19 */
        -0.02466906185782241  /* m=20 */
    };

    /* xi=0.1 */
    double v2[] = {
        -7.744627647800944, /* m=0 */
        -6.437861787521742, /* m=1 */
        -4.603981083799886, /* m=2 */
        -3.336933059473016, /* m=3 */
        -2.445590617656273, /* m=4 */
        -1.806588305649722  /* m=5 */
    };

    /* xi=10 */
    double v3[] = {
        -4.993485411083642, /* m=0 */
        -4.987803023655219, /* m=1 */
        -3.763305479295732, /* m=2 */
        -2.847070956568539, /* m=3 */
        -2.154476036611317, /* m=4 */
        -1.629918416194979  /* m=5 */
    };

    unittest_init(&test, "casimir_logdetD", "Computation of determinant", 1e-10);

    casimir = casimir_init(0.01); /* L/R doesn't matter */
    casimir_set_ldim(casimir, 500);

    for(size_t m=0; m < sizeof(v1)/sizeof(double); m++)
        AssertAlmostEqual(&test, casimir_logdetD(casimir, 1, m), v1[m]);

    for(size_t m=0; m < sizeof(v2)/sizeof(double); m++)
        AssertAlmostEqual(&test, casimir_logdetD(casimir, 0.1, m), v2[m]);

    for(size_t m=0; m < sizeof(v3)/sizeof(double); m++)
        AssertAlmostEqual(&test, casimir_logdetD(casimir, 10, m), v3[m]);

    casimir_free(casimir);

    return test_results(&test, stderr);
}

int test_logdetD0()
{
    const double eps = 1e-13;
    unittest_t test;
    casimir_t *casimir;

    unittest_init(&test, "casimir_logdetD0", "Computation of determinant for xi=0", 1e-10);

    casimir = casimir_init(0.01); /* L/R doesn't matter */
    casimir_set_ldim(casimir, 500);

    AssertAlmostEqual(&test, casimir_logdetD0_plasma(casimir, 0.01, eps), -14.569722716816960073);
    AssertAlmostEqual(&test, casimir_logdetD0_plasma(casimir, 0.1, eps), -14.569723099773346675);
    AssertAlmostEqual(&test, casimir_logdetD0_plasma(casimir, 1, eps), -14.571258410941549499);
    AssertAlmostEqual(&test, casimir_logdetD0_plasma(casimir, 10, eps), -14.868792120441892024);
    AssertAlmostEqual(&test, casimir_logdetD0_plasma(casimir, 100, eps), -18.461022336857670467);
    AssertAlmostEqual(&test, casimir_logdetD0_plasma(casimir, 1000, eps), -25.441820432119136797);

    casimir_free(casimir);

    return test_results(&test, stderr);
}

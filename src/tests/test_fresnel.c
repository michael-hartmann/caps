#include "sfunc.h"
#include "unittest.h"

#include "test_fresnel.h"

int test_fresnel()
{
    float80 r_TE, r_TM, T;
    double omegap, gamma_;
    unittest_t test;
    casimir_t casimir;

    unittest_init(&test, "Fresnel", "Test Fresnel coefficients");

    T = 1;
    omegap = 1.32e2;
    gamma_ = 6.9e-1;
    casimir_init(&casimir, 1/0.5-1, T);

    casimir_set_drude(&casimir, omegap, gamma_, omegap, gamma_);

    casimir_rp(&casimir, 1*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.97252954726278);
    AssertAlmostEqual(&test, r_TM, +0.98616846109802);

    casimir_rp(&casimir, 10*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.85446954163643);
    AssertAlmostEqual(&test, r_TM, +0.85579839473205);

    casimir_rp(&casimir, 100*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.24595396364878);
    AssertAlmostEqual(&test, r_TM, +0.24598373253191);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}

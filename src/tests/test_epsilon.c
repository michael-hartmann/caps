#include "libcasimir.h"
#include "unittest.h"

#include "test_epsilon.h"

int test_epsilon()
{
    double omegap, gamma_;
    unittest_t test;

    unittest_init(&test, "Epsilon", "Test dielectric function");

    /* parameters of gold */
    omegap = 1.32e16;
    gamma_ = 6.9e13;

    AssertAlmostEqual(&test, casimir_epsilon(1e-10, omegap, gamma_), 2.5252173913043478e28);
    AssertAlmostEqual(&test, casimir_epsilon(1e-9,  omegap, gamma_), 2.5252173913043477e27);
    AssertAlmostEqual(&test, casimir_epsilon(1e-8,  omegap, gamma_), 2.5252173913043479e26);
    AssertAlmostEqual(&test, casimir_epsilon(1e-7,  omegap, gamma_), 2.5252173913043479e25);
    AssertAlmostEqual(&test, casimir_epsilon(1e-6,  omegap, gamma_), 2.5252173913043478e24);
    AssertAlmostEqual(&test, casimir_epsilon(1e-5,  omegap, gamma_), 2.5252173913043478e23);
    AssertAlmostEqual(&test, casimir_epsilon(1e-4,  omegap, gamma_), 2.525217391304348e22);
    AssertAlmostEqual(&test, casimir_epsilon(1e-3,  omegap, gamma_), 2.525217391304348e21);
    AssertAlmostEqual(&test, casimir_epsilon(1e-2,  omegap, gamma_), 2.5252173913043475e20);
    AssertAlmostEqual(&test, casimir_epsilon(1e-1,  omegap, gamma_), 2.525217391304344e19);
    AssertAlmostEqual(&test, casimir_epsilon(1e+0,  omegap, gamma_), 2.5252173913043113e18);
    AssertAlmostEqual(&test, casimir_epsilon(1e+1,  omegap, gamma_), 2.5252173913039818e17);
    AssertAlmostEqual(&test, casimir_epsilon(1e+2,  omegap, gamma_), 2.525217391300688e16);
    AssertAlmostEqual(&test, casimir_epsilon(1e+3,  omegap, gamma_), 2.5252173912677515e15);
    AssertAlmostEqual(&test, casimir_epsilon(1e+4,  omegap, gamma_), 2.5252173909383844e14);
    AssertAlmostEqual(&test, casimir_epsilon(1e+5,  omegap, gamma_), 2.5252173876447125e13);
    AssertAlmostEqual(&test, casimir_epsilon(1e+6,  omegap, gamma_), 2.5252173547079951e12);
    AssertAlmostEqual(&test, casimir_epsilon(1e+7,  omegap, gamma_), 2.5252170253408661e11);
    AssertAlmostEqual(&test, casimir_epsilon(1e+8,  omegap, gamma_), 2.5252137316743023e10);
    AssertAlmostEqual(&test, casimir_epsilon(1e+9,  omegap, gamma_), 2.5251807954812393e9);
    AssertAlmostEqual(&test, casimir_epsilon(1e+10, omegap, gamma_), 2.5248514808013332e8);
    AssertAlmostEqual(&test, casimir_epsilon(1e+11, omegap, gamma_), 2.5215630522431258e7);
    AssertAlmostEqual(&test, casimir_epsilon(1e+12, omegap, gamma_), 2489143.857142857);
    AssertAlmostEqual(&test, casimir_epsilon(1e+13, omegap, gamma_), 220557.9620253165);
    AssertAlmostEqual(&test, casimir_epsilon(1e+14, omegap, gamma_), 10311.05917159763);
    AssertAlmostEqual(&test, casimir_epsilon(1e+15, omegap, gamma_), 163.9934518241347);
    AssertAlmostEqual(&test, casimir_epsilon(1e+16, omegap, gamma_), 2.730459827192373);
    AssertAlmostEqual(&test, casimir_epsilon(1e+17, omegap, gamma_), 1.017411985729847);
    AssertAlmostEqual(&test, casimir_epsilon(1e+18, omegap, gamma_), 1.00017422797827);

    AssertAlmostEqual(&test, casimir_lnepsilon(1e-10, omegap, gamma_), 65.39870975842067);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e-8,  omegap, gamma_), 60.79353957243258);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e-6,  omegap, gamma_), 56.18836938644449);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e-2,  omegap, gamma_), 46.97802901446831);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+0,  omegap, gamma_), 42.3728588284802);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+2,  omegap, gamma_), 37.76768864249068);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+4,  omegap, gamma_), 33.16251845635911);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+6,  omegap, gamma_), 28.55734825602358);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+8,  omegap, gamma_), 23.95217663529314);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+10, omegap, gamma_), 19.34686298546513);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+12, omegap, gamma_), 14.7274493768442);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+14, omegap, gamma_), 9.240972304188087);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+16, omegap, gamma_), 1.00447002988231);
    AssertAlmostEqual(&test, casimir_lnepsilon(1e+18, omegap, gamma_), 1.7421280233797525e-4);

    return test_results(&test, stderr);
}

#include "sfunc.h"
#include "libcasimir.h"
#include "unittest.h"

#include "test_mie.h"

static double _mie_lna_perf(int l, double arg, sign_t *sign)
{
    casimir_t self;
    casimir_init(&self, 0.5, 2*arg);
    return casimir_lna_perf(&self, l, 1, sign);
}

static double _mie_lnb_perf(int l, double arg, sign_t *sign)
{
    double result;
    casimir_t self;

    casimir_init(&self, 0.5, 2*arg);
    result = casimir_lnb_perf(&self, l, 1, sign);
    casimir_free(&self);

    return result;
}

int test_mie(void)
{
    sign_t sign;
    casimir_t self;
    unittest_t test;
    unittest_init(&test, "Mie", "Test Mie functions al,bl for various parameters");

    casimir_init(&self, 0.5,2);

    /* b_l */
    AssertAlmostEqual(&test, _mie_lnb_perf(5,3,&sign), -3.206110089012862);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, _mie_lnb_perf(6,3,&sign), -6.093433624873396);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lnb_perf(50,1,&sign),   -365.8137152817732);
    AssertAlmostEqual(&test, _mie_lnb_perf(50,100,&sign),  174.3104974165916);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,2,&sign),  -726.3166073149845);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,100,&sign), 104.9919945452843);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,200,&sign), 349.7964954441692);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,300,&sign), 565.9447715085943);
    AssertAlmostEqual(&test, _mie_lnb_perf(40,0.01,&sign), -648.6664276814638);
    AssertAlmostEqual(&test, _mie_lnb_perf(4,0.01,&sign),  -52.95166526324419);

    /* a_l */
    AssertAlmostEqual(&test, _mie_lna_perf(3,3,&sign), 1.692450306201961);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lna_perf(4,3,&sign), -0.50863950281017);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,1,&sign),     -557.4493819729695);
    AssertAlmostEqual(&test, _mie_lna_perf(70,2,&sign),     -459.6943641467319);
    AssertAlmostEqual(&test, _mie_lna_perf(70,70,&sign),     73.02602649528605);
    AssertAlmostEqual(&test, _mie_lna_perf(70,100,&sign),    151.4135590544529);
    AssertAlmostEqual(&test, _mie_lna_perf(7,0.2,&sign),    -50.34157726932342);
    AssertAlmostEqual(&test, _mie_lna_perf(20,0.1,&sign),   -206.3146872637107);
    AssertAlmostEqual(&test, _mie_lna_perf(20,0.01,&sign),  -300.7209163862779);
    AssertAlmostEqual(&test, _mie_lna_perf(30,0.01,&sign),  -471.3445070668955);
    AssertAlmostEqual(&test, _mie_lna_perf(30,0.001,&sign), -611.8021993589887);

    return test_results(&test, stderr);
}

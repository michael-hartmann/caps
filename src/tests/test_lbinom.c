#include "sfunc.h"
#include "unittest.h"

#include "test_lbinom.h"

int test_lbinom()
{
    unittest_t test;

    unittest_init(&test, "Binomial", "Test lbinom function");

    AssertAlmostEqual(&test, lbinom(3,1),   log(3));
    AssertAlmostEqual(&test, lbinom(4,1),   log(4));
    AssertAlmostEqual(&test, lbinom(40,1),  log(40));
    AssertAlmostEqual(&test, lbinom(400,1), log(400));

    AssertAlmostEqual(&test, lbinom(49,6), 16.45341121889613);
    AssertAlmostEqual(&test, lbinom(1000,500), 689.4672615678519);

    return test_results(&test, stderr);
}

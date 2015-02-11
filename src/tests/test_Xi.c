#include "unittest.h"
#include "libcasimir.h"

#include "test_Xi.h"

int test_Xi(void)
{
    int sign;
    unittest_t test;
    unittest_init(&test, "Xi", "Test Xi function for various parameters");

    AssertAlmostEqual(&test, casimir_lnXi(1,1,0,&sign), -3.060270794691562);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnXi(1,1,1,&sign), -3.753417975251507);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnXi(3,2,1,&sign), -1.817138914330164);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, casimir_lnXi(4,3,2,&sign), -0.1101206735572);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, casimir_lnXi(4,2,2,&sign), -2.394730234408415);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, casimir_lnXi(11,1,0,&sign), 2.60283664295575);
    AssertAlmostEqual(&test, casimir_lnXi(11,7,0,&sign), 19.22557931884024);
    AssertAlmostEqual(&test, casimir_lnXi(11,7,5,&sign), 16.28731202862324);
    AssertAlmostEqual(&test, casimir_lnXi(201,7,5,&sign), 623.3839523251071);

    AssertAlmostEqual(&test, casimir_lnXi(100,10,0,&sign), 269.8159771440838);
    AssertAlmostEqual(&test, casimir_lnXi(100,10,1,&sign), 269.7633468887551);
    AssertAlmostEqual(&test, casimir_lnXi(100,10,10,&sign), 263.2542489687728);

    AssertAlmostEqual(&test, casimir_lnXi(100,100,100,&sign), 587.0039751538028);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,50,&sign),  696.7380895450116);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,0,&sign),  722.7572112350813);
    AssertAlmostEqual(&test, casimir_lnXi(100,100,1,&sign),  722.747260904228);
    AssertAlmostEqual(&test, casimir_lnXi(17,14,10,&sign),    45.8135805997528);

    return test_results(&test, stderr);
}

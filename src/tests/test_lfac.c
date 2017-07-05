#include "sfunc.h"
#include "unittest.h"

#include "test_lfac.h"

int test_lfac(void)
{
    unittest_t test;
    unittest_init(&test, "lfac", "Test log factorial", 1e-15);

    for(int i = 0; i < 200000; i++)
        AssertAlmostEqual(&test, lfac(i), lgamma(1+i));

    return test_results(&test, stderr);
}

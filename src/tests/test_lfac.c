#include <math.h>

#include "logfac.h"
#include "unittest.h"

#include "test_lfac.h"

int test_lfac(void)
{
    unittest_t test;
    unittest_init(&test, "lfac", "logarithm of factorial", 1e-15);

    for(int i = 0; i < 100000001; i++)
        AssertAlmostEqual(&test, lfac(i), lgamma(1+i));

    return test_results(&test, stderr);
}

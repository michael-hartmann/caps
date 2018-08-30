#include <math.h>

#include "logfac.h"
#include "unittest.h"

#include "test_logi.h"

int test_logi(void)
{
    unittest_t test;
    unittest_init(&test, "logi", "Test log for integer arguments", 1e-15);

    for(int i = 1; i < 100000001; i++)
        AssertAlmostEqual(&test, logi(i), log(i));

    return test_results(&test, stderr);
}

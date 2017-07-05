#include "sfunc.h"
#include "unittest.h"

#include "test_logi.h"

int test_logi(void)
{
    unittest_t test;
    unittest_init(&test, "logi", "Test log for integer arguments");

    for(int i = 1; i < 200001; i++)
        AssertAlmostEqual(&test, logi(i), log(i));

    return test_results(&test, stderr);
}

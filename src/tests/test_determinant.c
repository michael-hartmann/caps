#include "matrix.h"
#include "unittest.h"

#include "test_givens.h"

int test_givens()
{
    matrix_edouble_t *M;
    unittest_t test;
    unittest_init(&test, "QR decomposition", "Test QR decomposition using givens rotation");

    {
        M = matrix_edouble_alloc(2);
        matrix_set(M, 0,0, 20e100);
        matrix_set(M, 0,1, 2);
        matrix_set(M, 1,0, 1);
        matrix_set(M, 1,1, 1e-100);

        matrix_edouble_balance(M);

        AssertAlmostEqual(&test, matrix_edouble_logdet(M), log(18));

        matrix_edouble_free(M);
    }

    return test_results(&test, stderr);
}

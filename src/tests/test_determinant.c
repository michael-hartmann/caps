#include <stdbool.h>
#include <math.h>

#include "floattypes.h"
#include "matrix.h"
#include "libcasimir.h"
#include "unittest.h"

#include "test_determinant.h"

int test_determinant()
{
    const char *methods[] = { "QR_FLOAT80", "LU_FLOAT80", NULL };
    casimir_t casimir;
    unittest_t test;
    unittest_init(&test, "logdet", "Test computation of determinant");

    casimir_init(&casimir, 1, 1);

    matrix_float80 *M      = matrix_float80_alloc(2);
    matrix_sign_t  *M_sign = matrix_sign_alloc(2);

    matrix_set(M_sign, 0,0, +1);
    matrix_set(M_sign, 0,1, +1);
    matrix_set(M_sign, 1,0, +1);
    matrix_set(M_sign, 1,1, +1);

    for(int i = 0; methods[i] != NULL; i++)
    {
        const char *method = methods[i];

        casimir_set_detalg(&casimir, method);

        matrix_set(M, 0,0, log80(2L));
        matrix_set(M, 0,1, log80(20e1000L));
        matrix_set(M, 1,0, log80(1e-1000L));
        matrix_set(M, 1,1, log80(1L));

        AssertAlmostEqual(&test, matrix_logdet1mM(&casimir, M, M_sign), log(14));
    }

    matrix_float80_free(M);
    matrix_sign_free(M_sign);

    return test_results(&test, stderr);
}

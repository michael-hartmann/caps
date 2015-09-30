#include <stdio.h>

#include "unittest.h"

#include "test_bessels.h"
#include "test_casimirF.h"
#include "test_doublefact.h"
#include "test_epsilon.h"
#include "test_fresnel.h"
#include "test_determinant.h"
#include "test_integration_drude.h"
#include "test_integration_perf.h"
#include "test_Lambda.h"
#include "test_lbinom.h"
#include "test_logdet.h"
#include "test_mie.h"
#include "test_mie_drude.h"
#include "test_Plm.h"
#include "test_Xi.h"

/* This makes perfectly sense even on a single core machine: This maybe won't
 * speedup things if you have less than CORES cores, but it will test, if
 * locking and multithreading is correct.
 */
#define CORES 10


int main(int argc, char *argv[])
{
    test_Lambda();
    test_lbinom();
    //test_mie_drude();

    test_fresnel();
    //test_integration_drude();
    test_Plm();
    test_doublefact();
    test_epsilon();
    test_Xi();
    test_integration_perf();
    test_mie();
    test_besselI();
    test_besselK();
    test_determinant();

    test_logdet();

    return 0;
}

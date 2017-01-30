#include <stdio.h>

#include "unittest.h"

#include "test_bessels.h"
#include "test_casimirF.h"
#include "test_doublefact.h"
#include "test_fresnel.h"
#include "test_integration_drude.h"
#include "test_integration_perf.h"
#include "test_Lambda.h"
#include "test_logdet.h"
#include "test_logdet_HT.h"
#include "test_mie.h"
#include "test_mie_drude.h"
#include "test_Plm.h"

int main(int argc, char *argv[])
{
    test_Lambda();

    test_fresnel();
    test_doublefact();

    test_besselI();
    test_besselK();

    test_mie();
    test_mie_drude();

    test_Plm();


    //test_integration_perf();
    //test_integration_drude();

    //test_logdet_HT();
    //test_logdet();

	return 0;
}

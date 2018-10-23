#include <stdio.h>

#include "unittest.h"

#include "test_lfac.h"
#include "test_logi.h"
#include "test_bessels.h"
#include "test_fresnel.h"
#include "test_Lambda.h"
#include "test_mie.h"
#include "test_mie_drude.h"
#include "test_Plm.h"

int main(int argc, char *argv[])
{
    test_lfac();
    test_logi();
    test_Lambda();

    test_fresnel();

    test_besselI0();
    test_besselI0e();

    test_besselK0();
    test_besselK0e();

    test_besselI1();
    test_besselI1e();

    test_besselK1();
    test_besselK1e();

    test_besselI();
    test_besselK();

    test_log_besselKn();

    test_mie();
    test_mie_drude();

    test_Plm();

	return 0;
}

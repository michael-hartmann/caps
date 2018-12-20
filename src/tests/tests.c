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
#include "test_logdetD.h"

int main(int argc, char *argv[])
{
    test_lfac();
    test_logi();
    test_casimir_lnLambda();

    test_casimir_fresnel();

    test_bessel_I0();
    test_bessel_logI0();

    test_bessel_K0();
    test_bessel_logK0();

    test_bessel_I1();
    test_bessel_logI1();

    test_bessel_K1();
    test_bessel_logK1();

    test_bessel_logIn_half();
    test_bessel_logKn_half();

    test_bessel_logIn();
    test_bessel_logKn();

    test_bessel_ratioI();

    test_casimir_lnab_perf();
    test_casimir_lnab();

    test_lnPlm();

    test_logdetD();
    test_logdetD0();

	return 0;
}

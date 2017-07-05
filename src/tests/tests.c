#include <stdio.h>

#include "unittest.h"

#include "test_lfac.h"
#include "test_logi.h"
#include "test_bessels.h"
#include "test_casimirF.h"
#include "test_fresnel.h"
#include "test_Lambda.h"
#include "test_logdet.h"
#include "test_logdet_HT.h"
#include "test_mie.h"
#include "test_mie_drude.h"
#include "test_Plm.h"

int main(int argc, char *argv[])
{
    test_lfac();
    test_logi();
    test_Lambda();

    test_fresnel();

    test_besselI();
    test_besselK();

    test_mie();
    test_mie_drude();

    test_Plm();

	return 0;
}

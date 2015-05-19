#include <math.h>

#include "unittest.h"
#include "sfunc.h"
#include "libcasimir.h"

#include "test_Plm.h"

int test_Plm()
{
    sign_t sign;
    plm_combination_t comb;
    unittest_t test;

    unittest_init(&test, "Plm", "Test associated Legendre polynomials");

    /* Plm */
    AssertAlmostEqual(&test, plm_Plm(2,0,3), 13);
    AssertAlmostEqual(&test, plm_Plm(2,1,3), -25.455844122715710878430397035774565414);
    AssertAlmostEqual(&test, plm_Plm(2,2,3), -24);

    AssertAlmostEqual(&test, plm_Plm(3,0,2), 17);
    AssertAlmostEqual(&test, plm_Plm(3,1,2), -49.363448015713002865532220732917362457);
    AssertAlmostEqual(&test, plm_Plm(3,2,2), -90);
    AssertAlmostEqual(&test, plm_Plm(3,3,2), 77.9422863405994782087350853677642565124);


    AssertAlmostEqual(&test, plm_Plm(7,0,2), 2199.125);
    AssertAlmostEqual(&test, plm_Plm(7,1,2), -15209.246394437784569295351790159409371);
    AssertAlmostEqual(&test, plm_Plm(7,2,2), -88026.75);
    AssertAlmostEqual(&test, plm_Plm(7,3,2), 414721.162832537248613903297356202638370);
    AssertAlmostEqual(&test, plm_Plm(7,4,2), 1528065.0);
    AssertAlmostEqual(&test, plm_Plm(7,5,2), -4132071.3392037110374969861832293381768);
    AssertAlmostEqual(&test, plm_Plm(7,6,2), -7297290.0);
    AssertAlmostEqual(&test, plm_Plm(7,7,2), 6319638.51878214629264244945670369368228);

    AssertAlmostEqual(&test, plm_Plm(8,0,2), 7691.1484375);
    AssertAlmostEqual(&test, plm_Plm(8,1,2), -60890.462646434827363619067066431631293);
    AssertAlmostEqual(&test, plm_Plm(8,2,2), -413142.1875);
    AssertAlmostEqual(&test, plm_Plm(8,3,2), 2354110.35991671119020797796713392578050);
    AssertAlmostEqual(&test, plm_Plm(8,4,2), 1.0957629375e7);
    AssertAlmostEqual(&test, plm_Plm(8,5,2), -40024377.285620259853402179892456726654);
    AssertAlmostEqual(&test, plm_Plm(8,6,2), -1.076350275e8);
    AssertAlmostEqual(&test, plm_Plm(8,7,2), 189589155.563464388779273483701110810468);
    AssertAlmostEqual(&test, plm_Plm(8,8,2), 1.64189025e8);

    AssertAlmostEqual(&test, plm_Plm(3,2,200), -119997000);

    AssertAlmostEqual(&test, plm_Plm(30,30,10), -2.512712634569675964e70);
    AssertAlmostEqual(&test, plm_Plm(30,0,10), 1.022858298005579512e38);

    AssertAlmostEqual(&test, plm_lnPlm(300,0,100, &sign),    1586.0630493580697);
    AssertAlmostEqual(&test, plm_lnPlm(300,200,100, &sign),  2637.2261846173691);

    AssertAlmostEqual(&test, plm_lnPlm(300,1,100, &sign),    1591.7668317492472);
    AssertAlmostEqual(&test, plm_lnPlm(300,201,100, &sign),  2641.8313213287660);
    AssertAlmostEqual(&test, plm_lnPlm(300,201,1000, &sign), 3332.6176009928417);
    AssertAlmostEqual(&test, plm_lnPlm(300,201,5000, &sign), 3815.4490789776808);

    /* dPlm */
    AssertAlmostEqual(&test, plm_dPlm(3,0,3), 66);
    AssertAlmostEqual(&test, plm_dPlm(3,1,3), -197.2827919510467);
    AssertAlmostEqual(&test, plm_dPlm(3,2,3), -390);
    AssertAlmostEqual(&test, plm_dPlm(3,3,3), 381.8376618407357);


    plm_PlmPlm(4, 3, 2, 2, &comb);
    AssertAlmostEqual(&test, comb.lnPl1mPl2m, log(54675));
    AssertEqual(&test, comb.sign_Pl1mPl2m, +1);

    AssertAlmostEqual(&test, comb.lnPl1mdPl2m, log(100237.5));
    AssertEqual(&test, comb.sign_Pl1mdPl2m, +1);

    AssertAlmostEqual(&test, comb.lndPl1mPl2m, log(129600.0));
    AssertEqual(&test, comb.sign_dPl1mPl2m, +1);

    AssertAlmostEqual(&test, comb.lndPl1mdPl2m, log(237600.0));
    AssertEqual(&test, comb.sign_dPl1mdPl2m, +1);


    plm_PlmPlm(4, 3, 1, 2, &comb);
    AssertAlmostEqual(&test, comb.lnPl1mPl2m, log(10687.5));
    AssertEqual(&test, comb.sign_Pl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lnPl1mdPl2m, log(18375.0));
    AssertEqual(&test, comb.sign_Pl1mdPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mPl2m, log(24438.75));
    AssertEqual(&test, comb.sign_dPl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mdPl2m, log(42017.5));
    AssertEqual(&test, comb.sign_dPl1mdPl2m, -1);


    plm_PlmPlm(7, 5, 3, 2, &comb);
    AssertAlmostEqual(&test, comb.lnPl1mPl2m, 22.09944134068912);
    AssertEqual(&test, comb.sign_Pl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lnPl1mdPl2m, 23.20753237331177);
    AssertEqual(&test, comb.sign_Pl1mdPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mPl2m, 23.51706034623926);
    AssertEqual(&test, comb.sign_dPl1mPl2m, -1);

    AssertAlmostEqual(&test, comb.lndPl1mdPl2m, 24.62515137886192);
    AssertEqual(&test, comb.sign_dPl1mdPl2m, -1);

    return test_results(&test, stderr);
}

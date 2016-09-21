#include <math.h>
#include "unittest.h"
#include "sfunc.h"
#include "libcasimir.h"

#include "test_Plm.h"

int test_Plm()
{
    sign_t sign;
    unittest_t test;

    unittest_init(&test, "Plm", "Test associated Legendre polynomials");

    AssertAlmostEqual(&test, plm_lnPlm(1,0,1.1, &sign), 0.09531017980432493);
    AssertEqual(&test, sign, 1);
    //AssertAlmostEqual(&test, plm_lndPlm(1,0,1.1, &sign), 0.0);
    //AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,0,2, &sign), 0.6931471805599453);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1,0,2, &sign), 0.0);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,0,5, &sign), 1.6094379124341003);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1,0,5, &sign), 0.0);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,0,10, &sign), 2.3025850929940455);
    AssertEqual(&test, sign, 1);
    //AssertAlmostEqual(&test, plm_lndPlm(1,0,10, &sign), 0.0);
    //AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,0,100, &sign), 4.605170185988091);
    AssertEqual(&test, sign, 1);
    //AssertAlmostEqual(&test, plm_lndPlm(1,0,100, &sign), 0.0);
    //AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,0,1000, &sign), 6.907755278982137);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1,0,1000, &sign), 0.0);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,0,1e+06, &sign), 13.815510557964274);
    AssertEqual(&test, sign, 1);
    //AssertAlmostEqual(&test, plm_lndPlm(1,0,1e+06, &sign), 0.0);
    //AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,1.1, &sign), -0.7803238741323337);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,1.1, &sign), 0.8756340539366586);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,2, &sign), 0.5493061443340548);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,2, &sign), 0.14384103622589045);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,5, &sign), 1.5890269151739727);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,5, &sign), 0.020410997260127562);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,10, &sign), 2.297559925067295);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,10, &sign), 0.00502516792675072);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,100, &sign), 4.605120183487924);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,100, &sign), 5.0002500166679165e-05);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,1000, &sign), 6.907754778981887);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,1000, &sign), 5.000002500001666e-07);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1,1,1e+06, &sign), 13.815510557963773);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1,1,1e+06, &sign), 5.000000000002499e-13);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,1.1, &sign), 1.1310847224188405);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,1.1, &sign), 3.323647917439114);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,2, &sign), 5.224401683597868);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,2, &sign), 6.26696332875612);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,5, &sign), 10.06581896445358);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,5, &sign), 10.083932349317221);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,10, &sign), 13.5654694258405);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,10, &sign), 12.876787274823545);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,100, &sign), 25.08943299974896);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,100, &sign), 22.093745172685395);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,1000, &sign), 36.602468468510885);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,1000, &sign), 31.304151546407503);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,0,1e+06, &sign), 71.14124597453196);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5,0,1e+06, &sign), 58.935173329002225);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,1.1, &sign), 2.54332404330678);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,1.1, &sign), 4.914416354260865);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,2, &sign), 6.816269473090174);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,2, &sign), 7.8666845415730196);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,5, &sign), 11.672959264491194);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,5, &sign), 11.692005490359081);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,10, &sign), 15.17434719989084);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,10, &sign), 14.48588992235948);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,100, &sign), 26.69886535617332);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,100, &sign), 23.70317975159409);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,1000, &sign), 38.211906325389386);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,1000, &sign), 32.91358942550825);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,1,1e+06, &sign), 72.75068388696599);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,1,1e+06, &sign), 60.54461124143629);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,1.1, &sign), 2.949565556832074);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,1.1, &sign), 6.214961397335166);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,2, &sign), 9.597715649164016);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,2, &sign), 10.801688453489952);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,5, &sign), 14.796319503363605);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,5, &sign), 14.837141497883861);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,10, &sign), 18.338984552830215);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,10, &sign), 17.65588770812377);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,100, &sign), 29.876785844933366);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,100, &sign), 26.881153576379706);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,1000, &sign), 41.38995882240317);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,1000, &sign), 36.091642455855634);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5,5,1e+06, &sign), 75.9287377173126);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5,5,1e+06, &sign), 63.72266507178344);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,1.1, &sign), 2.9853524385054198);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,1.1, &sign), 5.987173686157405);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,2, &sign), 11.47273368360779);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,2, &sign), 13.217796987090871);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,5, &sign), 21.193564845688325);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,5, &sign), 21.906036437175313);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,10, &sign), 28.19740034187603);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,10, &sign), 28.2021602850787);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,100, &sign), 51.24678451823574);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,100, &sign), 48.944246795957994);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,1000, &sign), 74.27286993040246);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,1000, &sign), 69.6677002180988);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,0,1e+06, &sign), 143.35042508864336);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(10,0,1e+06, &sign), 131.8374996236736);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,1.1, &sign), 5.206849812025071);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,1.1, &sign), 8.24148507283957);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,2, &sign), 13.767103131424927);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,2, &sign), 15.513965921474947);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,5, &sign), 23.495063352349288);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,5, &sign), 24.207754976774506);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,10, &sign), 30.499720210145995);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,10, &sign), 30.50453335996177);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,100, &sign), 53.54936697944592);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,100, &sign), 51.24682978354086);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,1000, &sign), 76.5754549970807);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,1000, &sign), 71.97028529004021);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,1,1e+06, &sign), 145.65301018163737);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,1,1e+06, &sign), 134.14008471666762);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,1.1, &sign), 11.570910101939365);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,1.1, &sign), 15.032510219688945);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,2, &sign), 21.586765038605602);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,2, &sign), 23.374797829993906);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,5, &sign), 31.483362977337464);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,5, &sign), 32.201303959458265);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,10, &sign), 38.50769304102551);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,10, &sign), 38.51378130364866);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,100, &sign), 61.56363955417693);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,100, &sign), 59.26111499103548);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,1000, &sign), 84.5897901028007);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,1000, &sign), 79.98462052207613);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,5,1e+06, &sign), 153.66734591893618);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,5,1e+06, &sign), 142.15442045396654);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,1.1, &sign), 12.496493340755178);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,1.1, &sign), 16.455036361818216);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,2, &sign), 25.792793525419064);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,2, &sign), 27.689913510304944);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,5, &sign), 36.19000123381824);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,5, &sign), 36.92397040889844);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,10, &sign), 43.27533133275146);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,10, &sign), 43.28538166860496);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,100, &sign), 66.35093391695776);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,100, &sign), 64.04844882896404);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,1000, &sign), 89.37727987189739);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,1000, &sign), 84.7721106859098);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(10,10,1e+06, &sign), 158.45483766171625);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(10,10,1e+06, &sign), 146.941912196747);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,1.1, &sign), 19.916497983164188);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,1.1, &sign), 24.594490557528676);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,2, &sign), 63.3546620129539);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,2, &sign), 66.71581375537707);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,5, &sign), 112.09588794463903);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,5, &sign), 114.41867570065685);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,10, &sign), 147.13153649500967);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,10, &sign), 148.7459486857179);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,100, &sign), 262.3837545036006);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,100, &sign), 261.69065682045215);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,1000, &sign), 377.514234199385);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,1000, &sign), 374.51850242078075);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,0,1e+06, &sign), 722.9020105222214);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(50,0,1e+06, &sign), 712.9985229696858);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,1.1, &sign), 23.81416668339634);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,1.1, &sign), 28.493163894891097);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,2, &sign), 67.26511989971114);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,2, &sign), 70.62633930039378);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,5, &sign), 116.00770261583082);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,5, &sign), 118.33049879463907);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,10, &sign), 151.0435086107852);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,10, &sign), 152.65792284241425);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,100, &sign), 266.2957770039401);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,100, &sign), 265.60267934099573);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,1000, &sign), 381.42625719976263);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,1000, &sign), 378.4305254213604);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,1,1e+06, &sign), 726.8140335276496);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,1,1e+06, &sign), 716.9105459751139);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,1.1, &sign), 38.91341611149594);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,1.1, &sign), 43.61591507391466);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,2, &sign), 82.66938336265711);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,2, &sign), 86.03222377157388);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,5, &sign), 131.44451317185312);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,5, &sign), 133.7675114544048);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,10, &sign), 166.4840975823886);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,10, &sign), 168.09856079357897);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,100, &sign), 281.7375751890318);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,100, &sign), 281.0444780109849);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,1000, &sign), 396.8680673857681);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,1000, &sign), 393.8723356122144);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,5,1e+06, &sign), 742.2558438348672);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,5,1e+06, &sign), 732.3523562823315);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,1.1, &sign), 56.6687112091108);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,1.1, &sign), 61.43804136051512);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,2, &sign), 101.35554051734911);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,2, &sign), 104.72341236898295);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,5, &sign), 150.23218621891817);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,5, &sign), 152.5558155413511);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,10, &sign), 185.28357503153336);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,10, &sign), 186.89819127245508);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,100, &sign), 300.54083123177764);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,100, &sign), 299.84773556903224);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,1000, &sign), 415.67136093134974);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,1000, &sign), 412.67562917294754);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,10,1e+06, &sign), 761.0591377592366);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,10,1e+06, &sign), 751.1556502067009);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,1.1, &sign), 141.5880558691765);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,1.1, &sign), 147.15603680267364);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,2, &sign), 208.06955679249592);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,2, &sign), 211.5761146898159);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,5, &sign), 260.0555953344918);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,5, &sign), 262.3990024220061);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,10, &sign), 295.48224582915793);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,10, &sign), 297.10173407744554);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,100, &sign), 410.8602587501894);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,100, &sign), 410.1672115746298);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,1000, &sign), 525.9919885248875);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,1000, &sign), 522.996257251334);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(50,50,1e+06, &sign), 871.3797774739818);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(50,50,1e+06, &sign), 861.4762899214467);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,1.1, &sign), 41.74776988801881);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,1.1, &sign), 47.12617621116364);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,2, &sign), 128.85703644825924);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,2, &sign), 132.91212249045202);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,5, &sign), 226.37212162047132);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,5, &sign), 229.38816125872464);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,10, &sign), 296.44734878125274);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,10, &sign), 298.7549337258732);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,100, &sign), 526.9530490939379);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,100, &sign), 526.9530988451629);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,1000, &sign), 757.2140210491165);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,1000, &sign), 754.9114364536101);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,0,1e+06, &sign), 1447.9895738216865);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,0,1e+06, &sign), 1438.7792334497108);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,1.1, &sign), 46.3458523370313);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,1.1, &sign), 51.72450311083945);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,2, &sign), 133.46142863478607);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,2, &sign), 137.51653146661278);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,5, &sign), 230.97718817389858);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,5, &sign), 233.99322990660644);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,10, &sign), 301.0524936509405);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,10, &sign), 303.360079103188);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,100, &sign), 531.5582190286508);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,100, &sign), 531.5582687849014);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,1000, &sign), 761.819191232592);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,1000, &sign), 759.5166066371359);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,1,1e+06, &sign), 1452.5947440076745);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,1,1e+06, &sign), 1443.3844036356988);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,1.1, &sign), 64.49512261401553);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,1.1, &sign), 69.87960395506906);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,2, &sign), 151.76190508973625);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,2, &sign), 155.8174107022392);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,5, &sign), 249.29384752416067);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,5, &sign), 252.30993952112544);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,10, &sign), 319.3710325604637);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,10, &sign), 321.67863019560497);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,100, &sign), 549.8773594968228);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,100, &sign), 549.8774093736886);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,1000, &sign), 780.1383376710668);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,1000, &sign), 777.8357530768168);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,5,1e+06, &sign), 1470.9138905064508);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,5,1e+06, &sign), 1461.703550134475);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,1.1, &sign), 86.62947342742189);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,1.1, &sign), 92.03174127335247);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,2, &sign), 174.36600447369366);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,2, &sign), 178.42276667019004);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,5, &sign), 271.94849484328086);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,5, &sign), 274.96474388323514);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,10, &sign), 342.03155311883495);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,10, &sign), 344.33918882359166);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,100, &sign), 572.5397599017762);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,100, &sign), 572.5398101555642);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,1000, &sign), 802.8007567332141);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,1000, &sign), 800.4981721427329);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,10,1e+06, &sign), 1493.5763097570402);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,10,1e+06, &sign), 1484.3659693850645);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,1.1, &sign), 241.06115080063847);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,1.1, &sign), 246.83729225018888);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,2, &sign), 342.19333262236023);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,2, &sign), 346.28871093455);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,5, &sign), 441.3749906552409);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,5, &sign), 444.39623900441694);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,10, &sign), 511.6456869778622);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,10, &sign), 513.9545393707148);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,100, &sign), 742.2140295117537);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,100, &sign), 742.2140918268997);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,1000, &sign), 972.4756233714953);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,1000, &sign), 970.1730389016172);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,50,1e+06, &sign), 1663.2511824254705);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(100,50,1e+06, &sign), 1654.040842053495);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,1.1, &sign), 352.14550616761403);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,1.1, &sign), 358.40663428167113);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,2, &sign), 485.1085080142529);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,2, &sign), 489.30821309213286);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,5, &sign), 589.0805850982447);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,5, &sign), 592.117139366319);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,10, &sign), 659.9338860875769);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,10, &sign), 662.2465215164244);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,100, &sign), 890.6899119296398);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,100, &sign), 890.6900119346402);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,1000, &sign), 1120.9533714790362);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,1000, &sign), 1118.6507873860426);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(100,100,1e+06, &sign), 1811.7289493772248);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(100,100,1e+06, &sign), 1802.5186090052496);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,1.1, &sign), 218.369931617909);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,1.1, &sign), 225.363459818782);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,2, &sign), 654.8363204761886);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,2, &sign), 660.5014675508598);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,5, &sign), 1142.5410497358569);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,5, &sign), 1147.1666102773108);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,10, &sign), 1492.9327632020668);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,10, &sign), 1496.849806332538);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,100, &sign), 2645.4662763264428);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,100, &sign), 2647.0757641913233);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,1000, &sign), 3796.771185904025);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,1000, &sign), 3796.078039222965);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,0,1e+06, &sign), 7250.6489502698905);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,0,1e+06, &sign), 7243.048047810349);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,1.1, &sign), 224.58313594464965);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,1.1, &sign), 231.57667371912143);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,2, &sign), 661.0507736951938);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,2, &sign), 666.7159214375093);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,5, &sign), 1148.7556371924848);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,5, &sign), 1153.3811978173608);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,10, &sign), 1499.1473662576052);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,10, &sign), 1503.064409408299);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,100, &sign), 2651.6808843748113);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,100, &sign), 2653.290372239892);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,1000, &sign), 3802.985794001947);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,1000, &sign), 3802.2926473208886);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,1,1e+06, &sign), 7256.863558368313);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,1,1e+06, &sign), 7249.262655908771);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,1.1, &sign), 249.38781929785907);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,1.1, &sign), 256.38158678352545);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,2, &sign), 685.8854287299744);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,2, &sign), 691.550592495484);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,5, &sign), 1173.5935139152787);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,5, &sign), 1178.2190765422765);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,10, &sign), 1523.9856173540122);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,10, &sign), 1527.902660990047);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,100, &sign), 2676.5192552991175);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,100, &sign), 2678.1287431690034);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,1000, &sign), 3827.8241661155325);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,1000, &sign), 3827.1310194345224);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,5,1e+06, &sign), 7281.70193049391);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,5,1e+06, &sign), 7274.101028034369);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,1.1, &sign), 280.2850888002344);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,1.1, &sign), 287.27957345146046);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,2, &sign), 716.876338307201);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,2, &sign), 722.5415521418774);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,5, &sign), 1204.5944910827982);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,5, &sign), 1209.220059966375);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,10, &sign), 1554.987764436067);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,10, &sign), 1558.9048095887892);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,100, &sign), 2707.521776843168);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,100, &sign), 2709.1312647280706);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,1000, &sign), 3858.8266913760813);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,1000, &sign), 3858.133544695221);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,10,1e+06, &sign), 7312.704455791996);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,10,1e+06, &sign), 7305.103553332455);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,1.1, &sign), 523.073734528164);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,1.1, &sign), 530.0906398386494);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,2, &sign), 962.6445786135785);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,2, &sign), 968.3113920154349);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,5, &sign), 1450.6847466448119);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,5, &sign), 1455.3105156975275);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,10, &sign), 1801.1154548591724);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,10, &sign), 1805.0325485434562);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,100, &sign), 2953.6614498983918);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,100, &sign), 2955.2709382638227);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,1000, &sign), 4104.966483359222);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,1000, &sign), 4104.273336683168);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,50,1e+06, &sign), 7558.844248976339);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(500,50,1e+06, &sign), 7551.243346516798);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,1.1, &sign), 815.4370579952206);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,1.1, &sign), 822.5181575097513);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,2, &sign), 1264.1198067814141);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,2, &sign), 1269.7915860453404);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,5, &sign), 1753.1644352752814);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,5, &sign), 1757.790829339767);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,10, &sign), 2103.7120973975634);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,10, &sign), 2107.629342712578);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,100, &sign), 3256.2955362655425);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,100, &sign), 3257.9050261326215);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,1000, &sign), 4407.600941375931);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,1000, &sign), 4406.907794714891);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(500,100,1e+06, &sign), 7861.4787107468);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(500,100,1e+06, &sign), 7853.877808287259);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,1.1, &sign), 439.8074345368616);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,1.1, &sign), 447.4948126495353);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,2, &sign), 1312.9688009771617);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,2, &sign), 1319.3271727168515);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,5, &sign), 2288.4104333445875);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,5, &sign), 2293.729151392768);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,10, &sign), 2989.1977370442632);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,10, &sign), 2993.807929878004);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,100, &sign), 5294.266010535109);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,100, &sign), 5296.568645605588);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,1000, &sign), 7596.87584208464);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,1000, &sign), 7596.875842584391);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,0,1e+06, &sign), 14504.631370941557);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,0,1e+06, &sign), 14497.723615662577);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,1.1, &sign), 446.71448877540297);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,1.1, &sign), 454.4018692752367);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,2, &sign), 1319.8764788611857);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,2, &sign), 1326.2348507676643);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,5, &sign), 2295.318178307942);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,5, &sign), 2300.6368963769673);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,10, &sign), 2996.105489803071);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,10, &sign), 3000.715682641865);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,100, &sign), 5301.173765789076);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,100, &sign), 5303.476400859606);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,1000, &sign), 7603.783597363372);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,1000, &sign), 7603.783597863123);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,1,1e+06, &sign), 14511.53912622054);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,1,1e+06, &sign), 14504.631370941559);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,1.1, &sign), 474.31867010444955);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,1.1, &sign), 482.00610789270064);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,2, &sign), 1347.4956274665913);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,2, &sign), 1353.854003375984);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,5, &sign), 2322.9389368154248);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,5, &sign), 2328.2576553847152);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,10, &sign), 3023.726435401402);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,10, &sign), 3028.3366283614696);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,100, &sign), 5328.794771271238);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,100, &sign), 5331.0974063429685);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,1000, &sign), 7631.404603439875);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,1000, &sign), 7631.404603939638);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,5,1e+06, &sign), 14539.160132303046);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,5,1e+06, &sign), 14532.252377024064);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,1.1, &sign), 508.76974368580466);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,1.1, &sign), 516.4573604579883);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,2, &sign), 1381.993471123552);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,2, &sign), 1388.3518595418434);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,5, &sign), 2357.4418113933493);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,5, &sign), 2362.7605315259657);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,10, &sign), 3058.2298946378514);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,10, &sign), 3062.840087976899);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,100, &sign), 5363.298417644634);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,100, &sign), 5365.601052720117);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,1000, &sign), 7665.908251670591);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,1000, &sign), 7665.908252170391);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,10,1e+06, &sign), 14573.663780552522);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,10,1e+06, &sign), 14566.75602527354);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,1.1, &sign), 782.1991216934603);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,1.1, &sign), 789.8924323223661);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,2, &sign), 1656.9173691658418);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,2, &sign), 1663.2761577036067);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,5, &sign), 2632.526680489959);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,5, &sign), 2637.845450646424);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,10, &sign), 3333.3334725075015);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,10, &sign), 3337.943677973757);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,100, &sign), 5638.407983877662);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,100, &sign), 5640.7106190732175);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,1000, &sign), 7941.017877337837);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,1000, &sign), 7941.017877838837);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,50,1e+06, &sign), 14848.773406820066);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,50,1e+06, &sign), 14841.865651541086);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,1.1, &sign), 1118.4853446625614);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,1.1, &sign), 1126.19604114139);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,2, &sign), 1997.8478295053453);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,2, &sign), 2004.2078663549348);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,5, &sign), 2973.9599454242243);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,5, &sign), 2979.278871872943);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,10, &sign), 3674.8251986109312);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,10, &sign), 3679.435441972815);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,100, &sign), 5979.918423380138);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,100, &sign), 5982.221058950918);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,1000, &sign), 8282.528502572215);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,1000, &sign), 8282.528503076968);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(1000,100,1e+06, &sign), 15190.284033930384);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(1000,100,1e+06, &sign), 15183.376278651402);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,1.1, &sign), 883.0290901488979);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,1.1, &sign), 891.4099661730719);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,2, &sign), 2629.5801771345295);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,2, &sign), 2636.6317347634345);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,5, &sign), 4580.4955905256975);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,5, &sign), 4586.507460913569);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,10, &sign), 5982.074071765259);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,10, &sign), 5987.377413039962);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,100, &sign), 10592.21186505218);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,100, &sign), 10595.20764731573);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,1000, &sign), 15197.431540536312);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,1000, &sign), 15198.124688216747);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,0,1e+06, &sign), 29012.94259837524);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,0,1e+06, &sign), 29006.72799027682);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,1.1, &sign), 890.6296422989395);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,1.1, &sign), 899.0105189191265);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,2, &sign), 2637.1810409077684);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,2, &sign), 2644.2325985783555);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,5, &sign), 4588.096487828742);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,5, &sign), 4594.108358221823);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,10, &sign), 5989.67497296503);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,10, &sign), 5994.978314240998);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,100, &sign), 10599.81276749922);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,100, &sign), 10602.808549762782);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,1000, &sign), 15205.032442995729);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,1000, &sign), 15205.725590676166);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,1,1e+06, &sign), 29020.543500834785);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,1,1e+06, &sign), 29014.32889273636);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,1.1, &sign), 921.019840981977);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,1.1, &sign), 929.4007319062644);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,2, &sign), 2667.578718520729);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,2, &sign), 2674.6302771916817);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,5, &sign), 4618.494970156808);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,5, &sign), 4624.506840674922);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,10, &sign), 6020.073548814538);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,10, &sign), 6025.376890120816);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,100, &sign), 10630.211373283117);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,100, &sign), 10633.20715554698);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,1000, &sign), 15235.431049076726);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,1000, &sign), 15236.124196757164);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,5,1e+06, &sign), 29050.94210691878);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,5,1e+06, &sign), 29044.727498820357);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,1.1, &sign), 958.9805485310268);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,1.1, &sign), 967.3614841529879);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,2, &sign), 2705.5627973939354);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,2, &sign), 2712.614359191016);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,5, &sign), 4656.481563761844);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,5, &sign), 4662.493434670686);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,10, &sign), 6058.060434674034);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,10, &sign), 6063.363776075033);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,100, &sign), 10668.19835268758);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,100, &sign), 10671.19413495238);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,1000, &sign), 15273.418029409615);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,1000, &sign), 15274.111177090062);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,10,1e+06, &sign), 29088.929087261047);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,10,1e+06, &sign), 29082.714479162627);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,1.1, &sign), 1261.5811042164298);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,1.1, &sign), 1269.9634680563722);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,2, &sign), 3008.9109702192545);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,2, &sign), 3015.962632042111);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,5, &sign), 4959.910205706487);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,5, &sign), 4965.92208911848);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,10, &sign), 6361.498428723967);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,10, &sign), 6366.8017731560285);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,100, &sign), 10971.639340174106);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,100, &sign), 10974.635122468919);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,1000, &sign), 15576.85904660582);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,1000, &sign), 15577.552194286567);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,50,1e+06, &sign), 29392.370104757327);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,50,1e+06, &sign), 29386.155496658906);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,1.1, &sign), 1637.1036366485025);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,1.1, &sign), 1645.4904375321241);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,2, &sign), 3386.7665032145783);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,2, &sign), 3393.81847748902);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,5, &sign), 5338.017175967623);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,5, &sign), 5344.0290984499425);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,10, &sign), 6739.634623846346);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,10, &sign), 6744.937977750358);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,100, &sign), 11349.784889756313);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,100, &sign), 11352.780672144907);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,1000, &sign), 15955.004689030764);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,1000, &sign), 15955.69783671245);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(2000,100,1e+06, &sign), 29770.515748120008);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(2000,100,1e+06, &sign), 29764.301140021584);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,1.1, &sign), 2213.27569287925);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,1.1, &sign), 2222.573069871491);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,2, &sign), 6579.995754239426);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,2, &sign), 6587.9636258146675);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,5, &sign), 11457.33249056974);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,5, &sign), 11464.260654783699);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,10, &sign), 14961.284502089478);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,10, &sign), 14967.504134851995);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,100, &sign), 26486.63085401582);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,100, &sign), 26490.54292701875);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,1000, &sign), 37999.68006129632);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,1000, &sign), 38001.289499708706);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,0,1e+06, &sign), 72538.4577060812);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,0,1e+06, &sign), 72533.15938871466);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,1.1, &sign), 2221.792745997359);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,1.1, &sign), 2231.090123084887);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,2, &sign), 6588.512931959001);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,2, &sign), 6596.480803540911);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,5, &sign), 11465.849681698872);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,5, &sign), 11472.777845913664);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,10, &sign), 14969.801694777063);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,10, &sign), 14976.021327539782);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,100, &sign), 26495.148047202234);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,100, &sign), 26499.060120205166);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,1000, &sign), 38008.197254487684);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,1000, &sign), 38009.80669290007);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,1,1e+06, &sign), 72546.97489927262);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,1,1e+06, &sign), 72541.67658190608);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,1.1, &sign), 2255.8561564051015);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,1.1, &sign), 2265.1535357795283);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,2, &sign), 6622.579332800245);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,2, &sign), 6630.5472045421775);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,5, &sign), 11499.916404369476);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,5, &sign), 11506.84456860427);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,10, &sign), 15003.868454850488);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,10, &sign), 15010.088087618056);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,100, &sign), 26529.214819247612);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,100, &sign), 26533.12689225059);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,1000, &sign), 38042.26402665189);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,1000, &sign), 38043.873465064265);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,5,1e+06, &sign), 72581.04167143803);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,5,1e+06, &sign), 72575.74335407147);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,1.1, &sign), 2298.4246117857556);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,1.1, &sign), 2307.7219983066725);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,2, &sign), 6665.15713326439);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,2, &sign), 6673.125005506396);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,5, &sign), 11542.49521055019);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,5, &sign), 11549.42337484749);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,10, &sign), 15046.44737791502);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,10, &sign), 15052.66701069774);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,100, &sign), 26571.79377972449);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,100, &sign), 26575.705852727617);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,1000, &sign), 38084.84298750007);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,1000, &sign), 38086.45242591246);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,10,1e+06, &sign), 72623.62063228997);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,10,1e+06, &sign), 72618.32231492341);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,1.1, &sign), 2638.539374646917);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,1.1, &sign), 2647.836989801569);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,2, &sign), 7005.570921840936);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,2, &sign), 7013.538810085007);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,5, &sign), 11882.941181909891);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,5, &sign), 11889.869348207401);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,10, &sign), 15386.89708955445);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,10, &sign), 15393.11672282207);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,100, &sign), 26912.244688558854);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,100, &sign), 26916.15676156678);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,1000, &sign), 38425.293908216525);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,1000, &sign), 38426.90334662896);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,50,1e+06, &sign), 72964.07155312643);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,50,1e+06, &sign), 72958.77323575989);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,1.1, &sign), 3062.597933636185);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,1.1, &sign), 3071.896262597925);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,2, &sign), 7430.563724364048);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,2, &sign), 7438.531662611273);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,5, &sign), 12308.034553792588);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,5, &sign), 12314.962726340695);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,10, &sign), 15812.002149781383);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,10, &sign), 15818.221784564303);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,100, &sign), 27337.35349001806);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,100, &sign), 27341.265563040994);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,1000, &sign), 38850.40274680726);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,1000, &sign), 38852.01218521984);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, plm_lnPlm(5000,100,1e+06, &sign), 73389.1803920922);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, plm_lndPlm(5000,100,1e+06, &sign), 73383.88207472565);
    AssertEqual(&test, sign, 1);


    return test_results(&test, stderr);
}

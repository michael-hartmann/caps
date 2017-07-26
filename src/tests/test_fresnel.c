#include "libcasimir.h"
#include "misc.h"
#include "unittest.h"

#include "test_fresnel.h"

static double epsilonm1(double xi, void *userdata)
{
    double *param = (double *)userdata;
    double omegap_ = param[0];
    double gamma_ = param[1];

    return omegap_*omegap_/(xi*(xi+gamma_));
}

int test_fresnel()
{
    unittest_t test;
    casimir_t *casimir;
    double r_TE, r_TM;
    double userdata[2];

    unittest_init(&test, "Fresnel", "Test Fresnel coefficients", 1e-10);

    casimir = casimir_init(1); /* L/R doesn't matter */

    /* omegap = 500, gamma = 1 */
    userdata[0] = 500;
    userdata[1] = 1;
    casimir_set_epsilonm1(casimir, epsilonm1, userdata);

    casimir_rp(casimir, 1e-05, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9999821115267384);
    AssertAlmostEqual(&test, r_TM, 0.9999910557233689);

    casimir_rp(casimir, 1e-05, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9987358192020817);
    AssertAlmostEqual(&test, r_TM, 0.9999998735145679);

    casimir_rp(casimir, 1e-05, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9874305676854896);
    AssertAlmostEqual(&test, r_TM, 0.9999999873505795);

    casimir_rp(casimir, 1e-05, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8812556070593281);
    AssertAlmostEqual(&test, r_TM, 0.9999999987325552);

    casimir_rp(casimir, 1e-05, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3033354238790456);
    AssertAlmostEqual(&test, r_TM, 0.9999999998503327);

    casimir_rp(casimir, 1e-05, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.006173013760962271);
    AssertAlmostEqual(&test, r_TM, 0.9999999999190053);

    casimir_rp(casimir, 1e-05, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249156388294992e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999999999199891);

    casimir_rp(casimir, 1e-05, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249929688293448e-07);
    AssertAlmostEqual(&test, r_TM, 0.999999999919999);

    casimir_rp(casimir, 1e-05, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937499843761e-11);
    AssertAlmostEqual(&test, r_TM, 0.9999999999199991);

    casimir_rp(casimir, 0.001, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9998734473448324);
    AssertAlmostEqual(&test, r_TM, 0.9998734599980317);

    casimir_rp(casimir, 0.001, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9998210411567139);
    AssertAlmostEqual(&test, r_TM, 0.999910516574536);

    casimir_rp(casimir, 0.001, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9987289532131425);
    AssertAlmostEqual(&test, r_TM, 0.999987407449685);

    casimir_rp(casimir, 0.001, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9874237682392432);
    AssertAlmostEqual(&test, r_TM, 0.9999987344953795);

    casimir_rp(casimir, 0.001, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8812004995226734);
    AssertAlmostEqual(&test, r_TM, 0.9999998731926268);

    casimir_rp(casimir, 0.001, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3031749985875467);
    AssertAlmostEqual(&test, r_TM, 0.9999999850237502);

    casimir_rp(casimir, 0.001, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.006166983421224391);
    AssertAlmostEqual(&test, r_TM, 0.9999999918926167);

    casimir_rp(casimir, 0.001, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.242976675592933e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999999919910001);

    casimir_rp(casimir, 0.001, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243756165787261e-09);
    AssertAlmostEqual(&test, r_TM, 0.999999991992);

    casimir_rp(casimir, 0.01, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9995980855661183);
    AssertAlmostEqual(&test, r_TM, 0.9995980859679514);

    casimir_rp(casimir, 0.01, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9995960815997577);
    AssertAlmostEqual(&test, r_TM, 0.9996000799919992);

    casimir_rp(casimir, 0.01, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9994316547608904);
    AssertAlmostEqual(&test, r_TM, 0.9997157869861878);

    casimir_rp(casimir, 0.01, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9959681525576005);
    AssertAlmostEqual(&test, r_TM, 0.9999600007199952);

    casimir_rp(casimir, 0.01, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9605984476932515);
    AssertAlmostEqual(&test, r_TM, 0.9999959794469903);

    casimir_rp(casimir, 0.01, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707649424797745);
    AssertAlmostEqual(&test, r_TM, 0.9999995899653723);

    casimir_rp(casimir, 0.01, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05523407456316204);
    AssertAlmostEqual(&test, r_TM, 0.9999999097523485);

    casimir_rp(casimir, 0.01, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180472075911461);
    AssertAlmostEqual(&test, r_TM, 0.9999999191000682);

    casimir_rp(casimir, 0.01, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188118046024957e-08);
    AssertAlmostEqual(&test, r_TM, 0.9999999191999964);

    casimir_rp(casimir, 0.1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9986742297853705);
    AssertAlmostEqual(&test, r_TM, 0.9986742297986193);

    casimir_rp(casimir, 0.1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9986741635491135);
    AssertAlmostEqual(&test, r_TM, 0.9986742960315687);

    casimir_rp(casimir, 0.1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9986676218387484);
    AssertAlmostEqual(&test, r_TM, 0.9986808049946104);

    casimir_rp(casimir, 0.1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9981255928705577);
    AssertAlmostEqual(&test, r_TM, 0.9990623566417538);

    casimir_rp(casimir, 0.1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9867559171014392);
    AssertAlmostEqual(&test, r_TM, 0.9998679992371108);

    casimir_rp(casimir, 0.1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8758376689648629);
    AssertAlmostEqual(&test, r_TM, 0.9999867051026231);

    casimir_rp(casimir, 0.1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.2880199414204031);
    AssertAlmostEqual(&test, r_TM, 0.999998408022621);

    casimir_rp(casimir, 0.1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.005618154796490814);
    AssertAlmostEqual(&test, r_TM, 0.9999991100569692);

    casimir_rp(casimir, 0.1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681811725210102e-07);
    AssertAlmostEqual(&test, r_TM, 0.9999991199997743);

    casimir_rp(casimir, 1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9943591231228547);
    AssertAlmostEqual(&test, r_TM, 0.9943591231234171);

    casimir_rp(casimir, 1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9943591203106755);
    AssertAlmostEqual(&test, r_TM, 0.9943591259355948);

    casimir_rp(casimir, 1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9943588418841);
    AssertAlmostEqual(&test, r_TM, 0.9943594043481903);

    casimir_rp(casimir, 1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9943310688709504);
    AssertAlmostEqual(&test, r_TM, 0.994387038931944);

    casimir_rp(casimir, 1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9920319360002561);
    AssertAlmostEqual(&test, r_TM, 0.9960079840320636);

    casimir_rp(casimir, 1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9447423553480173);
    AssertAlmostEqual(&test, r_TM, 0.9994370576088932);

    casimir_rp(casimir, 1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5721068911608871);
    AssertAlmostEqual(&test, r_TM, 0.9999412171656573);

    casimir_rp(casimir, 1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02943722376917871);
    AssertAlmostEqual(&test, r_TM, 0.9999830297179542);

    casimir_rp(casimir, 1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124980468590091e-06);
    AssertAlmostEqual(&test, r_TM, 0.9999840001559985);

    casimir_rp(casimir, 10, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9589184175703337);
    AssertAlmostEqual(&test, r_TM, 0.9589184175703738);

    casimir_rp(casimir, 10, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9589184173692535);
    AssertAlmostEqual(&test, r_TM, 0.9589184177714538);

    casimir_rp(casimir, 10, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9589183974603404);
    AssertAlmostEqual(&test, r_TM, 0.9589184376803574);

    casimir_rp(casimir, 10, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9589164066208867);
    AssertAlmostEqual(&test, r_TM, 0.9589204284234534);

    casimir_rp(casimir, 10, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9587178388366593);
    AssertAlmostEqual(&test, r_TM, 0.9591180420875826);

    casimir_rp(casimir, 10, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9424043128251387);
    AssertAlmostEqual(&test, r_TM, 0.9707687779679319);

    casimir_rp(casimir, 10, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6579975863417028);
    AssertAlmostEqual(&test, r_TM, 0.9957447864986312);

    casimir_rp(casimir, 10, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05114929733954139);
    AssertAlmostEqual(&test, r_TM, 0.9990260277159961);

    casimir_rp(casimir, 10, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681753559802708e-06);
    AssertAlmostEqual(&test, r_TM, 0.9991207637323686);

    casimir_rp(casimir, 100, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707650746582664);
    AssertAlmostEqual(&test, r_TM, 0.670765074658269);

    casimir_rp(casimir, 100, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707650746450499);
    AssertAlmostEqual(&test, r_TM, 0.6707650746714855);

    casimir_rp(casimir, 100, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707650733364823);
    AssertAlmostEqual(&test, r_TM, 0.670765075980053);

    casimir_rp(casimir, 100, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707649424797745);
    AssertAlmostEqual(&test, r_TM, 0.6707652068367183);

    casimir_rp(casimir, 100, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707518572776928);
    AssertAlmostEqual(&test, r_TM, 0.670778291612797);

    casimir_rp(casimir, 100, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6694479999223603);
    AssertAlmostEqual(&test, r_TM, 0.6720779323589846);

    casimir_rp(casimir, 100, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5705716081269869);
    AssertAlmostEqual(&test, r_TM, 0.7512835664224575);

    casimir_rp(casimir, 100, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05474395800245473);
    AssertAlmostEqual(&test, r_TM, 0.9169408801095464);

    casimir_rp(casimir, 100, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188036039477654e-06);
    AssertAlmostEqual(&test, r_TM, 0.9252396718977455);

    casimir_rp(casimir, 1000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0556782899629886);
    AssertAlmostEqual(&test, r_TM, 0.0556782899629886);

    casimir_rp(casimir, 1000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0556782899629388);
    AssertAlmostEqual(&test, r_TM, 0.0556782899630384);

    casimir_rp(casimir, 1000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05567828995800809);
    AssertAlmostEqual(&test, r_TM, 0.05567828996796911);

    casimir_rp(casimir, 1000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05567828946493709);
    AssertAlmostEqual(&test, r_TM, 0.05567829046104011);

    casimir_rp(casimir, 1000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05567824015788102);
    AssertAlmostEqual(&test, r_TM, 0.05567833976809591);

    casimir_rp(casimir, 1000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05567330989462119);
    AssertAlmostEqual(&test, r_TM, 0.05568327002858566);

    casimir_rp(casimir, 1000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05518466758304567);
    AssertAlmostEqual(&test, r_TM, 0.05617188512657054);

    casimir_rp(casimir, 1000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02940952398039355);
    AssertAlmostEqual(&test, r_TM, 0.08187020111318079);

    casimir_rp(casimir, 1000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243053978387245e-06);
    AssertAlmostEqual(&test, r_TM, 0.1110062672721893);

    casimir_rp(casimir, 100000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249859378408116e-06);
    AssertAlmostEqual(&test, r_TM, 6.249859378408115e-06);

    casimir_rp(casimir, 100000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249859378408116e-06);
    AssertAlmostEqual(&test, r_TM, 6.249859378408116e-06);

    casimir_rp(casimir, 100000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249859378408054e-06);
    AssertAlmostEqual(&test, r_TM, 6.249859378408177e-06);

    casimir_rp(casimir, 100000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249859378401866e-06);
    AssertAlmostEqual(&test, r_TM, 6.249859378414365e-06);

    casimir_rp(casimir, 100000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249859377783137e-06);
    AssertAlmostEqual(&test, r_TM, 6.249859379033093e-06);

    casimir_rp(casimir, 100000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249859315910304e-06);
    AssertAlmostEqual(&test, r_TM, 6.249859440905927e-06);

    casimir_rp(casimir, 100000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249853128633108e-06);
    AssertAlmostEqual(&test, r_TM, 6.249865628183123e-06);

    casimir_rp(casimir, 100000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249234462773157e-06);
    AssertAlmostEqual(&test, r_TM, 6.250484294043074e-06);

    casimir_rp(casimir, 100000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124949219605698e-06);
    AssertAlmostEqual(&test, r_TM, 9.374769537088471e-06);

    /* omegap = 500, gamma = 1 */
    userdata[0] = 250;
    userdata[1] = 1;
    casimir_set_epsilonm1(casimir, epsilonm1, userdata);

    casimir_rp(casimir, 1e-05, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9999642233734758);
    AssertAlmostEqual(&test, r_TM, 0.9999821115267377);

    casimir_rp(casimir, 1e-05, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.997473237062022);
    AssertAlmostEqual(&test, r_TM, 0.9999997470290001);

    casimir_rp(casimir, 1e-05, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9750196192885491);
    AssertAlmostEqual(&test, r_TM, 0.9999999746996413);

    casimir_rp(casimir, 1e-05, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.7770009847510861);
    AssertAlmostEqual(&test, r_TM, 0.9999999974500067);

    casimir_rp(casimir, 1e-05, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1207834416131368);
    AssertAlmostEqual(&test, r_TM, 0.9999999995920751);

    casimir_rp(casimir, 1e-05, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.001557620650169676);
    AssertAlmostEqual(&test, r_TM, 0.9999999996789983);

    casimir_rp(casimir, 1e-05, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562435549914989e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999999996799868);

    casimir_rp(casimir, 1e-05, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562483886884955e-07);
    AssertAlmostEqual(&test, r_TM, 0.9999999996799966);

    casimir_rp(casimir, 1e-05, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562484375107422e-11);
    AssertAlmostEqual(&test, r_TM, 0.9999999996799968);

    casimir_rp(casimir, 0.001, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.999746910705746);
    AssertAlmostEqual(&test, r_TM, 0.9997469360089422);

    casimir_rp(casimir, 0.001, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9996421143411279);
    AssertAlmostEqual(&test, r_TM, 0.9998210411559972);

    casimir_rp(casimir, 0.001, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9974595224992556);
    AssertAlmostEqual(&test, r_TM, 0.9999748150429667);

    casimir_rp(casimir, 0.001, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9750061921718663);
    AssertAlmostEqual(&test, r_TM, 0.9999974688403741);

    casimir_rp(casimir, 0.001, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.7769043872341478);
    AssertAlmostEqual(&test, r_TM, 0.9999997448726469);

    casimir_rp(casimir, 0.001, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1206897124769097);
    AssertAlmostEqual(&test, r_TM, 0.9999999591748986);

    casimir_rp(casimir, 0.001, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.001556084932933882);
    AssertAlmostEqual(&test, r_TM, 0.9999999678681566);

    casimir_rp(casimir, 0.001, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.560890332224016e-05);
    AssertAlmostEqual(&test, r_TM, 0.999999967967001);

    casimir_rp(casimir, 0.001, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.560939056065999e-09);
    AssertAlmostEqual(&test, r_TM, 0.9999999679680008);

    casimir_rp(casimir, 0.01, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9991963326836761);
    AssertAlmostEqual(&test, r_TM, 0.9991963334870195);

    casimir_rp(casimir, 0.01, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9991923263660608);
    AssertAlmostEqual(&test, r_TM, 0.9992003199359936);

    casimir_rp(casimir, 0.01, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9988636325839745);
    AssertAlmostEqual(&test, r_TM, 0.9994316547379356);

    casimir_rp(casimir, 0.01, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9919525772610052);
    AssertAlmostEqual(&test, r_TM, 0.9999200025599899);

    casimir_rp(casimir, 0.01, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9227643524708381);
    AssertAlmostEqual(&test, r_TM, 0.9999919540411953);

    casimir_rp(casimir, 0.01, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4566790624587621);
    AssertAlmostEqual(&test, r_TM, 0.9999991334801211);

    casimir_rp(casimir, 0.01, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0150093828311961);
    AssertAlmostEqual(&test, r_TM, 0.9999996669502031);

    casimir_rp(casimir, 0.01, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0001546551227680541);
    AssertAlmostEqual(&test, r_TM, 0.9999996767001199);

    casimir_rp(casimir, 0.01, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.547029655104265e-08);
    AssertAlmostEqual(&test, r_TM, 0.9999996768000944);

    casimir_rp(casimir, 0.1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9973502178195816);
    AssertAlmostEqual(&test, r_TM, 0.9973502178460442);

    casimir_rp(casimir, 0.1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9973500855227875);
    AssertAlmostEqual(&test, r_TM, 0.9973503501362407);

    casimir_rp(casimir, 0.1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9973370194999864);
    AssertAlmostEqual(&test, r_TM, 0.9973633508380725);

    casimir_rp(casimir, 0.1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9962547007880417);
    AssertAlmostEqual(&test, r_TM, 0.9981255920465945);

    casimir_rp(casimir, 0.1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9736878167904921);
    AssertAlmostEqual(&test, r_TM, 0.9997359986522194);

    casimir_rp(casimir, 0.1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.7675352031965532);
    AssertAlmostEqual(&test, r_TM, 0.9999732362499939);

    casimir_rp(casimir, 0.1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1120071700908499);
    AssertAlmostEqual(&test, r_TM, 0.9999955920255312);

    casimir_rp(casimir, 0.1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.001416433422538172);
    AssertAlmostEqual(&test, r_TM, 0.9999964700266075);

    casimir_rp(casimir, 0.1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.420454141915045e-07);
    AssertAlmostEqual(&test, r_TM, 0.9999964800113903);

    casimir_rp(casimir, 1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9887501104825681);
    AssertAlmostEqual(&test, r_TM, 0.9887501104836867);

    casimir_rp(casimir, 1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.988750104890003);
    AssertAlmostEqual(&test, r_TM, 0.988750116076249);

    casimir_rp(casimir, 1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9887495511846907);
    AssertAlmostEqual(&test, r_TM, 0.9887506697539167);

    casimir_rp(casimir, 1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9886943199372402);
    AssertAlmostEqual(&test, r_TM, 0.9888056272671327);

    casimir_rp(casimir, 1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9841274880081918);
    AssertAlmostEqual(&test, r_TM, 0.9920318725140232);

    casimir_rp(casimir, 1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8925790435922286);
    AssertAlmostEqual(&test, r_TM, 0.9988730975679015);

    casimir_rp(casimir, 1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3401371081726173);
    AssertAlmostEqual(&test, r_TM, 0.9998700329142739);

    casimir_rp(casimir, 1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.007692755337015324);
    AssertAlmostEqual(&test, r_TM, 0.999935011885231);

    casimir_rp(casimir, 1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -7.812487792211345e-07);
    AssertAlmostEqual(&test, r_TM, 0.9999360039957474);

    casimir_rp(casimir, 10, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9195414884627832);
    AssertAlmostEqual(&test, r_TM, 0.9195414884628601);

    casimir_rp(casimir, 10, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9195414880773914);
    AssertAlmostEqual(&test, r_TM, 0.9195414888482518);

    casimir_rp(casimir, 10, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9195414499198051);
    AssertAlmostEqual(&test, r_TM, 0.9195415270058205);

    casimir_rp(casimir, 10, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9195376342647189);
    AssertAlmostEqual(&test, r_TM, 0.9195453424840438);

    casimir_rp(casimir, 10, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9191570990317911);
    AssertAlmostEqual(&test, r_TM, 0.9199241264695532);

    casimir_rp(casimir, 10, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8881721675136626);
    AssertAlmostEqual(&test, r_TM, 0.9423797926117482);

    casimir_rp(casimir, 10, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4404055424593053);
    AssertAlmostEqual(&test, r_TM, 0.9909960977887158);

    casimir_rp(casimir, 10, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01381344702557553);
    AssertAlmostEqual(&test, r_TM, 0.9963942611849895);

    casimir_rp(casimir, 10, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.42045049588218e-06);
    AssertAlmostEqual(&test, r_TM, 0.996492336991343);

    casimir_rp(casimir, 100, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4566792327935799);
    AssertAlmostEqual(&test, r_TM, 0.4566792327935832);

    casimir_rp(casimir, 100, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4566792327765481);
    AssertAlmostEqual(&test, r_TM, 0.456679232810615);

    casimir_rp(casimir, 100, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4566792310902326);
    AssertAlmostEqual(&test, r_TM, 0.4566792344969305);

    casimir_rp(casimir, 100, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4566790624587621);
    AssertAlmostEqual(&test, r_TM, 0.4566794031283676);

    casimir_rp(casimir, 100, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4566622001062995);
    AssertAlmostEqual(&test, r_TM, 0.4566962651460687);

    casimir_rp(casimir, 100, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4549838657866092);
    AssertAlmostEqual(&test, r_TM, 0.4583712892547404);

    casimir_rp(casimir, 100, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3384883210246955);
    AssertAlmostEqual(&test, r_TM, 0.5606841732543006);

    casimir_rp(casimir, 100, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01486512834749351);
    AssertAlmostEqual(&test, r_TM, 0.749296287526434);

    casimir_rp(casimir, 100, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.547023369368423e-06);
    AssertAlmostEqual(&test, r_TM, 0.7557429883092063);

    casimir_rp(casimir, 1000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01514030680688736);
    AssertAlmostEqual(&test, r_TM, 0.01514030680688736);

    casimir_rp(casimir, 1000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01514030680687267);
    AssertAlmostEqual(&test, r_TM, 0.01514030680690205);

    casimir_rp(casimir, 1000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01514030680541849);
    AssertAlmostEqual(&test, r_TM, 0.01514030680835623);

    casimir_rp(casimir, 1000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0151403066600005);
    AssertAlmostEqual(&test, r_TM, 0.01514030695377423);

    casimir_rp(casimir, 1000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01514029211821492);
    AssertAlmostEqual(&test, r_TM, 0.01514032149555979);

    casimir_rp(casimir, 1000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01513883808074223);
    AssertAlmostEqual(&test, r_TM, 0.01514177553296716);

    casimir_rp(casimir, 1000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.01499483161853911);
    AssertAlmostEqual(&test, r_TM, 0.01528578135426197);

    casimir_rp(casimir, 1000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.007685195055814629);
    AssertAlmostEqual(&test, r_TM, 0.02259373559522342);

    casimir_rp(casimir, 1000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.56077811057277e-06);
    AssertAlmostEqual(&test, r_TM, 0.03027211466533983);

    casimir_rp(casimir, 100000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562479492460476e-06);
    AssertAlmostEqual(&test, r_TM, 1.562479492460476e-06);

    casimir_rp(casimir, 100000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562479492460476e-06);
    AssertAlmostEqual(&test, r_TM, 1.562479492460476e-06);

    casimir_rp(casimir, 100000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.56247949246046e-06);
    AssertAlmostEqual(&test, r_TM, 1.562479492460492e-06);

    casimir_rp(casimir, 100000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562479492458914e-06);
    AssertAlmostEqual(&test, r_TM, 1.562479492462038e-06);

    casimir_rp(casimir, 100000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562479492304229e-06);
    AssertAlmostEqual(&test, r_TM, 1.562479492616723e-06);

    casimir_rp(casimir, 100000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.56247947683573e-06);
    AssertAlmostEqual(&test, r_TM, 1.562479508085222e-06);

    casimir_rp(casimir, 100000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562477929987429e-06);
    AssertAlmostEqual(&test, r_TM, 1.562481054933523e-06);

    casimir_rp(casimir, 100000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.562323260622633e-06);
    AssertAlmostEqual(&test, r_TM, 1.562635724298319e-06);

    casimir_rp(casimir, 100000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -7.812409669017971e-07);
    AssertAlmostEqual(&test, r_TM, 2.343718018017248e-06);

    /* omegap = 50, gamma = 1 */
    userdata[0] = 50;
    userdata[1] = 1;
    casimir_set_epsilonm1(casimir, epsilonm1, userdata);

    casimir_rp(casimir, 1e-05, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9998211296668196);
    AssertAlmostEqual(&test, r_TM, 0.9999105608335486);

    casimir_rp(casimir, 1e-05, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9874299494497019);
    AssertAlmostEqual(&test, r_TM, 0.9999987351213578);

    casimir_rp(casimir, 1e-05, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8812555519909525);
    AssertAlmostEqual(&test, r_TM, 0.9999998732555992);

    casimir_rp(casimir, 1e-05, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3033354222738578);
    AssertAlmostEqual(&test, r_TM, 0.9999999850332744);

    casimir_rp(casimir, 1e-05, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.006173013760358641);
    AssertAlmostEqual(&test, r_TM, 0.9999999919005373);

    casimir_rp(casimir, 1e-05, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249156388288805e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999999919989201);

    casimir_rp(casimir, 1e-05, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249929688293387e-07);
    AssertAlmostEqual(&test, r_TM, 0.99999999199991);

    casimir_rp(casimir, 1e-05, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937422501558e-09);
    AssertAlmostEqual(&test, r_TM, 0.99999999199992);

    casimir_rp(casimir, 1e-05, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937500617182e-13);
    AssertAlmostEqual(&test, r_TM, 0.99999999199992);

    casimir_rp(casimir, 0.001, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9987351939895296);
    AssertAlmostEqual(&test, r_TM, 0.9987353203774624);

    casimir_rp(casimir, 0.001, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9982118522976891);
    AssertAlmostEqual(&test, r_TM, 0.9991055259280784);

    casimir_rp(casimir, 0.001, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9873620702540485);
    AssertAlmostEqual(&test, r_TM, 0.9998740791616632);

    casimir_rp(casimir, 0.001, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8811949904706274);
    AssertAlmostEqual(&test, r_TM, 0.9999873199724475);

    casimir_rp(casimir, 0.001, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3031748380970555);
    AssertAlmostEqual(&test, r_TM, 0.999998502377256);

    casimir_rp(casimir, 0.001, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.006166983360919668);
    AssertAlmostEqual(&test, r_TM, 0.9999991892623231);

    casimir_rp(casimir, 0.001, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.242976674974956e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999991991006476);

    casimir_rp(casimir, 0.001, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243748446863764e-07);
    AssertAlmostEqual(&test, r_TM, 0.9999991991996412);

    casimir_rp(casimir, 0.001, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243756242976555e-11);
    AssertAlmostEqual(&test, r_TM, 0.9999991992006412);

    casimir_rp(casimir, 0.01, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9959881196293536);
    AssertAlmostEqual(&test, r_TM, 0.9959881236331661);

    casimir_rp(casimir, 0.01, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9959681525576005);
    AssertAlmostEqual(&test, r_TM, 0.9960079919992112);

    casimir_rp(casimir, 0.01, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9943310688709504);
    AssertAlmostEqual(&test, r_TM, 0.9971615001855431);

    casimir_rp(casimir, 0.01, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9604078384326366);
    AssertAlmostEqual(&test, r_TM, 0.9996000000326364);

    casimir_rp(casimir, 0.01, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6707518572776928);
    AssertAlmostEqual(&test, r_TM, 0.9999589994840637);

    casimir_rp(casimir, 0.01, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0552340256058798);
    AssertAlmostEqual(&test, r_TM, 0.9999909753121529);

    casimir_rp(casimir, 0.01, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180472014800374);
    AssertAlmostEqual(&test, r_TM, 0.9999919100715876);

    casimir_rp(casimir, 0.01, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188042226818314e-06);
    AssertAlmostEqual(&test, r_TM, 0.9999919199652877);

    casimir_rp(casimir, 0.01, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188118804222564e-10);
    AssertAlmostEqual(&test, r_TM, 0.9999919200652758);

    casimir_rp(casimir, 0.1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9868212089133505);
    AssertAlmostEqual(&test, r_TM, 0.9868212090442641);

    casimir_rp(casimir, 0.1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9868205544266518);
    AssertAlmostEqual(&test, r_TM, 0.9868218634986703);

    casimir_rp(casimir, 0.1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9867559171014392);
    AssertAlmostEqual(&test, r_TM, 0.986886181098994);

    casimir_rp(casimir, 0.1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9814135114656931);
    AssertAlmostEqual(&test, r_TM, 0.9906629630745358);

    casimir_rp(casimir, 0.1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8752654100535393);
    AssertAlmostEqual(&test, r_TM, 0.998677936772859);

    casimir_rp(casimir, 0.1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.2880041806789192);
    AssertAlmostEqual(&test, r_TM, 0.9998408272104711);

    casimir_rp(casimir, 0.1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.005618149296669908);
    AssertAlmostEqual(&test, r_TM, 0.9999110134941994);

    casimir_rp(casimir, 0.1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681172555554531e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999119977452066);

    casimir_rp(casimir, 0.1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681818117246386e-09);
    AssertAlmostEqual(&test, r_TM, 0.9999120077423187);

    casimir_rp(casimir, 1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9450088346090806);
    AssertAlmostEqual(&test, r_TM, 0.9450088346144241);

    casimir_rp(casimir, 1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9450088078935583);
    AssertAlmostEqual(&test, r_TM, 0.9450088613299337);

    casimir_rp(casimir, 1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9450061628622682);
    AssertAlmostEqual(&test, r_TM, 0.9450115062351053);

    casimir_rp(casimir, 1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9447423553480173);
    AssertAlmostEqual(&test, r_TM, 0.9452740649458414);

    casimir_rp(casimir, 1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9231360255795406);
    AssertAlmostEqual(&test, r_TM, 0.9607843260087915);

    casimir_rp(casimir, 1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5705716081269869);
    AssertAlmostEqual(&test, r_TM, 0.9941699558071257);

    casimir_rp(casimir, 1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02943447641403658);
    AssertAlmostEqual(&test, r_TM, 0.9983057461694517);

    casimir_rp(casimir, 1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.000312304527844981);
    AssertAlmostEqual(&test, r_TM, 0.9984015586184849);

    casimir_rp(casimir, 1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124999804375016e-08);
    AssertAlmostEqual(&test, r_TM, 0.9984025558107825);

    casimir_rp(casimir, 10, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6593462936120627);
    AssertAlmostEqual(&test, r_TM, 0.6593462936123332);

    casimir_rp(casimir, 10, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6593462922586002);
    AssertAlmostEqual(&test, r_TM, 0.6593462949657958);

    casimir_rp(casimir, 10, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6593461582524615);
    AssertAlmostEqual(&test, r_TM, 0.6593464289718917);

    casimir_rp(casimir, 10, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6593327581252062);
    AssertAlmostEqual(&test, r_TM, 0.6593598286717967);

    casimir_rp(casimir, 10, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6579975863417028);
    AssertAlmostEqual(&test, r_TM, 0.6606907706475974);

    casimir_rp(casimir, 10, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5571494526139611);
    AssertAlmostEqual(&test, r_TM, 0.7418686812987689);

    casimir_rp(casimir, 10, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05069630437282881);
    AssertAlmostEqual(&test, r_TM, 0.9108638336346566);

    casimir_rp(casimir, 10, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0005674803890383209);
    AssertAlmostEqual(&test, r_TM, 0.9190295152238632);

    casimir_rp(casimir, 10, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681817479338949e-08);
    AssertAlmostEqual(&test, r_TM, 0.9191176382393497);

    casimir_rp(casimir, 100, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05523407505767999);
    AssertAlmostEqual(&test, r_TM, 0.05523407505768097);

    casimir_rp(casimir, 100, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0552340750527353);
    AssertAlmostEqual(&test, r_TM, 0.05523407506262566);

    casimir_rp(casimir, 100, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05523407456316203);
    AssertAlmostEqual(&test, r_TM, 0.05523407555219892);

    casimir_rp(casimir, 100, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0552340256058798);
    AssertAlmostEqual(&test, r_TM, 0.05523412450948088);

    casimir_rp(casimir, 100, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.055229130317236);
    AssertAlmostEqual(&test, r_TM, 0.0552390197954157);

    casimir_rp(casimir, 100, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05474395800245473);
    AssertAlmostEqual(&test, r_TM, 0.05572416549707099);

    casimir_rp(casimir, 100, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02916230886533675);
    AssertAlmostEqual(&test, r_TM, 0.08123073914173493);

    casimir_rp(casimir, 100, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006119354129736071);
    AssertAlmostEqual(&test, r_TM, 0.1095276046553069);

    casimir_rp(casimir, 100, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188111857913927e-08);
    AssertAlmostEqual(&test, r_TM, 0.1101320974597512);

    casimir_rp(casimir, 1000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006235971494595939);
    AssertAlmostEqual(&test, r_TM, 0.0006235971494595939);

    casimir_rp(casimir, 1000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006235971494589712);
    AssertAlmostEqual(&test, r_TM, 0.0006235971494602167);

    casimir_rp(casimir, 1000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006235971493973119);
    AssertAlmostEqual(&test, r_TM, 0.0006235971495218758);

    casimir_rp(casimir, 1000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006235971432313952);
    AssertAlmostEqual(&test, r_TM, 0.0006235971556877927);

    casimir_rp(casimir, 1000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006235965266403287);
    AssertAlmostEqual(&test, r_TM, 0.0006235977722788592);

    casimir_rp(casimir, 1000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006235348736906763);
    AssertAlmostEqual(&test, r_TM, 0.0006236594252285067);

    casimir_rp(casimir, 1000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006174305398443036);
    AssertAlmostEqual(&test, r_TM, 0.0006297637590274571);

    casimir_rp(casimir, 1000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0003119930417258387);
    AssertAlmostEqual(&test, r_TM, 0.0009352011360943747);

    casimir_rp(casimir, 1000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243131151029384e-08);
    AssertAlmostEqual(&test, r_TM, 0.001247131382704279);

    casimir_rp(casimir, 100000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249936719390742e-08);
    AssertAlmostEqual(&test, r_TM, 6.249936719390741e-08);

    casimir_rp(casimir, 100000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249936719390741e-08);
    AssertAlmostEqual(&test, r_TM, 6.249936719390741e-08);

    casimir_rp(casimir, 100000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249936719390678e-08);
    AssertAlmostEqual(&test, r_TM, 6.249936719390803e-08);

    casimir_rp(casimir, 100000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249936719384491e-08);
    AssertAlmostEqual(&test, r_TM, 6.24993671939699e-08);

    casimir_rp(casimir, 100000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249936718765748e-08);
    AssertAlmostEqual(&test, r_TM, 6.249936720015733e-08);

    casimir_rp(casimir, 100000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249936656891382e-08);
    AssertAlmostEqual(&test, r_TM, 6.249936781890099e-08);

    casimir_rp(casimir, 100000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249930469461053e-08);
    AssertAlmostEqual(&test, r_TM, 6.249942969320428e-08);

    casimir_rp(casimir, 100000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249311788290028e-08);
    AssertAlmostEqual(&test, r_TM, 6.250561650491453e-08);

    casimir_rp(casimir, 100000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124968555003919e-08);
    AssertAlmostEqual(&test, r_TM, 9.37490488377755e-08);

    /* omegap = 5, gamma = 1 */
    userdata[0] = 5;
    userdata[1] = 1;
    casimir_set_epsilonm1(casimir, epsilonm1, userdata);

    casimir_rp(casimir, 1e-05, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9982127359741984);
    AssertAlmostEqual(&test, r_TM, 0.9991059681620672);

    casimir_rp(casimir, 1e-05, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8812500453089638);
    AssertAlmostEqual(&test, r_TM, 0.9999873262692934);

    casimir_rp(casimir, 1e-05, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3033352617551718);
    AssertAlmostEqual(&test, r_TM, 0.9999985033296771);

    casimir_rp(casimir, 1e-05, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.006173013699995678);
    AssertAlmostEqual(&test, r_TM, 0.9999991900543822);

    casimir_rp(casimir, 1e-05, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249156387670216e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999991998926463);

    casimir_rp(casimir, 1e-05, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249929688287199e-07);
    AssertAlmostEqual(&test, r_TM, 0.99999919999164);

    casimir_rp(casimir, 1e-05, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937422501496e-09);
    AssertAlmostEqual(&test, r_TM, 0.99999919999263);

    casimir_rp(casimir, 1e-05, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.24993749984376e-11);
    AssertAlmostEqual(&test, r_TM, 0.9999991999926399);

    casimir_rp(casimir, 1e-05, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937500624916e-15);
    AssertAlmostEqual(&test, r_TM, 0.99999919999264);

    casimir_rp(casimir, 0.001, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9874237682392432);
    AssertAlmostEqual(&test, r_TM, 0.9874250177801339);

    casimir_rp(casimir, 0.001, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9822619575420547);
    AssertAlmostEqual(&test, r_TM, 0.9910911187228267);

    casimir_rp(casimir, 0.001, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8806456340171707);
    AssertAlmostEqual(&test, r_TM, 0.9987390396369088);

    casimir_rp(casimir, 0.001, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.3031587899928083);
    AssertAlmostEqual(&test, r_TM, 0.999850260054616);

    casimir_rp(casimir, 0.001, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.006166977330453288);
    AssertAlmostEqual(&test, r_TM, 0.9999189327001241);

    casimir_rp(casimir, 0.001, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.242976613177203e-05);
    AssertAlmostEqual(&test, r_TM, 0.9999199164141183);

    casimir_rp(casimir, 0.001, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243748446245635e-07);
    AssertAlmostEqual(&test, r_TM, 0.999919926312305);

    casimir_rp(casimir, 0.001, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243756165781017e-09);
    AssertAlmostEqual(&test, r_TM, 0.999919926411293);

    casimir_rp(casimir, 0.001, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243756243748447e-13);
    AssertAlmostEqual(&test, r_TM, 0.9999199264122928);

    casimir_rp(casimir, 0.01, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9606003587320746);
    AssertAlmostEqual(&test, r_TM, 0.9606003973399145);

    casimir_rp(casimir, 0.01, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9604078384326366);
    AssertAlmostEqual(&test, r_TM, 0.9607920000316735);

    casimir_rp(casimir, 0.01, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9447423553480173);
    AssertAlmostEqual(&test, r_TM, 0.9719729999917939);

    casimir_rp(casimir, 0.01, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.6694479999223603);
    AssertAlmostEqual(&test, r_TM, 0.9959291586684536);

    casimir_rp(casimir, 0.01, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.055229130317236);
    AssertAlmostEqual(&test, r_TM, 0.9990983035211199);

    casimir_rp(casimir, 0.01, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180465903697762);
    AssertAlmostEqual(&test, r_TM, 0.9991916541664022);

    casimir_rp(casimir, 0.01, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188042165557455e-06);
    AssertAlmostEqual(&test, r_TM, 0.9991926423490807);

    casimir_rp(casimir, 0.01, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188118045406207e-08);
    AssertAlmostEqual(&test, r_TM, 0.9991926522370328);

    casimir_rp(casimir, 0.01, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188118811804541e-12);
    AssertAlmostEqual(&test, r_TM, 0.9991926523369017);

    casimir_rp(casimir, 0.1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8758434651694456);
    AssertAlmostEqual(&test, r_TM, 0.8758434663288353);

    casimir_rp(casimir, 0.1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8758376689648629);
    AssertAlmostEqual(&test, r_TM, 0.8758492622806944);

    casimir_rp(casimir, 0.1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8752654100535393);
    AssertAlmostEqual(&test, r_TM, 0.8764190191089246);

    casimir_rp(casimir, 0.1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.82915966461503);
    AssertAlmostEqual(&test, r_TM, 0.910395374100604);

    casimir_rp(casimir, 0.1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.2864375441243015);
    AssertAlmostEqual(&test, r_TM, 0.9843282962781196);

    casimir_rp(casimir, 0.1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.005617599368953802);
    AssertAlmostEqual(&test, r_TM, 0.9911786382328077);

    casimir_rp(casimir, 0.1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681166931838347e-05);
    AssertAlmostEqual(&test, r_TM, 0.991275777576993);

    casimir_rp(casimir, 0.1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681811668397731e-07);
    AssertAlmostEqual(&test, r_TM, 0.9912767546031145);

    casimir_rp(casimir, 0.1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681818181166839e-11);
    AssertAlmostEqual(&test, r_TM, 0.9912767644716538);

    casimir_rp(casimir, 1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5721224617164661);
    AssertAlmostEqual(&test, r_TM, 0.5721224617476084);

    casimir_rp(casimir, 1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5721223060200891);
    AssertAlmostEqual(&test, r_TM, 0.5721226174439442);

    casimir_rp(casimir, 1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5721068911608871);
    AssertAlmostEqual(&test, r_TM, 0.5721380318907949);

    casimir_rp(casimir, 1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5705716081269869);
    AssertAlmostEqual(&test, r_TM, 0.5736692348636366);

    casimir_rp(casimir, 1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4583736308584794);
    AssertAlmostEqual(&test, r_TM, 0.6674301434496239);

    casimir_rp(casimir, 1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02916230886533675);
    AssertAlmostEqual(&test, r_TM, 0.8543858497913782);

    casimir_rp(casimir, 1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0003122736320901164);
    AssertAlmostEqual(&test, r_TM, 0.8619887404668809);

    casimir_rp(casimir, 1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124977343944774e-06);
    AssertAlmostEqual(&test, r_TM, 0.8620681629049647);

    casimir_rp(casimir, 1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124999997734375e-10);
    AssertAlmostEqual(&test, r_TM, 0.8620689654369798);

    casimir_rp(casimir, 10, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0511539144368014);
    AssertAlmostEqual(&test, r_TM, 0.05115391443689374);

    casimir_rp(casimir, 10, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05115391397509605);
    AssertAlmostEqual(&test, r_TM, 0.05115391489859908);

    casimir_rp(casimir, 10, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05115386826173764);
    AssertAlmostEqual(&test, r_TM, 0.05115396061195727);

    casimir_rp(casimir, 10, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05114929733954139);
    AssertAlmostEqual(&test, r_TM, 0.05115853153196707);

    casimir_rp(casimir, 10, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05069630437282881);
    AssertAlmostEqual(&test, r_TM, 0.05161150302169435);

    casimir_rp(casimir, 10, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02690117400976085);
    AssertAlmostEqual(&test, r_TM, 0.07534646971399397);

    casimir_rp(casimir, 10, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0005619242053038956);
    AssertAlmostEqual(&test, r_TM, 0.1014847111738468);

    casimir_rp(casimir, 10, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.681185504524927e-06);
    AssertAlmostEqual(&test, r_TM, 0.1020351942921349);

    casimir_rp(casimir, 10, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -5.68181811854339e-10);
    AssertAlmostEqual(&test, r_TM, 0.1020408157642649);

    casimir_rp(casimir, 100, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180472076528684);
    AssertAlmostEqual(&test, r_TM, 0.0006180472076528806);

    casimir_rp(casimir, 100, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180472075911461);
    AssertAlmostEqual(&test, r_TM, 0.0006180472077146028);

    casimir_rp(casimir, 100, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180472014800374);
    AssertAlmostEqual(&test, r_TM, 0.0006180472138257115);

    casimir_rp(casimir, 100, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006180465903697762);
    AssertAlmostEqual(&test, r_TM, 0.0006180478249359728);

    casimir_rp(casimir, 100, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006179854854459903);
    AssertAlmostEqual(&test, r_TM, 0.000618108929859754);

    casimir_rp(casimir, 100, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0006119354129736071);
    AssertAlmostEqual(&test, r_TM, 0.0006241590022859687);

    casimir_rp(casimir, 100, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0003092146244940272);
    AssertAlmostEqual(&test, r_TM, 0.0009268796729160469);

    casimir_rp(casimir, 100, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.126775233353781e-06);
    AssertAlmostEqual(&test, r_TM, 0.001229967177221198);

    casimir_rp(casimir, 100, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.188112616110017e-10);
    AssertAlmostEqual(&test, r_TM, 0.001236093324329362);

    casimir_rp(casimir, 1000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243678275989208e-06);
    AssertAlmostEqual(&test, r_TM, 6.243678275989208e-06);

    casimir_rp(casimir, 1000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243678275982965e-06);
    AssertAlmostEqual(&test, r_TM, 6.243678275995451e-06);

    casimir_rp(casimir, 1000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243678275364849e-06);
    AssertAlmostEqual(&test, r_TM, 6.243678276613567e-06);

    casimir_rp(casimir, 1000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243678213553206e-06);
    AssertAlmostEqual(&test, r_TM, 6.24367833842521e-06);

    casimir_rp(casimir, 1000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243672032395143e-06);
    AssertAlmostEqual(&test, r_TM, 6.243684519583273e-06);

    casimir_rp(casimir, 1000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243053978387245e-06);
    AssertAlmostEqual(&test, r_TM, 6.244302573591171e-06);

    casimir_rp(casimir, 1000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.181860443500364e-06);
    AssertAlmostEqual(&test, r_TM, 6.305496108478004e-06);

    casimir_rp(casimir, 1000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.121858629784236e-06);
    AssertAlmostEqual(&test, r_TM, 9.365497922072481e-06);

    casimir_rp(casimir, 1000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.243131922767849e-10);
    AssertAlmostEqual(&test, r_TM, 1.248673223829943e-05);

    casimir_rp(casimir, 100000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.24993749281265e-10);
    AssertAlmostEqual(&test, r_TM, 6.249937492812649e-10);

    casimir_rp(casimir, 100000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.24993749281265e-10);
    AssertAlmostEqual(&test, r_TM, 6.24993749281265e-10);

    casimir_rp(casimir, 100000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937492812588e-10);
    AssertAlmostEqual(&test, r_TM, 6.249937492812712e-10);

    casimir_rp(casimir, 100000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937492806401e-10);
    AssertAlmostEqual(&test, r_TM, 6.2499374928189e-10);

    casimir_rp(casimir, 100000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937492187656e-10);
    AssertAlmostEqual(&test, r_TM, 6.249937493437643e-10);

    casimir_rp(casimir, 100000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249937430313276e-10);
    AssertAlmostEqual(&test, r_TM, 6.249937555312024e-10);

    casimir_rp(casimir, 100000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249931242881416e-10);
    AssertAlmostEqual(&test, r_TM, 6.249943742743885e-10);

    casimir_rp(casimir, 100000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -6.249312561557276e-10);
    AssertAlmostEqual(&test, r_TM, 6.250562424068024e-10);

    casimir_rp(casimir, 100000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -3.124968748359411e-10);
    AssertAlmostEqual(&test, r_TM, 9.374906237265888e-10);

    /* omegap = 1, gamma = 1 */
    userdata[0] = 1;
    userdata[1] = 1;
    casimir_set_epsilonm1(casimir, epsilonm1, userdata);

    casimir_rp(casimir, 1e-05, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9910955943251396);
    AssertAlmostEqual(&test, r_TM, 0.9955378194734733);

    casimir_rp(casimir, 1e-05, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5366572430987437);
    AssertAlmostEqual(&test, r_TM, 0.9999336730202796);

    casimir_rp(casimir, 1e-05, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02382278674124819);
    AssertAlmostEqual(&test, r_TM, 0.9999790240535241);

    casimir_rp(casimir, 1e-05, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0002498725780954969);
    AssertAlmostEqual(&test, r_TM, 0.9999799902027987);

    casimir_rp(casimir, 1e-05, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499962500328121e-06);
    AssertAlmostEqual(&test, r_TM, 0.9999800001000032);

    casimir_rp(casimir, 1e-05, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499974875250006e-08);
    AssertAlmostEqual(&test, r_TM, 0.999980000199);

    casimir_rp(casimir, 1e-05, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499974998999998e-10);
    AssertAlmostEqual(&test, r_TM, 0.9999800001999899);

    casimir_rp(casimir, 1e-05, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499975000237498e-12);
    AssertAlmostEqual(&test, r_TM, 0.9999800001999999);

    casimir_rp(casimir, 1e-05, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499975000249997e-16);
    AssertAlmostEqual(&test, r_TM, 0.9999800002);

    casimir_rp(casimir, 0.001, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.938690201292677);
    AssertAlmostEqual(&test, r_TM, 0.9386961378094758);

    casimir_rp(casimir, 0.001, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.9144270385869224);
    AssertAlmostEqual(&test, r_TM, 0.9562354025693161);

    casimir_rp(casimir, 0.001, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5349017311524594);
    AssertAlmostEqual(&test, r_TM, 0.9934213840780934);

    casimir_rp(casimir, 0.001, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02379807471143736);
    AssertAlmostEqual(&test, r_TM, 0.997904686648065);

    casimir_rp(casimir, 0.001, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0002496253277116806);
    AssertAlmostEqual(&test, r_TM, 0.9980010032413836);

    casimir_rp(casimir, 0.001, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497489997568163e-06);
    AssertAlmostEqual(&test, r_TM, 0.998001990025983);

    casimir_rp(casimir, 0.001, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497502372502381e-08);
    AssertAlmostEqual(&test, r_TM, 0.9980019998963074);

    casimir_rp(casimir, 0.001, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497502496252497e-10);
    AssertAlmostEqual(&test, r_TM, 0.9980019999950109);

    casimir_rp(casimir, 0.001, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497502497502373e-14);
    AssertAlmostEqual(&test, r_TM, 0.9980019999960078);

    casimir_rp(casimir, 0.01, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.8181899184325874);
    AssertAlmostEqual(&test, r_TM, 0.8181900820624044);

    casimir_rp(casimir, 0.01, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.817374313208272);
    AssertAlmostEqual(&test, r_TM, 0.8190024069066566);

    casimir_rp(casimir, 0.01, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.7532899862422072);
    AssertAlmostEqual(&test, r_TM, 0.8673093705872301);

    casimir_rp(casimir, 0.01, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1691695160099745);
    AssertAlmostEqual(&test, r_TM, 0.9722486742041565);

    casimir_rp(casimir, 0.01, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.002462824057661107);
    AssertAlmostEqual(&test, r_TM, 0.9801031646237987);

    casimir_rp(casimir, 0.01, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475122520328511e-05);
    AssertAlmostEqual(&test, r_TM, 0.9801989903215066);

    casimir_rp(casimir, 0.01, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475246274630722e-07);
    AssertAlmostEqual(&test, r_TM, 0.9801999510870444);

    casimir_rp(casimir, 0.01, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.47524751225125e-09);
    AssertAlmostEqual(&test, r_TM, 0.980199960694952);

    casimir_rp(casimir, 0.01, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475247524751226e-13);
    AssertAlmostEqual(&test, r_TM, 0.9801999607919918);

    casimir_rp(casimir, 0.1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.521143789973193);
    AssertAlmostEqual(&test, r_TM, 0.5211437932543151);

    casimir_rp(casimir, 0.1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.5211273867112726);
    AssertAlmostEqual(&test, r_TM, 0.521160196131157);

    casimir_rp(casimir, 0.1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.519510281845536);
    AssertAlmostEqual(&test, r_TM, 0.5227734921124086);

    casimir_rp(casimir, 0.1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.4038532922409105);
    AssertAlmostEqual(&test, r_TM, 0.6215781778486741);

    casimir_rp(casimir, 0.1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.02154315476997078);
    AssertAlmostEqual(&test, r_TM, 0.8124759394183629);

    casimir_rp(casimir, 0.1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.0002271467757869858);
    AssertAlmostEqual(&test, r_TM, 0.8195975818702496);

    casimir_rp(casimir, 0.1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.272714669503118e-06);
    AssertAlmostEqual(&test, r_TM, 0.8196713853830204);

    casimir_rp(casimir, 0.1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.272727146694223e-08);
    AssertAlmostEqual(&test, r_TM, 0.8196721236898685);

    casimir_rp(casimir, 0.1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.27272727271467e-12);
    AssertAlmostEqual(&test, r_TM, 0.8196721311467952);

    casimir_rp(casimir, 1, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1010205144253955);
    AssertAlmostEqual(&test, r_TM, 0.1010205144418921);

    casimir_rp(casimir, 1, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1010204319508072);
    AssertAlmostEqual(&test, r_TM, 0.101020596916479);

    casimir_rp(casimir, 1, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1010122668235371);
    AssertAlmostEqual(&test, r_TM, 0.1010287620298654);

    casimir_rp(casimir, 1, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.1002024332974939);
    AssertAlmostEqual(&test, r_TM, 0.1018384589811415);

    casimir_rp(casimir, 1, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.05572809000084122);
    AssertAlmostEqual(&test, r_TM, 0.1458980337503154);

    casimir_rp(casimir, 1, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.001234569782916727);
    AssertAlmostEqual(&test, r_TM, 0.198814520296916);

    casimir_rp(casimir, 1, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.249843769724048e-05);
    AssertAlmostEqual(&test, r_TM, 0.1999880014698181);

    casimir_rp(casimir, 1, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.249998437501973e-07);
    AssertAlmostEqual(&test, r_TM, 0.199999880000147);

    casimir_rp(casimir, 1, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.24999999984375e-11);
    AssertAlmostEqual(&test, r_TM, 0.199999999988);

    casimir_rp(casimir, 10, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.002262455019687151);
    AssertAlmostEqual(&test, r_TM, 0.002262455019691655);

    casimir_rp(casimir, 10, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.002262454997166996);
    AssertAlmostEqual(&test, r_TM, 0.00226245504221181);

    casimir_rp(casimir, 10, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.002262452767450922);
    AssertAlmostEqual(&test, r_TM, 0.002262457271927884);

    casimir_rp(casimir, 10, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.002262229818035647);
    AssertAlmostEqual(&test, r_TM, 0.002262680221342929);

    casimir_rp(casimir, 10, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.002240154610870701);
    AssertAlmostEqual(&test, r_TM, 0.002284755426257819);

    casimir_rp(casimir, 10, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.001133788305528256);
    AssertAlmostEqual(&test, r_TM, 0.003391115969619387);

    casimir_rp(casimir, 10, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.250123757945875e-05);
    AssertAlmostEqual(&test, r_TM, 0.004502386098661533);

    casimir_rp(casimir, 10, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.272498989874316e-07);
    AssertAlmostEqual(&test, r_TM, 0.004524659632581684);

    casimir_rp(casimir, 10, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.272727249896695e-11);
    AssertAlmostEqual(&test, r_TM, 0.004524886855101246);

    casimir_rp(casimir, 100, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475124995328485e-05);
    AssertAlmostEqual(&test, r_TM, 2.475124995328534e-05);

    casimir_rp(casimir, 100, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475124995081009e-05);
    AssertAlmostEqual(&test, r_TM, 2.475124995576009e-05);

    casimir_rp(casimir, 100, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475124970578485e-05);
    AssertAlmostEqual(&test, r_TM, 2.475125020078534e-05);

    casimir_rp(casimir, 100, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475122520328511e-05);
    AssertAlmostEqual(&test, r_TM, 2.475127470328508e-05);

    casimir_rp(casimir, 100, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.474877519827486e-05);
    AssertAlmostEqual(&test, r_TM, 2.475372470829533e-05);

    casimir_rp(casimir, 100, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.450620008333465e-05);
    AssertAlmostEqual(&test, r_TM, 2.499629982323256e-05);

    casimir_rp(casimir, 100, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.237593129072503e-05);
    AssertAlmostEqual(&test, r_TM, 3.712656860826392e-05);

    casimir_rp(casimir, 100, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.450738922292608e-07);
    AssertAlmostEqual(&test, r_TM, 4.925742598461207e-05);

    casimir_rp(casimir, 100, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.475245049384889e-11);
    AssertAlmostEqual(&test, r_TM, 4.950247512379332e-05);

    casimir_rp(casimir, 1000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497501249999531e-07);
    AssertAlmostEqual(&test, r_TM, 2.497501249999531e-07);

    casimir_rp(casimir, 1000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497501249997034e-07);
    AssertAlmostEqual(&test, r_TM, 2.497501250002029e-07);

    casimir_rp(casimir, 1000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497501249749782e-07);
    AssertAlmostEqual(&test, r_TM, 2.497501250249281e-07);

    casimir_rp(casimir, 1000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497501225024532e-07);
    AssertAlmostEqual(&test, r_TM, 2.497501274974531e-07);

    casimir_rp(casimir, 1000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497498752502027e-07);
    AssertAlmostEqual(&test, r_TM, 2.497503747497036e-07);

    casimir_rp(casimir, 1000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.497251524971772e-07);
    AssertAlmostEqual(&test, r_TM, 2.49775097502729e-07);

    casimir_rp(casimir, 1000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.472773527080236e-07);
    AssertAlmostEqual(&test, r_TM, 2.522228972918827e-07);

    casimir_rp(casimir, 1000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.24875093687541e-07);
    AssertAlmostEqual(&test, r_TM, 3.746251563123575e-07);

    casimir_rp(casimir, 1000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.49725277210055e-11);
    AssertAlmostEqual(&test, r_TM, 4.994752774721541e-07);

    casimir_rp(casimir, 100000, 1e-05, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499975000125e-11);
    AssertAlmostEqual(&test, r_TM, 2.499975000125e-11);

    casimir_rp(casimir, 100000, 0.001, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499975000125e-11);
    AssertAlmostEqual(&test, r_TM, 2.499975000125e-11);

    casimir_rp(casimir, 100000, 0.01, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499975000124975e-11);
    AssertAlmostEqual(&test, r_TM, 2.499975000125025e-11);

    casimir_rp(casimir, 100000, 0.1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.4999750001225e-11);
    AssertAlmostEqual(&test, r_TM, 2.4999750001275e-11);

    casimir_rp(casimir, 100000, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.499974999875003e-11);
    AssertAlmostEqual(&test, r_TM, 2.499975000374997e-11);

    casimir_rp(casimir, 100000, 10, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.49997497512525e-11);
    AssertAlmostEqual(&test, r_TM, 2.49997502512475e-11);

    casimir_rp(casimir, 100000, 100, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.4999725001525e-11);
    AssertAlmostEqual(&test, r_TM, 2.4999775000975e-11);

    casimir_rp(casimir, 100000, 1000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -2.49972502762225e-11);
    AssertAlmostEqual(&test, r_TM, 2.50022497262775e-11);

    casimir_rp(casimir, 100000, 100000, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -1.249987500093749e-11);
    AssertAlmostEqual(&test, r_TM, 3.74996250015625e-11);

    casimir_free(casimir);

    return test_results(&test, stderr);
}

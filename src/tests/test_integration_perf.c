#include <math.h>

#include "libcasimir.h"
#include "integration_perf.h"
#include "unittest.h"
#include "sfunc.h"

#include "test_integration_perf.h"

static void _integrals(int l1, int l2, int m, double nT, casimir_integrals_t *cint);
static void _integrals(int l1, int l2, int m, double nT, casimir_integrals_t *cint)
{
    integration_perf_t int_perf;
    int lmax = MAX(l1,l2);
    casimir_integrate_perf_init(&int_perf, nT, lmax);
    casimir_integrate_perf(&int_perf, l1, l2, m, cint);
    casimir_integrate_perf_free(&int_perf);
}

int test_integration_perf(void)
{
    double v;
    casimir_integrals_t cint;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    v = _I(2,2,2);
    AssertAlmostEqual(&test, v, 0.4054651081081643819780131163);

    v = _I(5,2,2);
    AssertAlmostEqual(&test, v, 5.813197260642067905077208633);

    v = _I(10,2,2);
    AssertAlmostEqual(&test, v, 15.95677753626913653914287233);

    v = _I(20,2,2);
    AssertAlmostEqual(&test, v, 42.90107196175527602778595691);

    v = _I(50,2,2);
    AssertAlmostEqual(&test, v, 148.6198399132092393065774679);

    v = _I(100,2,2);
    AssertAlmostEqual(&test, v, 363.5462729668148118110492452);

    v = _I(200,2,2);
    AssertAlmostEqual(&test, v, 862.6979732793756131251640376);

    v = _I(300,2,2);
    AssertAlmostEqual(&test, v, 1414.170985399411237979371551);

    v = _I(500,2,2);
    AssertAlmostEqual(&test, v, 2610.341684504699101672443343);

    v = _I(700,2,2);
    AssertAlmostEqual(&test, v, 3888.794466186317574388209390227783135);

    v = _I(1000,2,2);
    AssertAlmostEqual(&test, v, 5910.7939575865590314904753188706990739);

    v = _I(1500,2,2);
    AssertAlmostEqual(&test, v, 9472.8696066438596720790474385838167760);

    v = _I(3000,2,2);
    AssertAlmostEqual(&test, v, 21022.142076441619865987086582686107);

    // XXX
    _integrals(1050,1050,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 13940.125756903571190096123124106829015888756201576644);
    //AssertAlmostEqual(&test, cint.lnB_TM, 13967.951219063326784005854180747517807698817156658259); // ???
    //AssertAlmostEqual(&test, cint.lnC_TM, 13954.038134073832182440009136156538010794692114959778); // ???

    _integrals(1,1050,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 6244.4854968828077531090811790525063931088933444934980);

    _integrals(1,1,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -2.980829253011726);
    AssertAlmostEqual(&test, cint.lnB_TM, -2.0645385211375711716721);
    AssertAlmostEqual(&test, cint.lnC_TM, -2.5753641449035618548786);
    AssertAlmostEqual(&test, cint.lnD_TM, -2.5753641449035618548786);

    _integrals(241,73,1,30,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 406.63047665158294437064);
    AssertAlmostEqual(&test, cint.lnB_TM, 419.71230683599700819362);
    AssertAlmostEqual(&test, cint.lnC_TM, 412.57255550309976814896);
    AssertAlmostEqual(&test, cint.lnD_TM, 413.766977852356385781);

    _integrals(241,1,1,30,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 249.75276347175786475423);
    AssertAlmostEqual(&test, cint.lnB_TM, 258.05248402595679167552);
    AssertAlmostEqual(&test, cint.lnC_TM, 251.17334248392289626321);
    AssertAlmostEqual(&test, cint.lnD_TM, 256.62788585958419530558);

    _integrals(241,241,1,30,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 838.84852861683729524124);
    AssertAlmostEqual(&test, cint.lnB_TM, 853.98316452183914507246);
    AssertAlmostEqual(&test, cint.lnC_TM, 846.41479992430049881808);

    _integrals(3,2,1,2,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -4.094372316589062);
    AssertAlmostEqual(&test, cint.lnB_TM, -1.970116759119433);
    AssertAlmostEqual(&test, cint.lnC_TM, -3.298725852652321);

    _integrals(4,4,0,0.005,&cint);
    AssertAlmostEqual(&test, cint.lnB_TM, 56.977025325953406);

    _integrals(4,4,1,0.005,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +2.4806179125126554e17*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -2.2226323455151368e24*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.9457269656680333e20*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +6.9457269656680333e20*-2);

    _integrals(40,40,1,0.25,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +1.5754477603435539e159*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -6.3723632215476122e166*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -9.9568222699306801e162*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.9568222699306801e162*-2);

    _integrals(40,40,40,1,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +6.4140686579381969e91*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.0147301906459434e95*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -2.5352219594503741e93*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +2.5352219594503736e93*-2);

    _integrals(7,4,3,8.5,&cint);
    AssertAlmostEqual(&test, cint.signA_TM*exp(cint.lnA_TM), +4.8180365200137397e-9*-2);
    AssertAlmostEqual(&test, cint.signB_TM*exp(cint.lnB_TM), -1.3731640166794149e-8*-2);
    AssertAlmostEqual(&test, cint.signC_TM*exp(cint.lnC_TM), -6.7659079909128738e-9*-2);
    AssertAlmostEqual(&test, cint.signD_TM*exp(cint.lnD_TM), +9.44463292099617e-9*-2);

    _integrals(40,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 185.27722707813169721211989855051);

    _integrals(140,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 575.400220880156994701641252076629);

    _integrals(240,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 1025.59108523802829059595981668750);

    _integrals(540,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 2561.62734676892999652846813001817);

    _integrals(1540,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 8589.22040500894307686493465726593);

    _integrals(20, 20, 3, 10, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), 781.09859105555221959);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), -134864.11599496152485);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), -10117.073210379348893);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), 10117.073210379348893);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), -781.09859105555221959);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), 134864.11599496152485);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), 10117.073210379348893);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), -10117.073210379348893);

    _integrals(20, 20, 3, 0.005, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), 5.1118182099912547966e+132);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), -3.5441939544749868444e+141);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), -1.3290727331464213355e+137);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), 1.3290727331464213355e+137);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), -5.1118182099912547966e+132);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), 3.5441939544749868444e+141);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), 1.3290727331464213355e+137);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), -1.3290727331464213355e+137);

    _integrals(20, 10, 3, 0.1, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), 8.5978097063718598278e+56);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), -4.1556106748900681252e+62);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), -4.1556099454869528064e+59);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), 8.3112107466018794333e+59);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), -8.5978097063718598278e+56);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), 4.1556106748900681252e+62);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), 4.1556099454869528064e+59);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), -8.3112107466018794333e+59);

    _integrals(20, 15, 10, 5, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), -5431256655.4366550446);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), 204531605206.09912109);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), 28587590863.750656128);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), -37805365266.289375305);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), 5431256655.4366550446);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), -204531605206.09912109);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), -28587590863.750656128);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), 37805365266.289375305);

    _integrals(50, 15, 10, 10, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), -6.4110387204283793408e+19);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), 5.1666432341025142866e+21);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), 3.1581511099318560358e+20);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), -1.032906554654643454e+21);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), 6.4110387204283793408e+19);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), -5.1666432341025142866e+21);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), -3.1581511099318560358e+20);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), 1.032906554654643454e+21);

    _integrals(100, 25, 20, 20, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), -5.0294652559566390516e+36);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), 3.1964362960862542766e+38);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), 2.0270824787427330596e+37);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), -7.8698461944637818726e+37);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), 5.0294652559566390516e+36);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), -3.1964362960862542766e+38);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), -2.0270824787427330596e+37);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), 7.8698461944637818726e+37);

    _integrals(60, 55, 40, 0.11, &cint);
    AssertAlmostEqual(&test, cint.signA_TE * exp(cint.lnA_TE), -1.5670522864686087965e+280);
    AssertAlmostEqual(&test, cint.signB_TE * exp(cint.lnB_TE), 8.7546005438190589883e+285);
    AssertAlmostEqual(&test, cint.signC_TE * exp(cint.lnC_TE), 1.1165268702393727479e+283);
    AssertAlmostEqual(&test, cint.signD_TE * exp(cint.lnD_TE), -1.2180291188873130328e+283);
    AssertAlmostEqual(&test, cint.signA_TM * exp(cint.lnA_TM), 1.5670522864686087965e+280);
    AssertAlmostEqual(&test, cint.signB_TM * exp(cint.lnB_TM), -8.7546005438190589883e+285);
    AssertAlmostEqual(&test, cint.signC_TM * exp(cint.lnC_TM), -1.1165268702393727479e+283);
    AssertAlmostEqual(&test, cint.signD_TM * exp(cint.lnD_TM), 1.2180291188873130328e+283);

    return test_results(&test, stderr);
}

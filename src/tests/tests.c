#include <math.h>
#include <stdio.h>

#include "libcasimir.h"
#include "integration.h"
#include "integration_perf.h"
#include "sfunc.h"
#include "matrix.h"
#include "unittest.h"

#include "test_bessels.h"
#include "test_doublefact.h"
#include "test_epsilon.h"
#include "test_givens.h"
#include "test_integration_perf.h"
#include "test_Lambda.h"
#include "test_Plm.h"
#include "test_Xi.h"

#include "tests.h"

/* This makes perfectly sense even on a single core machine: This maybe won't
 * speedup things if you have less than CORES cores, but it will test, if
 * locking and multithreading is correct.
 */
#define CORES 10

int test_casimirF()
{
    unittest_t test;
    casimir_t casimir;
    double F;

    unittest_init(&test, "casimirF", "Compare free energies");

    casimir_init(&casimir, 0.85, 2.7);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, CORES);
    casimir_set_lmax(&casimir, 30);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -1.34361893570375);
    casimir_free(&casimir);

    casimir_init(&casimir, 0.7, 1);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, CORES);
    casimir_set_lmax(&casimir, 15);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -0.220709222562969);
    casimir_free(&casimir);

    casimir_init(&casimir, 0.85, 2.7);
    casimir_set_precision(&casimir, 1e-14);
    casimir_set_cores(&casimir, CORES);
    casimir_set_lmax(&casimir, 30);
    F = casimir_F(&casimir, NULL);
    AssertAlmostEqual(&test, F, -1.34361893570375);
    casimir_free(&casimir);

    return test_results(&test, stderr);
}

int test_logdet()
{
    unittest_t test;
    casimir_t casimir;
    integration_perf_t int_perf;
    const double RbyScriptL = 0.97;
    const double T = 0.1;
    const int lmax = 200;
    double logdet;

    unittest_init(&test, "logdet M", "calculate logdet");

    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_lmax(&casimir, lmax);
    casimir_set_cores(&casimir, CORES);


    logdet = casimir_logdetD(&casimir, 0, 0, NULL);
    AssertAlmostEqual(&test, logdet, -3.45236396285874);

    logdet = casimir_logdetD(&casimir, 0, 1, NULL);
    AssertAlmostEqual(&test, logdet, -2.63586999367158);

    logdet = casimir_logdetD(&casimir, 0, 10, NULL);
    AssertAlmostEqual(&test, logdet, -0.0276563864490425);

    casimir_integrate_perf_init(&int_perf, T, lmax);
    logdet = casimir_logdetD(&casimir, 1, 1, &int_perf);
    AssertAlmostEqual(&test, logdet, -2.63900987016801);
    casimir_integrate_perf_free(&int_perf);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}


static double _mie_lna_perf(int l, double arg, sign_t *sign)
{
    casimir_t self;
    casimir_init(&self, 0.5, 2*arg);
    return casimir_lna_perf(&self, l, 1, sign);
}

static double _mie_lnb_perf(int l, double arg, sign_t *sign)
{
    double result;
    casimir_t self;

    casimir_init(&self, 0.5, 2*arg);
    result = casimir_lnb_perf(&self, l, 1, sign);
    casimir_free(&self);

    return result;
}

int test_mie(void)
{
    sign_t sign;
    casimir_t self;
    unittest_t test;
    unittest_init(&test, "Mie", "Test Mie functions al,bl for various parameters");

    casimir_init(&self, 0.5,2);

    /* b_l */
    AssertAlmostEqual(&test, _mie_lnb_perf(5,3,&sign), -3.206110089012862);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, _mie_lnb_perf(6,3,&sign), -6.093433624873396);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lnb_perf(50,1,&sign),   -365.8137152817732);
    AssertAlmostEqual(&test, _mie_lnb_perf(50,100,&sign),  174.3104974165916);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,2,&sign),  -726.3166073149845);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,100,&sign), 104.9919945452843);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,200,&sign), 349.7964954441692);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,300,&sign), 565.9447715085943);
    AssertAlmostEqual(&test, _mie_lnb_perf(40,0.01,&sign), -648.6664276814638);
    AssertAlmostEqual(&test, _mie_lnb_perf(4,0.01,&sign),  -52.95166526324419);

    /* a_l */
    AssertAlmostEqual(&test, _mie_lna_perf(3,3,&sign), 1.692450306201961);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lna_perf(4,3,&sign), -0.50863950281017);
    AssertEqual(&test, sign, 1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,1,&sign),     -557.4493819729695);
    AssertAlmostEqual(&test, _mie_lna_perf(70,2,&sign),     -459.6943641467319);
    AssertAlmostEqual(&test, _mie_lna_perf(70,70,&sign),     73.02602649528605);
    AssertAlmostEqual(&test, _mie_lna_perf(70,100,&sign),    151.4135590544529);
    AssertAlmostEqual(&test, _mie_lna_perf(7,0.2,&sign),    -50.34157726932342);
    AssertAlmostEqual(&test, _mie_lna_perf(20,0.1,&sign),   -206.3146872637107);
    AssertAlmostEqual(&test, _mie_lna_perf(20,0.01,&sign),  -300.7209163862779);
    AssertAlmostEqual(&test, _mie_lna_perf(30,0.01,&sign),  -471.3445070668955);
    AssertAlmostEqual(&test, _mie_lna_perf(30,0.001,&sign), -611.8021993589887);

    return test_results(&test, stderr);
}

int test_mie_drude(void)
{
    double T, RbyScriptL, omegap, gamma_;
    sign_t sign_a, sign_b;
    double lna, lnb;
    casimir_t casimir;
    unittest_t test;
    unittest_init(&test, "Mie Drude", "Test Mie functions al,bl for various parameters");

    RbyScriptL = 0.85;
    T          = 2.7;
    omegap     = 1;
    gamma_     = 1;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir, gamma_);

    // void casimir_lnab(casimir_t *self, const int n, const int l, double *lna, double *lnb, int *sign_a, int *sign_b);

    casimir_lnab(&casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.4553021173541333720476206556765874126);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -5.8688410248499158590684691465517966950);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 1.77939374475276907218961292198964700248);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, 0.47334767009104316463610995593888380542);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 2, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -5.8888749734115470918000948427901346308);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -8.4332679619305632421065168698778165934);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    RbyScriptL = 0.95;
    T          = 0.1;
    omegap     = 0.1;
    gamma_     = 1.4;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir, gamma_);

    casimir_lnab(&casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -11.294081746210154580677515252479982979);
    AssertAlmostEqual(&test, lnb, -18.282872530972731935336084568403897476);

    casimir_lnab(&casimir, 1, 150, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -2128.7634085682857996403030143536620910);
    AssertAlmostEqual(&test, lnb, -2144.8884593589829670487681397791770641);

    casimir_lnab(&casimir, 100, 15, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -13.098866648021450586862321076066204163);
    AssertAlmostEqual(&test, lnb, -15.636661776753788796856285766794114116);

    casimir_lnab(&casimir, 200, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 5.85182666636940197335162042117627870470);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, 3.98653809954659615487558488229399701526);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    RbyScriptL = 0.5;
    T          = 1;
    omegap     = 1e-4;
    gamma_     = 1e-4;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir,  gamma_);

    casimir_lnab(&casimir, 1, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -55.208389631622379721466940745498367474);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -62.068501413438155604254192987355717706);
    AssertEqual(&test, sign_b, +1);

    casimir_free(&casimir);


    RbyScriptL = 0.5;
    T          = 1;
    omegap     = 1;
    gamma_     = 1e-4;
    casimir_init(&casimir, RbyScriptL, T);
    casimir_set_omegap_sphere(&casimir, omegap);
    casimir_set_gamma_sphere(&casimir,  gamma_);

    casimir_lnab(&casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.8540239333215465024615528927141405724);
    AssertAlmostEqual(&test, lnb, -7.2588880177134734386242016983196388100);

    casimir_lnab(&casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -7.6570049740166300703828407893184117471);
    AssertAlmostEqual(&test, lnb, -12.197169758793283148024686910238050626);

    casimir_lnab(&casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -12.482146722495979280724692478599547992);
    AssertAlmostEqual(&test, lnb, -17.727243821435497530589658157489862574);

    casimir_lnab(&casimir, 1, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -17.949100704846886740879599953685442860);
    AssertAlmostEqual(&test, lnb, -23.709990656446819196291385264231700122);

    casimir_lnab(&casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -23.891599427762872566501942669866619999);
    AssertAlmostEqual(&test, lnb, -30.060470310939164883544475705953387986);

    casimir_lnab(&casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -58.285652601186559317984617617473593988);
    AssertEqual(&test, sign_a, -1);
    AssertAlmostEqual(&test, lnb, -65.757440056322685520981556114527343549);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab(&casimir, 5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -0.7747752591093184484501684513821649046);
    AssertAlmostEqual(&test, lnb, -1.5892101743013427987322480836938505974);

    casimir_lnab(&casimir, 5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -1.8145275883758257249878124562886200521);
    AssertAlmostEqual(&test, lnb, -3.4950263771356228454682283405268413444);

    casimir_lnab(&casimir, 5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.6528777367390642776262971550796544139);
    AssertAlmostEqual(&test, lnb, -5.9236362448960108785813307258565097790);

    casimir_lnab(&casimir, 5, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -6.0449566310035356718674254329839576890);
    AssertAlmostEqual(&test, lnb, -8.7689157304390758571133982892467814487);

    casimir_lnab(&casimir, 5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -8.8673827346012321060116358013846191791);
    AssertAlmostEqual(&test, lnb, -11.960448909722605823940749038683258544);

    casimir_lnab(&casimir, 5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -27.402076826754931142375278917970072361);
    AssertAlmostEqual(&test, lnb, -31.720170045619969451321686623861774853);

    casimir_lnab(&casimir, 15, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 7.27860937372782962645683152471467237745);
    AssertAlmostEqual(&test, lnb, 7.19219783208384282132829359294763960883);

    casimir_lnab(&casimir, 15, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 7.02334016320151547642830308466323532804);
    AssertAlmostEqual(&test, lnb, 6.58397773056963416767295911583683756757);

    casimir_lnab(&casimir, 15, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 6.45115961754533443069229264656379407834);
    AssertAlmostEqual(&test, lnb, 5.69616122563189712145054916633333964902);

    casimir_lnab(&casimir, 15, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 5.58804392430302555818524298325189488817);
    AssertAlmostEqual(&test, lnb, 4.54964202124786251853718659079643076704);

    casimir_lnab(&casimir, 15, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 4.45913124438624191784643966629578158975);
    AssertAlmostEqual(&test, lnb, 3.16500182122169513370075408520449339632);

    casimir_lnab(&casimir, 15, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -4.4496934916915767898295519231813122808);
    AssertAlmostEqual(&test, lnb, -6.7270716450698327878555128172661915061);

    casimir_lnab(&casimir, 15, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -33.497220794437611205903569828087297055);
    AssertAlmostEqual(&test, lnb, -36.973404266432968593645620958704189729);

    casimir_lnab(&casimir, 15, 50, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -167.19982979934473385634246969571044075);
    AssertAlmostEqual(&test, lnb, -172.42020385124701110811733746273093103);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}


int test_integration_drude(void)
{
    casimir_integrals_t cint;
    casimir_t casimir;
    double omegap, gamma_;
    unittest_t test;

    unittest_init(&test, "Integration Drude", "Test integration for various parameters for Drude");

    omegap = 1.32e2;
    gamma_ = 6.9e-1;
    casimir_init(&casimir, 0.5, 1);
    casimir_set_omegap_plane(&casimir, omegap);
    casimir_set_gamma_plane(&casimir, gamma_);

    {
        casimir_integrate_drude(&casimir, &cint, 3, 2, 1, 1, 1);

        AssertAlmostEqual(&test, cint.lnA_TE, -0.62981145199252068602408);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -0.59589434712666196879313);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, 2.7383266248198112347328);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, 2.7945369735963442460566);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, 0.74158885587407677484842);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, 0.78654850992186874297884);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, 1.1256865694785886272321);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, 1.1711803659444285938208);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 3, 2, 2, 1, 1);

        AssertAlmostEqual(&test, cint.lnA_TE, -0.7132786835392315505014);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -0.67405454481613061169669);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, 1.5622997116874727152691);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, 1.6191154482067796685377);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, 0.18064782184885164701636);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, 0.22781988160375950931957);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, 0.52072721231985447793985);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, 0.56889321108023160164297);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 3, 2, 1, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, -4.1459191747317624709052);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -4.1197945729671869295841);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, -2.0359815567492122711267);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, -1.9903481981146134279244);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, -3.356496453224571521287);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, -3.3216912696961575288389);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, -3.0252864555803092122331);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, -2.9891216277980328174213);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 2, 1, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, -3.4410543260111500276103);
        AssertEqual(&test, cint.signA_TE, +1);
        AssertAlmostEqual(&test, cint.lnA_TM, -3.4078967115562010093212);
        AssertEqual(&test, cint.signA_TM, -1);

        AssertAlmostEqual(&test, cint.lnB_TE, -0.74990261955450098567219);
        AssertEqual(&test, cint.signB_TE, -1);
        AssertAlmostEqual(&test, cint.lnB_TM, -0.69545756637556400940693);
        AssertEqual(&test, cint.signB_TM, +1);

        AssertAlmostEqual(&test, cint.lnC_TE, -2.5076062999924391950641);
        AssertEqual(&test, cint.signC_TE, +1);
        AssertAlmostEqual(&test, cint.lnC_TM, -2.4646103122548521433085);
        AssertEqual(&test, cint.signC_TM, -1);

        AssertAlmostEqual(&test, cint.lnD_TE, -1.883875993266507818876);
        AssertEqual(&test, cint.signD_TE, -1);
        AssertAlmostEqual(&test, cint.lnD_TM, -1.8392724297362885494622);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 3, 1, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, -2.6626325469377015011493);
        AssertEqual(&test, cint.signA_TE, -1);
        AssertAlmostEqual(&test, cint.lnA_TM, -2.6217256954918375062335);
        AssertEqual(&test, cint.signA_TM, +1);

        AssertAlmostEqual(&test, cint.lnB_TE, 0.69662250934949833938181);
        AssertEqual(&test, cint.signB_TE, +1);
        AssertAlmostEqual(&test, cint.lnB_TM, 0.76016712598231202383933);
        AssertEqual(&test, cint.signB_TM, -1);

        AssertAlmostEqual(&test, cint.lnC_TE, -1.2220766351698472482591);
        AssertEqual(&test, cint.signC_TE, -1);
        AssertAlmostEqual(&test, cint.lnC_TM, -1.1696039758614060731918);
        AssertEqual(&test, cint.signC_TM, +1);

        AssertAlmostEqual(&test, cint.lnD_TE, -0.94454349314724581653695);
        AssertEqual(&test, cint.signD_TE, +1);
        AssertAlmostEqual(&test, cint.lnD_TM, -0.89174332617421184874047);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 3, 1, 1, 20);

        AssertAlmostEqual(&test, cint.lnA_TE, -43.137947441949791493356);
        AssertAlmostEqual(&test, cint.lnA_TM, -43.120372668298594245786);
        AssertAlmostEqual(&test, cint.lnB_TE, -42.712836955286222313897);
        AssertAlmostEqual(&test, cint.lnB_TM, -42.688853217438179634812);
        AssertAlmostEqual(&test, cint.lnC_TE, -42.981679865536465578416);
        AssertAlmostEqual(&test, cint.lnC_TM, -42.961747317757136470438);
        AssertAlmostEqual(&test, cint.lnD_TE, -42.896313862496345284576);
        AssertAlmostEqual(&test, cint.lnD_TM, -42.875331151064598797301);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 4, 3, 1, 1, 0.01);

        AssertAlmostEqual(&test, cint.lnA_TE, 29.361786303876121509058);
        AssertAlmostEqual(&test, cint.lnA_TM, 29.727644399301097042065);
        AssertAlmostEqual(&test, cint.lnB_TE, 43.288690290377248702015);
        AssertAlmostEqual(&test, cint.lnB_TM, 43.774264506257246488789);
        AssertAlmostEqual(&test, cint.lnC_TE, 36.104129446007204258962);
        AssertAlmostEqual(&test, cint.lnC_TM, 36.53003707484319590354);
        AssertAlmostEqual(&test, cint.lnD_TE, 36.39181148264103954446);
        AssertAlmostEqual(&test, cint.lnD_TM, 36.817719115542167950645);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 20, 20, 0, 1, 2);

        AssertAlmostEqual(&test, cint.lnB_TE, 80.261829579383622339087);
        AssertAlmostEqual(&test, cint.lnB_TM, 80.616659994373914035408);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 20, 20, 1, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, 69.659355648443184396324);
        AssertAlmostEqual(&test, cint.lnA_TM, 69.996776884623239716192);
        AssertAlmostEqual(&test, cint.lnB_TE, 80.212799262982493059119);
        AssertAlmostEqual(&test, cint.lnB_TM, 80.567633810049211725342);
        AssertAlmostEqual(&test, cint.lnC_TE, 74.92332975118821662198);
        AssertAlmostEqual(&test, cint.lnC_TM, 75.269464691472079467208);
        AssertAlmostEqual(&test, cint.lnD_TE, 74.92332975118821662198);
        AssertAlmostEqual(&test, cint.lnD_TM, 75.269464691472079467208);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 7, 7, 1, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, 6.6971469912051709882689);
        AssertAlmostEqual(&test, cint.lnA_TM, 6.8044855533996484279007);
        AssertAlmostEqual(&test, cint.lnB_TE, 12.986970305176281775914);
        AssertAlmostEqual(&test, cint.lnB_TM, 13.113563515622480422875);
        AssertAlmostEqual(&test, cint.lnC_TE, 9.8017901880975337522829);
        AssertAlmostEqual(&test, cint.lnC_TM, 9.9188754798338178698508);
        AssertAlmostEqual(&test, cint.lnD_TE, 9.8017901880975337522829);
        AssertAlmostEqual(&test, cint.lnD_TM, 9.9188754798338178698508);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 60, 7, 1, 1, 4);

        AssertAlmostEqual(&test, cint.lnA_TE, 111.142136572991446682);
        AssertAlmostEqual(&test, cint.lnA_TM, 111.67030961586370531194);
        AssertAlmostEqual(&test, cint.lnB_TE, 121.40133118469416742774);
        AssertAlmostEqual(&test, cint.lnB_TM, 121.94548929011007049376);
        AssertAlmostEqual(&test, cint.lnC_TE, 115.18974866771901373092);
        AssertAlmostEqual(&test, cint.lnC_TM, 115.72592401117753046786);
        AssertAlmostEqual(&test, cint.lnD_TE, 117.33853608403136877599);
        AssertAlmostEqual(&test, cint.lnD_TM, 117.87470585363882880168);
    }

    {
        casimir_integrate_drude(&casimir, &cint, 60, 50, 1, 1, 2);

        AssertAlmostEqual(&test, cint.lnA_TE, 317.54496357191706674104);
        AssertAlmostEqual(&test, cint.lnA_TM, 317.54496357191706674104);
        AssertAlmostEqual(&test, cint.lnB_TE, 332.17054025252939692152);
        AssertAlmostEqual(&test, cint.lnB_TM, 332.17054025252939692152);
        AssertAlmostEqual(&test, cint.lnC_TE, 324.76202098855843751046);
        AssertAlmostEqual(&test, cint.lnC_TM, 324.76202098855843751046);
        AssertAlmostEqual(&test, cint.lnD_TE, 324.94434361445935795232);
        AssertAlmostEqual(&test, cint.lnD_TM, 324.94434361445935795232);
    }

    casimir_free(&casimir);

    return test_results(&test, stderr);
}


int test_fresnel()
{
    //const double c = 299792458;
    edouble r_TE, r_TM, T;
    double omegap, gamma_;
    unittest_t test;
    casimir_t casimir;

    unittest_init(&test, "Fresnel", "Test Fresnel coefficients");

    T = 1;
    omegap = 1.32e2;
    gamma_ = 6.9e-1;
    casimir_init(&casimir, 0.5, T);

    casimir_set_omegap_plane(&casimir, omegap);
    casimir_set_gamma_plane(&casimir,  gamma_);

    casimir_rp(&casimir, 1*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.97252954726278);
    AssertAlmostEqual(&test, r_TM, +0.98616846109802);

    casimir_rp(&casimir, 10*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.85446954163643);
    AssertAlmostEqual(&test, r_TM, +0.85579839473205);

    casimir_rp(&casimir, 100*T, 1, &r_TE, &r_TM);
    AssertAlmostEqual(&test, r_TE, -0.24595396364878);
    AssertAlmostEqual(&test, r_TM, +0.24598373253191);

    casimir_free(&casimir);

    return test_results(&test, stderr);
}

int main(int argc, char *argv[])
{
    test_Lambda();
    test_mie_drude();

    test_fresnel();
    test_integration_drude();
    test_Plm();
    test_doublefact();
    test_epsilon();
    test_Xi();
    test_integration_perf();
    test_mie();
    test_besselI();
    test_besselK();
    test_givens();
    test_logdet();
    test_casimirF();
    
    return 0;
}

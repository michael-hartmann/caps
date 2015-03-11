#include "sfunc.h"
#include "libcasimir.h"
#include "unittest.h"

#include "test_mie_drude.h"

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

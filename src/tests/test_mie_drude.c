/* This code was created by test_Mie_drude.py */
#include "sfunc.h"
#include "libcasimir.h"
#include "unittest.h"

#include "test_mie_drude.h"

int test_mie_drude(void)
{
    double lna, lnb;
    double userdata[2];
    sign_t sign_a, sign_b;
    casimir_t *casimir;
    unittest_t test;

    unittest_init(&test, "Mie (Drude)", "Test Mie coefficients for various parameters");

    casimir = casimir_init(1);
    userdata[0] = 50000; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1159277664814);
    AssertAlmostEqual(&test, lnb, -30.8210995465918);

    casimir_lnab(casimir, 0.0001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -52.9186351440535);
    AssertAlmostEqual(&test, lnb, -53.3441410389290);

    casimir_lnab(casimir, 0.0001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -76.3987413462750);
    AssertAlmostEqual(&test, lnb, -76.7144800670882);

    casimir_lnab(casimir, 0.0001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -124.856306648681);
    AssertAlmostEqual(&test, lnb, -125.082715087215);

    casimir_lnab(casimir, 0.0001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -251.521915025219);
    AssertAlmostEqual(&test, lnb, -251.701372948978);

    casimir_lnab(casimir, 0.0001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2856.25013979865);
    AssertAlmostEqual(&test, lnb, -2857.04529806842);

    casimir_lnab(casimir, 0.0001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15828.7460563613);
    AssertAlmostEqual(&test, lnb, -15831.6388012656);

    casimir_lnab(casimir, 0.0001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -33026.9769783904);
    AssertAlmostEqual(&test, lnb, -33031.1694102899);

    casimir_lnab(casimir, 0.0001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -181158.086103576);
    AssertAlmostEqual(&test, lnb, -181165.465810366);

    casimir_lnab(casimir, 0.001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.2081725621256);
    AssertAlmostEqual(&test, lnb, -23.9051185501941);

    casimir_lnab(casimir, 0.001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -41.4057096498743);
    AssertAlmostEqual(&test, lnb, -41.8175064421236);

    casimir_lnab(casimir, 0.001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -60.2806456666204);
    AssertAlmostEqual(&test, lnb, -60.5771921146443);

    casimir_lnab(casimir, 0.001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -99.5278706049169);
    AssertAlmostEqual(&test, lnb, -99.7241218406189);

    casimir_lnab(casimir, 0.001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -203.167628060834);
    AssertAlmostEqual(&test, lnb, -203.289530698612);

    casimir_lnab(casimir, 0.001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2393.43053610562);
    AssertAlmostEqual(&test, lnb, -2393.69433832413);

    casimir_lnab(casimir, 0.001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13523.8583782740);
    AssertAlmostEqual(&test, lnb, -13525.0554629530);

    casimir_lnab(casimir, 0.001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -28419.5042073092);
    AssertAlmostEqual(&test, lnb, -28421.6211392231);

    casimir_lnab(casimir, 0.001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -158129.932588543);
    AssertAlmostEqual(&test, lnb, -158135.021726622);

    casimir_lnab(casimir, 0.01, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.3004246396561);
    AssertAlmostEqual(&test, lnb, -16.9947556693060);

    casimir_lnab(casimir, 0.01, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -29.8927812465623);
    AssertAlmostEqual(&test, lnb, -30.3002537538391);

    casimir_lnab(casimir, 0.01, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -44.1625471338384);
    AssertAlmostEqual(&test, lnb, -44.4530427613503);

    casimir_lnab(casimir, 0.01, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -74.1994324913207);
    AssertAlmostEqual(&test, lnb, -74.3861766420109);

    casimir_lnab(casimir, 0.01, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -154.813339952895);
    AssertAlmostEqual(&test, lnb, -154.917093666718);

    casimir_lnab(casimir, 0.01, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1930.61093229090);
    AssertAlmostEqual(&test, lnb, -1930.70167787637);

    casimir_lnab(casimir, 0.01, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11218.9707001623);
    AssertAlmostEqual(&test, lnb, -11219.3725063412);

    casimir_lnab(casimir, 0.01, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -23812.0314362157);
    AssertAlmostEqual(&test, lnb, -23812.8167219197);

    casimir_lnab(casimir, 0.01, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -135101.779073507);
    AssertAlmostEqual(&test, lnb, -135104.675788057);

    casimir_lnab(casimir, 0.1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.39333077003238);
    AssertAlmostEqual(&test, lnb, -10.0847476859746);

    casimir_lnab(casimir, 0.1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -18.3795611718255);
    AssertAlmostEqual(&test, lnb, -18.7853923166282);

    casimir_lnab(casimir, 0.1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -28.0441629055691);
    AssertAlmostEqual(&test, lnb, -28.3326764195014);

    casimir_lnab(casimir, 0.1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -48.8707871637774);
    AssertAlmostEqual(&test, lnb, -49.0545445714623);

    casimir_lnab(casimir, 0.1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -106.458937369474);
    AssertAlmostEqual(&test, lnb, -106.557030390201);

    casimir_lnab(casimir, 0.1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1467.79131629450);
    AssertAlmostEqual(&test, lnb, -1467.82793326036);

    casimir_lnab(casimir, 0.1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8914.08301960411);
    AssertAlmostEqual(&test, lnb, -8914.21772664695);

    casimir_lnab(casimir, 0.1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19204.5586638984);
    AssertAlmostEqual(&test, lnb, -19204.8243699985);

    casimir_lnab(casimir, 0.1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -112073.625558226);
    AssertAlmostEqual(&test, lnb, -112074.870654882);

    casimir_lnab(casimir, 0.5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -4.57382426354832);
    AssertAlmostEqual(&test, lnb, -5.22460667514423);

    casimir_lnab(casimir, 0.5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -10.3251084858909);
    AssertAlmostEqual(&test, lnb, -10.7236935347369);

    casimir_lnab(casimir, 0.5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.7710907971802);
    AssertAlmostEqual(&test, lnb, -17.0568491249318);

    casimir_lnab(casimir, 0.5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -31.1618950135075);
    AssertAlmostEqual(&test, lnb, -31.3443920783102);

    casimir_lnab(casimir, 0.5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -72.6579374838834);
    AssertAlmostEqual(&test, lnb, -72.7546203748580);

    casimir_lnab(casimir, 0.5, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1144.29399752672);
    AssertAlmostEqual(&test, lnb, -1144.31787378049);

    casimir_lnab(casimir, 0.5, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -7303.03560933398);
    AssertAlmostEqual(&test, lnb, -7303.10694715070);

    casimir_lnab(casimir, 0.5, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15984.0733711409);
    AssertAlmostEqual(&test, lnb, -15984.2128979848);

    casimir_lnab(casimir, 0.5, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -95977.6369899751);
    AssertAlmostEqual(&test, lnb, -95978.3169395785);

    casimir_lnab(casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -2.49510979505747);
    AssertAlmostEqual(&test, lnb, -3.05878951540333);

    casimir_lnab(casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6.83558165055565);
    AssertAlmostEqual(&test, lnb, -7.21439818618779);

    casimir_lnab(casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -11.8970652837077);
    AssertAlmostEqual(&test, lnb, -12.1757682018406);

    casimir_lnab(casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.5214195355106);
    AssertAlmostEqual(&test, lnb, -23.7020382743985);

    casimir_lnab(casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -58.0930865882382);
    AssertAlmostEqual(&test, lnb, -58.1892567188044);

    casimir_lnab(casimir, 1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1004.97048178196);
    AssertAlmostEqual(&test, lnb, -1004.99180214002);

    casimir_lnab(casimir, 1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6609.19509432153);
    AssertAlmostEqual(&test, lnb, -6609.25371144194);

/*
    casimir_lnab(casimir, 1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6609.19509432153);
    AssertAlmostEqual(&test, lnb, -6609.25371144194);
*/

    casimir_lnab(casimir, 1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -89045.4720184510);
    AssertAlmostEqual(&test, lnb, -89046.0306898940);

    casimir_lnab(casimir, 2, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -0.31270026008453);
    AssertAlmostEqual(&test, lnb, -0.69332817781021);

    casimir_lnab(casimir, 2, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -3.26537589276904);
    AssertAlmostEqual(&test, lnb, -3.58335490407145);

    casimir_lnab(casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -6.95596925338256);
    AssertAlmostEqual(&test, lnb, -7.20995773284596);

    casimir_lnab(casimir, 2, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -15.8334440481969);
    AssertAlmostEqual(&test, lnb, -16.0073226755596);

    casimir_lnab(casimir, 2, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -43.5019855279802);
    AssertAlmostEqual(&test, lnb, -43.5970245999405);

    casimir_lnab(casimir, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -865.644168636168);
    AssertAlmostEqual(&test, lnb, -865.663964152055);

    casimir_lnab(casimir, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35401747601);
    AssertAlmostEqual(&test, lnb, -5915.40505044155);

/*
    casimir_lnab(casimir, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35401747601);
    AssertAlmostEqual(&test, lnb, -5915.40505044155);
*/

/*
    casimir_lnab(casimir, 2, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35401747601);
    AssertAlmostEqual(&test, lnb, -5915.40505044155);
*/

    casimir_lnab(casimir, 10, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 8.916918236764934);
    AssertAlmostEqual(&test, lnb, 8.901018834225075);

    casimir_lnab(casimir, 10, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 8.146855595007360);
    AssertAlmostEqual(&test, lnb, 8.109982752630207);

    casimir_lnab(casimir, 10, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 7.017206589890090);
    AssertAlmostEqual(&test, lnb, 6.962263995064897);

    casimir_lnab(casimir, 10, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 3.795475872121321);
    AssertAlmostEqual(&test, lnb, 3.721288716606143);

    casimir_lnab(casimir, 10, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8.60869810168657);
    AssertAlmostEqual(&test, lnb, -8.67954430448111);

    casimir_lnab(casimir, 10, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -542.027827087426);
    AssertAlmostEqual(&test, lnb, -542.046162743147);

    casimir_lnab(casimir, 10, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -4304.28269560390);
    AssertAlmostEqual(&test, lnb, -4304.32668331238);

    casimir_lnab(casimir, 10, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9989.60063139424);
    AssertAlmostEqual(&test, lnb, -9989.68555359466);

/*
    casimir_lnab(casimir, 10, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9989.60063139424);
    AssertAlmostEqual(&test, lnb, -9989.68555359466);
*/

    casimir_lnab(casimir, 100, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99.26284514532585);
    AssertAlmostEqual(&test, lnb, 99.26282592859967);

    casimir_lnab(casimir, 100, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99.18288560962050);
    AssertAlmostEqual(&test, lnb, 99.18282808237047);

    casimir_lnab(casimir, 100, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99.06298604327061);
    AssertAlmostEqual(&test, lnb, 99.06287135575831);

    casimir_lnab(casimir, 100, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 98.70357250964279);
    AssertAlmostEqual(&test, lnb, 98.70328851132725);

    casimir_lnab(casimir, 100, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 97.11129188241470);
    AssertAlmostEqual(&test, lnb, 97.11029286127742);

    casimir_lnab(casimir, 100, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -67.2547109076145);
    AssertAlmostEqual(&test, lnb, -67.2690863267283);

    casimir_lnab(casimir, 100, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1996.92601006013);
    AssertAlmostEqual(&test, lnb, -1996.96801688693);

    casimir_lnab(casimir, 100, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5380.89156960429);
    AssertAlmostEqual(&test, lnb, -5380.97288401401);

/*
    casimir_lnab(casimir, 100, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5380.89156960429);
    AssertAlmostEqual(&test, lnb, -5380.97288401401);
*/

    casimir_lnab(casimir, 1000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 999.2628356653451);
    AssertAlmostEqual(&test, lnb, 999.2628353292358);

    casimir_lnab(casimir, 1000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 999.2548360227806);
    AssertAlmostEqual(&test, lnb, 999.2548350144614);

    casimir_lnab(casimir, 1000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 999.2428365989168);
    AssertAlmostEqual(&test, lnb, 999.2428345823050);

    casimir_lnab(casimir, 1000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 999.2068386151943);
    AssertAlmostEqual(&test, lnb, 999.2068335738637);

    casimir_lnab(casimir, 1000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 999.0468528000582);
    AssertAlmostEqual(&test, lnb, 999.0468343184185);

    casimir_lnab(casimir, 1000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 979.1348549582963);
    AssertAlmostEqual(&test, lnb, 979.1331939722822);

    casimir_lnab(casimir, 1000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 531.2375811370925);
    AssertAlmostEqual(&test, lnb, 531.2085415530016);

    casimir_lnab(casimir, 1000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -653.357239755226);
    AssertAlmostEqual(&test, lnb, -653.429548449970);

/*
    casimir_lnab(casimir, 1000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -653.357239755226);
    AssertAlmostEqual(&test, lnb, -653.429548449970);
*/

    casimir_lnab(casimir, 10000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9998.909053003338);
    AssertAlmostEqual(&test, lnb, 9998.909052971941);

    casimir_lnab(casimir, 10000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9998.908253034757);
    AssertAlmostEqual(&test, lnb, 9998.908252940565);

    casimir_lnab(casimir, 10000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9998.907053081925);
    AssertAlmostEqual(&test, lnb, 9998.907052893541);

    casimir_lnab(casimir, 10000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9998.903453223716);
    AssertAlmostEqual(&test, lnb, 9998.903452752757);

    casimir_lnab(casimir, 10000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9998.887453859127);
    AssertAlmostEqual(&test, lnb, 9998.887452132278);

    casimir_lnab(casimir, 10000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9996.889600239351);
    AssertAlmostEqual(&test, lnb, 9996.889441716154);

    casimir_lnab(casimir, 10000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9948.853112329297);
    AssertAlmostEqual(&test, lnb, 9948.849200139500);

/*
    casimir_lnab(casimir, 10000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9948.853112329297);
    AssertAlmostEqual(&test, lnb, 9948.849200139500);
*/

/*
    casimir_lnab(casimir, 10000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9948.853112329297);
    AssertAlmostEqual(&test, lnb, 9948.849200139500);
*/

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99996.41953292556);
    AssertAlmostEqual(&test, lnb, 99996.41953292413);

    casimir_lnab(casimir, 100000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99996.41945292699);
    AssertAlmostEqual(&test, lnb, 99996.41945292270);

    casimir_lnab(casimir, 100000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99996.41933292914);
    AssertAlmostEqual(&test, lnb, 99996.41933292055);

    casimir_lnab(casimir, 100000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99996.41897293558);
    AssertAlmostEqual(&test, lnb, 99996.41897291411);

    casimir_lnab(casimir, 100000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99996.41737296421);
    AssertAlmostEqual(&test, lnb, 99996.41737288550);

    casimir_lnab(casimir, 100000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99996.21757660643);
    AssertAlmostEqual(&test, lnb, 99996.21756937926);

    casimir_lnab(casimir, 100000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99991.40970436927);
    AssertAlmostEqual(&test, lnb, 99991.40952513683);

    casimir_lnab(casimir, 100000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99976.40059879121);
    AssertAlmostEqual(&test, lnb, 99976.39988277039);

/*
    casimir_lnab(casimir, 100000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99976.40059879121);
    AssertAlmostEqual(&test, lnb, 99976.39988277039);
*/


    userdata[0] = 5000; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1159277666167);
    AssertAlmostEqual(&test, lnb, -30.9314748316560);

    casimir_lnab(casimir, 0.0001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -52.9186351441289);
    AssertAlmostEqual(&test, lnb, -53.5278601163339);

    casimir_lnab(casimir, 0.0001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -76.3987413463280);
    AssertAlmostEqual(&test, lnb, -76.9711875129116);

    casimir_lnab(casimir, 0.0001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -124.856306648715);
    AssertAlmostEqual(&test, lnb, -125.483804575374);

    casimir_lnab(casimir, 0.0001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -251.521915025238);
    AssertAlmostEqual(&test, lnb, -252.448926555279);

    casimir_lnab(casimir, 0.0001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2856.25013979865);
    AssertAlmostEqual(&test, lnb, -2860.46862287780);

    casimir_lnab(casimir, 0.0001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15828.7460563613);
    AssertAlmostEqual(&test, lnb, -15836.1311495378);

    casimir_lnab(casimir, 0.0001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -33026.9769783904);
    AssertAlmostEqual(&test, lnb, -33035.7444413673);

    casimir_lnab(casimir, 0.0001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -181158.086103576);
    AssertAlmostEqual(&test, lnb, -181170.069745086);

    casimir_lnab(casimir, 0.001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.2081725663976);
    AssertAlmostEqual(&test, lnb, -23.9395260024476);

    casimir_lnab(casimir, 0.001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -41.4057096522484);
    AssertAlmostEqual(&test, lnb, -41.8748451820842);

    casimir_lnab(casimir, 0.001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -60.2806456682832);
    AssertAlmostEqual(&test, lnb, -60.6574516407742);

    casimir_lnab(casimir, 0.001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -99.5278706059638);
    AssertAlmostEqual(&test, lnb, -99.8501747361564);

    casimir_lnab(casimir, 0.001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -203.167628061383);
    AssertAlmostEqual(&test, lnb, -203.529594147250);

    casimir_lnab(casimir, 0.001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2393.43053610571);
    AssertAlmostEqual(&test, lnb, -2395.56900072445);

    casimir_lnab(casimir, 0.001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13523.8583782741);
    AssertAlmostEqual(&test, lnb, -13528.9528438759);

    casimir_lnab(casimir, 0.001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -28419.5042073093);
    AssertAlmostEqual(&test, lnb, -28425.9727785427);

    casimir_lnab(casimir, 0.001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -158129.932588543);
    AssertAlmostEqual(&test, lnb, -158139.614656763);

    casimir_lnab(casimir, 0.01, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.3004247753295);
    AssertAlmostEqual(&test, lnb, -17.0056336398328);

    casimir_lnab(casimir, 0.01, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -29.8927813219402);
    AssertAlmostEqual(&test, lnb, -30.3183833489943);

    casimir_lnab(casimir, 0.01, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -44.1625471866059);
    AssertAlmostEqual(&test, lnb, -44.4784237027239);

    casimir_lnab(casimir, 0.01, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -74.1994325244942);
    AssertAlmostEqual(&test, lnb, -74.4260587842684);

    casimir_lnab(casimir, 0.01, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -154.813339970179);
    AssertAlmostEqual(&test, lnb, -154.993213904547);

    casimir_lnab(casimir, 0.01, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1930.61093229286);
    AssertAlmostEqual(&test, lnb, -1931.40978440289);

    casimir_lnab(casimir, 0.01, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11218.9707001632);
    AssertAlmostEqual(&test, lnb, -11221.8722647737);

    casimir_lnab(casimir, 0.01, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -23812.0314362166);
    AssertAlmostEqual(&test, lnb, -23816.2334265737);

    casimir_lnab(casimir, 0.01, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -135101.779073507);
    AssertAlmostEqual(&test, lnb, -135109.168618412);

    casimir_lnab(casimir, 0.1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.39333524001735);
    AssertAlmostEqual(&test, lnb, -10.0883356916960);

    casimir_lnab(casimir, 0.1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -18.3795636585684);
    AssertAlmostEqual(&test, lnb, -18.7913680089945);

    casimir_lnab(casimir, 0.1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -28.0441646465997);
    AssertAlmostEqual(&test, lnb, -28.3410413128944);

    casimir_lnab(casimir, 0.1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -48.8707882582465);
    AssertAlmostEqual(&test, lnb, -49.0676884294258);

    casimir_lnab(casimir, 0.1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -106.458937939385);
    AssertAlmostEqual(&test, lnb, -106.582121770255);

    casimir_lnab(casimir, 0.1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1467.79131635449);
    AssertAlmostEqual(&test, lnb, -1468.06731345759);

    casimir_lnab(casimir, 0.1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8914.08301961869);
    AssertAlmostEqual(&test, lnb, -8915.33145958767);

    casimir_lnab(casimir, 0.1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19204.5586639088);
    AssertAlmostEqual(&test, lnb, -19206.7502909838);

    casimir_lnab(casimir, 0.1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -112073.625558235);
    AssertAlmostEqual(&test, lnb, -112078.807903467);

    casimir_lnab(casimir, 0.5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -4.57388074309667);
    AssertAlmostEqual(&test, lnb, -5.22651694760275);

    casimir_lnab(casimir, 0.5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -10.3251407256382);
    AssertAlmostEqual(&test, lnb, -10.7268307454283);

    casimir_lnab(casimir, 0.5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.7711134598172);
    AssertAlmostEqual(&test, lnb, -17.0612276374752);

    casimir_lnab(casimir, 0.5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -31.1619092870310);
    AssertAlmostEqual(&test, lnb, -31.3512609158360);

    casimir_lnab(casimir, 0.5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -72.6579449218891);
    AssertAlmostEqual(&test, lnb, -72.7677232861149);

    casimir_lnab(casimir, 0.5, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1144.29399830436);
    AssertAlmostEqual(&test, lnb, -1144.44314074805);

    casimir_lnab(casimir, 0.5, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -7303.03560949981);
    AssertAlmostEqual(&test, lnb, -7303.71813928759);

    casimir_lnab(casimir, 0.5, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15984.0733712375);
    AssertAlmostEqual(&test, lnb, -15985.3688142142);

    casimir_lnab(casimir, 0.5, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -95977.6369900358);
    AssertAlmostEqual(&test, lnb, -95981.5491894513);

    casimir_lnab(casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -2.49528098479370);
    AssertAlmostEqual(&test, lnb, -3.06042744468423);

    casimir_lnab(casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6.83568480157874);
    AssertAlmostEqual(&test, lnb, -7.21700284943256);

    casimir_lnab(casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -11.8971386258412);
    AssertAlmostEqual(&test, lnb, -12.1793722814230);

    casimir_lnab(casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.5214659969069);
    AssertAlmostEqual(&test, lnb, -23.7076641505091);

    casimir_lnab(casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -58.0931108593250);
    AssertAlmostEqual(&test, lnb, -58.1999636524308);

    casimir_lnab(casimir, 1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1004.97048431945);
    AssertAlmostEqual(&test, lnb, -1005.09410654636);

    casimir_lnab(casimir, 1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6609.19509485234);
    AssertAlmostEqual(&test, lnb, -6609.75618739131);

    casimir_lnab(casimir, 1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -14597.0857694536);
    AssertAlmostEqual(&test, lnb, -14598.1658608761);

    casimir_lnab(casimir, 1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -89045.4720186149);
    AssertAlmostEqual(&test, lnb, -89048.9980790879);

    casimir_lnab(casimir, 2, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -0.31319508331980);
    AssertAlmostEqual(&test, lnb, -0.69495743008139);

    casimir_lnab(casimir, 2, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -3.26570916730829);
    AssertAlmostEqual(&test, lnb, -3.58574798151290);

    casimir_lnab(casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -6.95621477668315);
    AssertAlmostEqual(&test, lnb, -7.21317680512191);

    casimir_lnab(casimir, 2, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -15.8336028998429);
    AssertAlmostEqual(&test, lnb, -16.0122559318776);

    casimir_lnab(casimir, 2, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -43.5020693147017);
    AssertAlmostEqual(&test, lnb, -43.6063283805850);

    casimir_lnab(casimir, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -865.644177422043);
    AssertAlmostEqual(&test, lnb, -865.752575201070);

    casimir_lnab(casimir, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35401929584);
    AssertAlmostEqual(&test, lnb, -5915.84173208558);

    casimir_lnab(casimir, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13210.0978871072);
    AssertAlmostEqual(&test, lnb, -13211.0438356659);

    casimir_lnab(casimir, 2, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -82113.3069911911);
    AssertAlmostEqual(&test, lnb, -82116.5641612067);

    casimir_lnab(casimir, 10, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 8.913293144845530);
    AssertAlmostEqual(&test, lnb, 8.897085997314026);

    casimir_lnab(casimir, 10, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 8.143478093533186);
    AssertAlmostEqual(&test, lnb, 8.105760419947395);

    casimir_lnab(casimir, 10, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 7.014115949677902);
    AssertAlmostEqual(&test, lnb, 6.957647786531710);

    casimir_lnab(casimir, 10, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 3.792937199440454);
    AssertAlmostEqual(&test, lnb, 3.715665136663746);

    casimir_lnab(casimir, 10, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8.61032250594057);
    AssertAlmostEqual(&test, lnb, -8.68833259961737);

    casimir_lnab(casimir, 10, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -542.028014890855);
    AssertAlmostEqual(&test, lnb, -542.122140858635);

    casimir_lnab(casimir, 10, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -4304.28273423467);
    AssertAlmostEqual(&test, lnb, -4304.70171098947);

    casimir_lnab(casimir, 10, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9989.60065203381);
    AssertAlmostEqual(&test, lnb, -9990.41837677595);

    casimir_lnab(casimir, 10, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -66017.3160384090);
    AssertAlmostEqual(&test, lnb, -66020.2894892743);

    casimir_lnab(casimir, 100, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99.22668277162390);
    AssertAlmostEqual(&test, lnb, 99.22663459916018);

    casimir_lnab(casimir, 100, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99.14675213365063);
    AssertAlmostEqual(&test, lnb, 99.14660780884361);

    casimir_lnab(casimir, 100, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99.02689578392091);
    AssertAlmostEqual(&test, lnb, 99.02660770963230);

    casimir_lnab(casimir, 100, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 98.66761097301383);
    AssertAlmostEqual(&test, lnb, 98.66689505906935);

    casimir_lnab(casimir, 100, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 97.07588618589847);
    AssertAlmostEqual(&test, lnb, 97.07332803331118);

    casimir_lnab(casimir, 100, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -67.2708373996868);
    AssertAlmostEqual(&test, lnb, -67.3502943083965);

    casimir_lnab(casimir, 100, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1996.92968567549);
    AssertAlmostEqual(&test, lnb, -1997.32935004681);

    casimir_lnab(casimir, 100, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5380.89353102557);
    AssertAlmostEqual(&test, lnb, -5381.67758181789);

    casimir_lnab(casimir, 100, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -42988.9159436673);
    AssertAlmostEqual(&test, lnb, -42991.8118455282);

    casimir_lnab(casimir, 1000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 998.9052781144618);
    AssertAlmostEqual(&test, lnb, 998.9052749578570);

    casimir_lnab(casimir, 1000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 998.8972812923338);
    AssertAlmostEqual(&test, lnb, 998.8972718225989);

    casimir_lnab(casimir, 1000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 998.8852860989925);
    AssertAlmostEqual(&test, lnb, 998.8852671597611);

    casimir_lnab(casimir, 1000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 998.8493008058834);
    AssertAlmostEqual(&test, lnb, 998.8492534595933);

    casimir_lnab(casimir, 1000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 998.6893713761701);
    AssertAlmostEqual(&test, lnb, 998.6891978022412);

    casimir_lnab(casimir, 1000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 978.7842000791890);
    AssertAlmostEqual(&test, lnb, 978.7685868171832);

    casimir_lnab(casimir, 1000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 530.9813315784118);
    AssertAlmostEqual(&test, lnb, 530.7061604040235);

    casimir_lnab(casimir, 1000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -653.528359431779);
    AssertAlmostEqual(&test, lnb, -654.207922744695);

    casimir_lnab(casimir, 1000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19936.1311856588);
    AssertAlmostEqual(&test, lnb, -19938.9427394680);

    casimir_lnab(casimir, 10000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9996.419092501966);
    AssertAlmostEqual(&test, lnb, 9996.419092358814);

    casimir_lnab(casimir, 10000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9996.418292645139);
    AssertAlmostEqual(&test, lnb, 9996.418292215684);

    casimir_lnab(casimir, 10000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9996.417092859937);
    AssertAlmostEqual(&test, lnb, 9996.417092001029);

    casimir_lnab(casimir, 10000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9996.413493504621);
    AssertAlmostEqual(&test, lnb, 9996.413491357353);

    casimir_lnab(casimir, 10000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9996.397496375094);
    AssertAlmostEqual(&test, lnb, 9996.397488501799);

    casimir_lnab(casimir, 10000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9994.399921688987);
    AssertAlmostEqual(&test, lnb, 9994.399199037270);

    casimir_lnab(casimir, 10000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9946.370045958873);
    AssertAlmostEqual(&test, lnb, 9946.352276085007);

/*
    casimir_lnab(casimir, 10000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9946.370045958873);
    AssertAlmostEqual(&test, lnb, 9946.352276085007);
*/

    casimir_lnab(casimir, 10000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 7.523163845262640);
    AssertAlmostEqual(&test, lnb, -1.50463276905252);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99991.92779509478);
    AssertAlmostEqual(&test, lnb, 99991.92779509318);

    casimir_lnab(casimir, 100000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99991.92771509638);
    AssertAlmostEqual(&test, lnb, 99991.92771509159);

    casimir_lnab(casimir, 100000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99991.92759509878);
    AssertAlmostEqual(&test, lnb, 99991.92759508919);

    casimir_lnab(casimir, 100000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99991.92723510597);
    AssertAlmostEqual(&test, lnb, 99991.92723508200);

    casimir_lnab(casimir, 100000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99991.92563513794);
    AssertAlmostEqual(&test, lnb, 99991.92563505005);

    casimir_lnab(casimir, 100000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99991.72583919705);
    AssertAlmostEqual(&test, lnb, 99991.72583112692);

    casimir_lnab(casimir, 100000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99986.91797698919);
    AssertAlmostEqual(&test, lnb, 99986.91777685325);

    casimir_lnab(casimir, 100000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99971.90890268711);
    AssertAlmostEqual(&test, lnb, 99971.90810318194);

/*
    casimir_lnab(casimir, 100000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99971.90890268711);
    AssertAlmostEqual(&test, lnb, 99971.90810318194);
*/


    userdata[0] = 500; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1159277683113);
    AssertAlmostEqual(&test, lnb, -32.1420093429172);

    casimir_lnab(casimir, 0.0001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -52.9186351453174);
    AssertAlmostEqual(&test, lnb, -55.2876247644402);

    casimir_lnab(casimir, 0.0001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -76.3987413473553);
    AssertAlmostEqual(&test, lnb, -79.1452856302028);

    casimir_lnab(casimir, 0.0001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -124.856306649626);
    AssertAlmostEqual(&test, lnb, -128.241309584109);

    casimir_lnab(casimir, 0.0001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -251.521915026078);
    AssertAlmostEqual(&test, lnb, -255.988191933668);

    casimir_lnab(casimir, 0.0001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2856.25013979945);
    AssertAlmostEqual(&test, lnb, -2865.04442281066);

    casimir_lnab(casimir, 0.0001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15828.7460563621);
    AssertAlmostEqual(&test, lnb, -15840.7350908866);

    casimir_lnab(casimir, 0.0001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -33026.9769783912);
    AssertAlmostEqual(&test, lnb, -33040.3493032076);

    casimir_lnab(casimir, 0.0001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -181158.086103577);
    AssertAlmostEqual(&test, lnb, -181174.674902906);

    casimir_lnab(casimir, 0.001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.2081726099791);
    AssertAlmostEqual(&test, lnb, -24.3042106448209);

    casimir_lnab(casimir, 0.001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -41.4057096773917);
    AssertAlmostEqual(&test, lnb, -42.4732386852593);

    casimir_lnab(casimir, 0.001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -60.2806456868102);
    AssertAlmostEqual(&test, lnb, -61.4774174506440);

    casimir_lnab(casimir, 0.001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -99.5278706191849);
    AssertAlmostEqual(&test, lnb, -101.070315021779);

    casimir_lnab(casimir, 0.001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -203.167628071068);
    AssertAlmostEqual(&test, lnb, -205.515194391509);

    casimir_lnab(casimir, 0.001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2393.43053611368);
    AssertAlmostEqual(&test, lnb, -2399.92585407154);

    casimir_lnab(casimir, 0.001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13523.8583782821);
    AssertAlmostEqual(&test, lnb, -13533.5458389058);

    casimir_lnab(casimir, 0.001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -28419.5042073172);
    AssertAlmostEqual(&test, lnb, -28430.5748745480);

    casimir_lnab(casimir, 0.001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -158129.932588551);
    AssertAlmostEqual(&test, lnb, -158144.219703408);

    casimir_lnab(casimir, 0.01, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.3004261345743);
    AssertAlmostEqual(&test, lnb, -17.1165666300263);

    casimir_lnab(casimir, 0.01, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -29.8927820798969);
    AssertAlmostEqual(&test, lnb, -30.5030268892021);

    casimir_lnab(casimir, 0.01, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -44.1625477201133);
    AssertAlmostEqual(&test, lnb, -44.7364174095809);

    casimir_lnab(casimir, 0.01, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -74.1994328653210);
    AssertAlmostEqual(&test, lnb, -74.8291342845191);

    casimir_lnab(casimir, 0.01, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -154.813340159805);
    AssertAlmostEqual(&test, lnb, -155.744296943635);

    casimir_lnab(casimir, 0.01, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1930.61093237432);
    AssertAlmostEqual(&test, lnb, -1934.83898105033);

    casimir_lnab(casimir, 0.01, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11218.9707002433);
    AssertAlmostEqual(&test, lnb, -11226.3656315191);

    casimir_lnab(casimir, 0.01, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -23812.0314362966);
    AssertAlmostEqual(&test, lnb, -23820.8087464763);

    casimir_lnab(casimir, 0.01, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -135101.779073587);
    AssertAlmostEqual(&test, lnb, -135113.772565230);

    casimir_lnab(casimir, 0.1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.39337994825345);
    AssertAlmostEqual(&test, lnb, -10.1244506787708);

    casimir_lnab(casimir, 0.1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -18.3795885405448);
    AssertAlmostEqual(&test, lnb, -18.8515087407516);

    casimir_lnab(casimir, 0.1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -28.0441820773457);
    AssertAlmostEqual(&test, lnb, -28.4252108384302);

    casimir_lnab(casimir, 0.1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -48.8707992350692);
    AssertAlmostEqual(&test, lnb, -49.1998660553876);

    casimir_lnab(casimir, 0.1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -106.458943699630);
    AssertAlmostEqual(&test, lnb, -106.833780959834);

    casimir_lnab(casimir, 0.1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1467.79131739075);
    AssertAlmostEqual(&test, lnb, -1470.00479648707);

    casimir_lnab(casimir, 0.1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8914.08302049356);
    AssertAlmostEqual(&test, lnb, -8919.27069815676);

    casimir_lnab(casimir, 0.1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19204.5586647807);
    AssertAlmostEqual(&test, lnb, -19211.1212669261);

    casimir_lnab(casimir, 0.1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -112073.625559106);
    AssertAlmostEqual(&test, lnb, -112083.401925895);

    casimir_lnab(casimir, 0.5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -4.57444551253063);
    AssertAlmostEqual(&test, lnb, -5.24568267933457);

    casimir_lnab(casimir, 0.5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -10.3254631627194);
    AssertAlmostEqual(&test, lnb, -10.7583080012426);

    casimir_lnab(casimir, 0.5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.7713401540298);
    AssertAlmostEqual(&test, lnb, -17.1051580714901);

    casimir_lnab(casimir, 0.5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -31.1620521347026);
    AssertAlmostEqual(&test, lnb, -31.4201669563367);

    casimir_lnab(casimir, 0.5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -72.6580195186850);
    AssertAlmostEqual(&test, lnb, -72.8990745304989);

    casimir_lnab(casimir, 0.5, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1144.29400795078);
    AssertAlmostEqual(&test, lnb, -1145.60644637030);

    casimir_lnab(casimir, 0.5, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -7303.03561556671);
    AssertAlmostEqual(&test, lnb, -7306.95299413723);

    casimir_lnab(casimir, 0.5, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15984.0733771972);
    AssertAlmostEqual(&test, lnb, -15989.3441704457);

    casimir_lnab(casimir, 0.5, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -95977.6369959764);
    AssertAlmostEqual(&test, lnb, -95986.1143773895);

    casimir_lnab(casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -2.49699249332108);
    AssertAlmostEqual(&test, lnb, -3.07684718196641);

    casimir_lnab(casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6.83671630736972);
    AssertAlmostEqual(&test, lnb, -7.24311856154593);

    casimir_lnab(casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -11.8978721586954);
    AssertAlmostEqual(&test, lnb, -12.2155095845009);

    casimir_lnab(casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.5219308445905);
    AssertAlmostEqual(&test, lnb, -23.7640695842190);

    casimir_lnab(casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -58.0933540392666);
    AssertAlmostEqual(&test, lnb, -58.3072630359163);

    casimir_lnab(casimir, 1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1004.97051389396);
    AssertAlmostEqual(&test, lnb, -1006.06615734291);

    casimir_lnab(casimir, 1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6609.19511124777);
    AssertAlmostEqual(&test, lnb, -6612.72624751850);

    casimir_lnab(casimir, 1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -14597.0857853842);
    AssertAlmostEqual(&test, lnb, -14601.9561974422);

    casimir_lnab(casimir, 1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -89045.4720344566);
    AssertAlmostEqual(&test, lnb, -89053.5441488068);

    casimir_lnab(casimir, 2, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -0.31814086475595);
    AssertAlmostEqual(&test, lnb, -0.71127728227742);

    casimir_lnab(casimir, 2, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -3.26904117557170);
    AssertAlmostEqual(&test, lnb, -3.60972780442086);

    casimir_lnab(casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -6.95866995108623);
    AssertAlmostEqual(&test, lnb, -7.24543779846586);

    casimir_lnab(casimir, 2, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -15.8351919071527);
    AssertAlmostEqual(&test, lnb, -16.0616977738451);

    casimir_lnab(casimir, 2, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -43.5029083732862);
    AssertAlmostEqual(&test, lnb, -43.6995445152265);

    casimir_lnab(casimir, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -865.644276371479);
    AssertAlmostEqual(&test, lnb, -866.604819255952);

    casimir_lnab(casimir, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35406915864);
    AssertAlmostEqual(&test, lnb, -5918.61619580173);

    casimir_lnab(casimir, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13210.0979350272);
    AssertAlmostEqual(&test, lnb, -13214.6856885286);

    casimir_lnab(casimir, 2, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -82113.3070387163);
    AssertAlmostEqual(&test, lnb, -82121.0916469729);

    casimir_lnab(casimir, 10, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 8.877050631636407);
    AssertAlmostEqual(&test, lnb, 8.857754273142655);

    casimir_lnab(casimir, 10, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 8.109716453097378);
    AssertAlmostEqual(&test, lnb, 8.063523213690927);

    casimir_lnab(casimir, 10, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 6.983224459722323);
    AssertAlmostEqual(&test, lnb, 6.911459009642089);

    casimir_lnab(casimir, 10, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 3.767560446950503);
    AssertAlmostEqual(&test, lnb, 3.659374951615984);

    casimir_lnab(casimir, 10, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8.62657535308101);
    AssertAlmostEqual(&test, lnb, -8.77632612980932);

    casimir_lnab(casimir, 10, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -542.030068683149);
    AssertAlmostEqual(&test, lnb, -542.860252350036);

    casimir_lnab(casimir, 10, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -4304.28366639578);
    AssertAlmostEqual(&test, lnb, -4307.26112768316);

    casimir_lnab(casimir, 10, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9989.60153394377);
    AssertAlmostEqual(&test, lnb, -9993.88558642200);

    casimir_lnab(casimir, 10, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -66017.3169093592);
    AssertAlmostEqual(&test, lnb, -66024.7908333921);

    casimir_lnab(casimir, 100, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 98.86768543510334);
    AssertAlmostEqual(&test, lnb, 98.86735283678003);

    casimir_lnab(casimir, 100, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 98.78803862349888);
    AssertAlmostEqual(&test, lnb, 98.78704173999789);

    casimir_lnab(casimir, 100, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 98.66860667056510);
    AssertAlmostEqual(&test, lnb, 98.66661562881138);

    casimir_lnab(casimir, 100, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 98.31058547610765);
    AssertAlmostEqual(&test, lnb, 98.30562814793536);

    casimir_lnab(casimir, 100, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 96.72430910023380);
    AssertAlmostEqual(&test, lnb, 96.70645220520135);

    casimir_lnab(casimir, 100, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -67.4422875019250);
    AssertAlmostEqual(&test, lnb, -68.1356429160195);

    casimir_lnab(casimir, 100, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1997.01236368717);
    AssertAlmostEqual(&test, lnb, -1999.83602203563);

    casimir_lnab(casimir, 100, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5380.97146759965);
    AssertAlmostEqual(&test, lnb, -5385.09577282544);

    casimir_lnab(casimir, 100, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -42988.9928466536);
    AssertAlmostEqual(&test, lnb, -42996.3046862140);

    casimir_lnab(casimir, 1000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 996.4146950237423);
    AssertAlmostEqual(&test, lnb, 996.4146806698208);

    casimir_lnab(casimir, 1000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 996.4067093984783);
    AssertAlmostEqual(&test, lnb, 996.4066663373365);

    casimir_lnab(casimir, 1000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 996.3947309994143);
    AssertAlmostEqual(&test, lnb, 996.3946448789989);

    casimir_lnab(casimir, 1000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 996.3587960818032);
    AssertAlmostEqual(&test, lnb, 996.3585807947760);

    casimir_lnab(casimir, 1000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 996.1990904101846);
    AssertAlmostEqual(&test, lnb, 996.1983012526403);

    casimir_lnab(casimir, 1000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 976.3202507348934);
    AssertAlmostEqual(&test, lnb, 976.2502878334325);

    casimir_lnab(casimir, 1000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 528.7636054933661);
    AssertAlmostEqual(&test, lnb, 527.7371502707579);

    casimir_lnab(casimir, 1000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -655.652466360678);
    AssertAlmostEqual(&test, lnb, -657.749061677611);

    casimir_lnab(casimir, 1000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19938.2488438479);
    AssertAlmostEqual(&test, lnb, -19943.4356723452);

    casimir_lnab(casimir, 10000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9991.927345291142);
    AssertAlmostEqual(&test, lnb, 9991.927345131294);

    casimir_lnab(casimir, 10000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9991.926545451011);
    AssertAlmostEqual(&test, lnb, 9991.926544971467);

    casimir_lnab(casimir, 10000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9991.925345690856);
    AssertAlmostEqual(&test, lnb, 9991.925344731766);

    casimir_lnab(casimir, 10000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9991.921746410675);
    AssertAlmostEqual(&test, lnb, 9991.921744012953);

    casimir_lnab(casimir, 10000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9991.905749615082);
    AssertAlmostEqual(&test, lnb, 9991.905740823462);

    casimir_lnab(casimir, 10000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9989.908216582450);
    AssertAlmostEqual(&test, lnb, 9989.907409674053);

    casimir_lnab(casimir, 10000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9941.879315993920);
    AssertAlmostEqual(&test, lnb, 9941.859492804306);

/*
    casimir_lnab(casimir, 10000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9941.879315993920);
    AssertAlmostEqual(&test, lnb, 9941.859492804306);
*/

    casimir_lnab(casimir, 10000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -3.95718418260814);
    AssertAlmostEqual(&test, lnb, 9.193306218910946);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99987.32386122631);
    AssertAlmostEqual(&test, lnb, 99987.32386122471);

    casimir_lnab(casimir, 100000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99987.32378122791);
    AssertAlmostEqual(&test, lnb, 99987.32378122311);

    casimir_lnab(casimir, 100000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99987.32366123031);
    AssertAlmostEqual(&test, lnb, 99987.32366122071);

    casimir_lnab(casimir, 100000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99987.32330123751);
    AssertAlmostEqual(&test, lnb, 99987.32330121351);

    casimir_lnab(casimir, 100000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99987.32170126952);
    AssertAlmostEqual(&test, lnb, 99987.32170118152);

    casimir_lnab(casimir, 100000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99987.12190533356);
    AssertAlmostEqual(&test, lnb, 99987.12189725346);

    casimir_lnab(casimir, 100000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99982.31404324445);
    AssertAlmostEqual(&test, lnb, 99982.31384286102);

    casimir_lnab(casimir, 100000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99967.30496931262);
    AssertAlmostEqual(&test, lnb, 99967.30416881909);

/*
    casimir_lnab(casimir, 100000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99967.30496931262);
    AssertAlmostEqual(&test, lnb, 99967.30416881909);
*/


    userdata[0] = 100; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1159277972141);
    AssertAlmostEqual(&test, lnb, -34.9270201732242);

    casimir_lnab(casimir, 0.0001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -52.9186351693441);
    AssertAlmostEqual(&test, lnb, -58.2768815405117);

    casimir_lnab(casimir, 0.0001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -76.3987413697664);
    AssertAlmostEqual(&test, lnb, -82.2224199985171);

    casimir_lnab(casimir, 0.0001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -124.856306670750);
    AssertAlmostEqual(&test, lnb, -131.390891360055);

    casimir_lnab(casimir, 0.0001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -251.521915046240);
    AssertAlmostEqual(&test, lnb, -259.184587886664);

    casimir_lnab(casimir, 0.0001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2856.25013981875);
    AssertAlmostEqual(&test, lnb, -2868.26300750585);

    casimir_lnab(casimir, 0.0001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15828.7460563814);
    AssertAlmostEqual(&test, lnb, -15843.9539547844);

    casimir_lnab(casimir, 0.0001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -33026.9769784104);
    AssertAlmostEqual(&test, lnb, -33043.5681760417);

    casimir_lnab(casimir, 0.0001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -181158.086103596);
    AssertAlmostEqual(&test, lnb, -181177.893778611);

    casimir_lnab(casimir, 0.001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.2081729321297);
    AssertAlmostEqual(&test, lnb, -25.9053113498129);

    casimir_lnab(casimir, 0.001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -41.4057099285141);
    AssertAlmostEqual(&test, lnb, -44.5556283354711);

    casimir_lnab(casimir, 0.001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -60.2806459158140);
    AssertAlmostEqual(&test, lnb, -63.8586224213267);

    casimir_lnab(casimir, 0.001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -99.5278708319838);
    AssertAlmostEqual(&test, lnb, -103.787449474881);

    casimir_lnab(casimir, 0.001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -203.167628273049);
    AssertAlmostEqual(&test, lnb, -208.537117172957);

    casimir_lnab(casimir, 0.001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2393.43053630684);
    AssertAlmostEqual(&test, lnb, -2403.14182730444);

    casimir_lnab(casimir, 0.001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13523.8583784744);
    AssertAlmostEqual(&test, lnb, -13536.7645955776);

    casimir_lnab(casimir, 0.001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -28419.5042075095);
    AssertAlmostEqual(&test, lnb, -28433.7937204932);

    casimir_lnab(casimir, 0.001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -158129.932588743);
    AssertAlmostEqual(&test, lnb, -158147.438578035);

    casimir_lnab(casimir, 0.01, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.3004325418233);
    AssertAlmostEqual(&test, lnb, -17.6510452758808);

    casimir_lnab(casimir, 0.01, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -29.8927860127745);
    AssertAlmostEqual(&test, lnb, -31.3528002751322);

    casimir_lnab(casimir, 0.01, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -44.1625508166686);
    AssertAlmostEqual(&test, lnb, -45.8602748664506);

    casimir_lnab(casimir, 0.01, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -74.1994353412410);
    AssertAlmostEqual(&test, lnb, -76.3899231164739);

    casimir_lnab(casimir, 0.01, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -154.813342271202);
    AssertAlmostEqual(&test, lnb, -157.968499588439);

    casimir_lnab(casimir, 0.01, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1930.61093432325);
    AssertAlmostEqual(&test, lnb, -1938.02966838043);

    casimir_lnab(casimir, 0.01, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11218.9707021844);
    AssertAlmostEqual(&test, lnb, -11229.5833274459);

    casimir_lnab(casimir, 0.01, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -23812.0314382368);
    AssertAlmostEqual(&test, lnb, -23824.0273262314);

    casimir_lnab(casimir, 0.01, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -135101.779075527);
    AssertAlmostEqual(&test, lnb, -135116.991429181);

    casimir_lnab(casimir, 0.1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.39357976416710);
    AssertAlmostEqual(&test, lnb, -10.2900797505772);

    casimir_lnab(casimir, 0.1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -18.3797009869866);
    AssertAlmostEqual(&test, lnb, -19.1262041141441);

    casimir_lnab(casimir, 0.1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -28.0442621322247);
    AssertAlmostEqual(&test, lnb, -28.8073456176828);

    casimir_lnab(casimir, 0.1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -48.8708519895825);
    AssertAlmostEqual(&test, lnb, -49.7896685902581);

    casimir_lnab(casimir, 0.1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -106.458976235942);
    AssertAlmostEqual(&test, lnb, -107.885501022547);

    casimir_lnab(casimir, 0.1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1467.79133864142);
    AssertAlmostEqual(&test, lnb, -1473.00256929417);

    casimir_lnab(casimir, 0.1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8914.08304163449);
    AssertAlmostEqual(&test, lnb, -8922.47882412415);

    casimir_lnab(casimir, 0.1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19204.5586859110);
    AssertAlmostEqual(&test, lnb, -19214.3374296195);

    casimir_lnab(casimir, 0.1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -112073.625580228);
    AssertAlmostEqual(&test, lnb, -112086.620692704);

    casimir_lnab(casimir, 0.5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -4.57695800789616);
    AssertAlmostEqual(&test, lnb, -5.33224446588576);

    casimir_lnab(casimir, 0.5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -10.3269022903522);
    AssertAlmostEqual(&test, lnb, -10.9003888853375);

    casimir_lnab(casimir, 0.5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.7723565100645);
    AssertAlmostEqual(&test, lnb, -17.3031618855088);

    casimir_lnab(casimir, 0.5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -31.1627009843122);
    AssertAlmostEqual(&test, lnb, -31.7293308476897);

    casimir_lnab(casimir, 0.5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -72.6583771451707);
    AssertAlmostEqual(&test, lnb, -73.4770136447419);

    casimir_lnab(casimir, 0.5, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1144.29415418778);
    AssertAlmostEqual(&test, lnb, -1148.23431826657);

    casimir_lnab(casimir, 0.5, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -7303.03575970417);
    AssertAlmostEqual(&test, lnb, -7310.13333408076);

    casimir_lnab(casimir, 0.5, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15984.0735212583);
    AssertAlmostEqual(&test, lnb, -15992.5531559335);

    casimir_lnab(casimir, 0.5, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -95977.6371399795);
    AssertAlmostEqual(&test, lnb, -95989.3328535840);

    casimir_lnab(casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -2.50459611473380);
    AssertAlmostEqual(&test, lnb, -3.15070799856665);

    casimir_lnab(casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6.84131135926316);
    AssertAlmostEqual(&test, lnb, -7.36063418403314);

    casimir_lnab(casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -11.9011502314995);
    AssertAlmostEqual(&test, lnb, -12.3780032345205);

    casimir_lnab(casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.5240266656151);
    AssertAlmostEqual(&test, lnb, -24.0169750380513);

    casimir_lnab(casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -58.0944917417081);
    AssertAlmostEqual(&test, lnb, -58.7821335164782);

    casimir_lnab(casimir, 1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1004.97090722753);
    AssertAlmostEqual(&test, lnb, -1008.52418038617);

    casimir_lnab(casimir, 1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6609.19549557944);
    AssertAlmostEqual(&test, lnb, -6615.88818046454);

    casimir_lnab(casimir, 1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -14597.0861694981);
    AssertAlmostEqual(&test, lnb, -14605.1602965571);

    casimir_lnab(casimir, 1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -89045.4724184152);
    AssertAlmostEqual(&test, lnb, -89056.7624252834);

    casimir_lnab(casimir, 2, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -0.34007648324806);
    AssertAlmostEqual(&test, lnb, -0.78439978700532);

    casimir_lnab(casimir, 2, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -3.28385673312994);
    AssertAlmostEqual(&test, lnb, -3.71732457658362);

    casimir_lnab(casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -6.96961692329500);
    AssertAlmostEqual(&test, lnb, -7.39019623647199);

    casimir_lnab(casimir, 2, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -15.8423264980465);
    AssertAlmostEqual(&test, lnb, -16.2831686082018);

    casimir_lnab(casimir, 2, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -43.5067837188489);
    AssertAlmostEqual(&test, lnb, -44.1131246158703);

    casimir_lnab(casimir, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -865.645467259695);
    AssertAlmostEqual(&test, lnb, -868.928574908340);

    casimir_lnab(casimir, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35522173441);
    AssertAlmostEqual(&test, lnb, -5921.76024915621);

    casimir_lnab(casimir, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13210.0990868934);
    AssertAlmostEqual(&test, lnb, -13217.8849396028);

    casimir_lnab(casimir, 2, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -82113.3081901130);
    AssertAlmostEqual(&test, lnb, -82124.3097237963);

    casimir_lnab(casimir, 10, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 8.716419076400930);
    AssertAlmostEqual(&test, lnb, 8.683207775500249);

    casimir_lnab(casimir, 10, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 7.960136213309104);
    AssertAlmostEqual(&test, lnb, 7.875909428714485);

    casimir_lnab(casimir, 10, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 6.846303821711579);
    AssertAlmostEqual(&test, lnb, 6.706128436925515);

    casimir_lnab(casimir, 10, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 3.654683332029487);
    AssertAlmostEqual(&test, lnb, 3.408990143110445);

    casimir_lnab(casimir, 10, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8.70035656294111);
    AssertAlmostEqual(&test, lnb, -9.16600047608385);

    casimir_lnab(casimir, 10, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -542.051974098285);
    AssertAlmostEqual(&test, lnb, -545.029955076633);

    casimir_lnab(casimir, 10, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -4304.30457272327);
    AssertAlmostEqual(&test, lnb, -4310.38006103160);

    casimir_lnab(casimir, 10, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9989.62242621709);
    AssertAlmostEqual(&test, lnb, -9997.07785351025);

    casimir_lnab(casimir, 10, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -66017.3377932929);
    AssertAlmostEqual(&test, lnb, -66028.0086199275);

    casimir_lnab(casimir, 100, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 97.49764674281867);
    AssertAlmostEqual(&test, lnb, 97.49648014059801);

    casimir_lnab(casimir, 100, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 97.41883103983317);
    AssertAlmostEqual(&test, lnb, 97.41533561581357);

    casimir_lnab(casimir, 100, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 97.30063927921260);
    AssertAlmostEqual(&test, lnb, 97.29366152357149);

    casimir_lnab(casimir, 100, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 96.94629274436243);
    AssertAlmostEqual(&test, lnb, 96.92894556189040);

    casimir_lnab(casimir, 100, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 95.37555765888724);
    AssertAlmostEqual(&test, lnb, 95.31347288199910);

    casimir_lnab(casimir, 100, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -68.4898431494614);
    AssertAlmostEqual(&test, lnb, -70.3714009188785);

    casimir_lnab(casimir, 100, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1998.03808034542);
    AssertAlmostEqual(&test, lnb, -2002.94774483492);

    casimir_lnab(casimir, 100, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5381.99861433345);
    AssertAlmostEqual(&test, lnb, -5388.28579294146);

    casimir_lnab(casimir, 100, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -42990.0204317418);
    AssertAlmostEqual(&test, lnb, -42999.5223758947);

    casimir_lnab(casimir, 1000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 993.3054203632264);
    AssertAlmostEqual(&test, lnb, 993.3054043948058);

    casimir_lnab(casimir, 1000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 993.2974363523421);
    AssertAlmostEqual(&test, lnb, 993.2973884478464);

    casimir_lnab(casimir, 1000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 993.2854603745790);
    AssertAlmostEqual(&test, lnb, 993.2853645678857);

    casimir_lnab(casimir, 1000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 993.2495327189377);
    AssertAlmostEqual(&test, lnb, 993.2492932194388);

    casimir_lnab(casimir, 1000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 993.0898592876880);
    AssertAlmostEqual(&test, lnb, 993.0889814035693);

    casimir_lnab(casimir, 1000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 973.2146326818015);
    AssertAlmostEqual(&test, lnb, 973.1370825785836);

    casimir_lnab(casimir, 1000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 525.6733945799852);
    AssertAlmostEqual(&test, lnb, 524.5753532956497);

    casimir_lnab(casimir, 1000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -658.749587036286);
    AssertAlmostEqual(&test, lnb, -660.944466981496);

    casimir_lnab(casimir, 1000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19941.3545193342);
    AssertAlmostEqual(&test, lnb, -19946.6533630820);

    casimir_lnab(casimir, 10000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9988.709668178240);
    AssertAlmostEqual(&test, lnb, 9988.709668018200);

    casimir_lnab(casimir, 10000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9988.708868338301);
    AssertAlmostEqual(&test, lnb, 9988.708867858181);

    casimir_lnab(casimir, 10000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9988.707668578433);
    AssertAlmostEqual(&test, lnb, 9988.707667618193);

    casimir_lnab(casimir, 10000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9988.704069299115);
    AssertAlmostEqual(&test, lnb, 9988.704066898518);

    casimir_lnab(casimir, 10000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9988.688072507355);
    AssertAlmostEqual(&test, lnb, 9988.688063705193);

    casimir_lnab(casimir, 10000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9986.690539952868);
    AssertAlmostEqual(&test, lnb, 9986.689732077249);

    casimir_lnab(casimir, 10000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9938.661650532715);
    AssertAlmostEqual(&test, lnb, 9938.641803806149);

/*
    casimir_lnab(casimir, 10000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9938.661650532715);
    AssertAlmostEqual(&test, lnb, 9938.641803806149);
*/

    casimir_lnab(casimir, 10000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 5.579178307646773);
    AssertAlmostEqual(&test, lnb, -3.33155787723610);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99984.10499740121);
    AssertAlmostEqual(&test, lnb, 99984.10499739961);

    casimir_lnab(casimir, 100000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99984.10491740281);
    AssertAlmostEqual(&test, lnb, 99984.10491739801);

    casimir_lnab(casimir, 100000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99984.10479740521);
    AssertAlmostEqual(&test, lnb, 99984.10479739561);

    casimir_lnab(casimir, 100000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99984.10443741241);
    AssertAlmostEqual(&test, lnb, 99984.10443738841);

    casimir_lnab(casimir, 100000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99984.10283744442);
    AssertAlmostEqual(&test, lnb, 99984.10283735641);

    casimir_lnab(casimir, 100000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99983.90304150851);
    AssertAlmostEqual(&test, lnb, 99983.90303342830);

    casimir_lnab(casimir, 100000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99979.09517942055);
    AssertAlmostEqual(&test, lnb, 99979.09497903471);

    casimir_lnab(casimir, 100000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99964.08610549231);
    AssertAlmostEqual(&test, lnb, 99964.08530498918);

/*
    casimir_lnab(casimir, 100000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99964.08610549231);
    AssertAlmostEqual(&test, lnb, 99964.08530498918);
*/


    userdata[0] = 50; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1159278872271);
    AssertAlmostEqual(&test, lnb, -36.2957462502624);

    casimir_lnab(casimir, 0.0001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -52.9186352443524);
    AssertAlmostEqual(&test, lnb, -59.6549100570271);

    casimir_lnab(casimir, 0.0001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -76.3987414397737);
    AssertAlmostEqual(&test, lnb, -83.6038684855173);

    casimir_lnab(casimir, 0.0001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -124.856306736757);
    AssertAlmostEqual(&test, lnb, -132.774918687818);

    casimir_lnab(casimir, 0.0001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -251.521915109247);
    AssertAlmostEqual(&test, lnb, -260.570168608970);

    casimir_lnab(casimir, 0.0001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2856.25013987905);
    AssertAlmostEqual(&test, lnb, -2869.64929276714);

    casimir_lnab(casimir, 0.0001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15828.7460564414);
    AssertAlmostEqual(&test, lnb, -15845.3402487728);

    casimir_lnab(casimir, 0.0001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -33026.9769784704);
    AssertAlmostEqual(&test, lnb, -33044.9544703094);

    casimir_lnab(casimir, 0.0001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -181158.086103656);
    AssertAlmostEqual(&test, lnb, -181179.280072968);

    casimir_lnab(casimir, 0.001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.2081738365482);
    AssertAlmostEqual(&test, lnb, -27.1379721711891);

    casimir_lnab(casimir, 0.001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -41.4057106800817);
    AssertAlmostEqual(&test, lnb, -45.8648377460133);

    casimir_lnab(casimir, 0.001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -60.2806466168046);
    AssertAlmostEqual(&test, lnb, -65.1985261552242);

    casimir_lnab(casimir, 0.001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -99.5278714927093);
    AssertAlmostEqual(&test, lnb, -105.151571872961);

    casimir_lnab(casimir, 0.001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -203.167628903686);
    AssertAlmostEqual(&test, lnb, -209.916332774643);

    casimir_lnab(casimir, 0.001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2393.43053691044);
    AssertAlmostEqual(&test, lnb, -2404.52803075831);

    casimir_lnab(casimir, 0.001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13523.8583790756);
    AssertAlmostEqual(&test, lnb, -13538.1508862149);

    casimir_lnab(casimir, 0.001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -28419.5042081104);
    AssertAlmostEqual(&test, lnb, -28435.1800139206);

    casimir_lnab(casimir, 0.001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -158129.932589343);
    AssertAlmostEqual(&test, lnb, -158148.824872359);

    casimir_lnab(casimir, 0.01, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.3004431995967);
    AssertAlmostEqual(&test, lnb, -18.3328195475876);

    casimir_lnab(casimir, 0.01, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -29.8927940618809);
    AssertAlmostEqual(&test, lnb, -32.2695613866148);

    casimir_lnab(casimir, 0.01, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -44.1625580832500);
    AssertAlmostEqual(&test, lnb, -46.9176111220742);

    casimir_lnab(casimir, 0.01, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -74.1994420595206);
    AssertAlmostEqual(&test, lnb, -77.5936086179280);

    casimir_lnab(casimir, 0.01, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -154.813348640259);
    AssertAlmostEqual(&test, lnb, -159.289241250861);

    casimir_lnab(casimir, 0.01, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1930.61094041352);
    AssertAlmostEqual(&test, lnb, -1939.41506266905);

    casimir_lnab(casimir, 0.01, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11218.9707082504);
    AssertAlmostEqual(&test, lnb, -11230.9695849017);

    casimir_lnab(casimir, 0.01, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -23812.0314442998);
    AssertAlmostEqual(&test, lnb, -23825.4136113382);

    casimir_lnab(casimir, 0.01, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -135101.779081587);
    AssertAlmostEqual(&test, lnb, -135118.377723171);

    casimir_lnab(casimir, 0.1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.39383696152267);
    AssertAlmostEqual(&test, lnb, -10.5081193417979);

    casimir_lnab(casimir, 0.1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -18.3798535515411);
    AssertAlmostEqual(&test, lnb, -19.4795860116321);

    casimir_lnab(casimir, 0.1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -28.0443782597089);
    AssertAlmostEqual(&test, lnb, -29.2839323144160);

    casimir_lnab(casimir, 0.1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -48.8709404989010);
    AssertAlmostEqual(&test, lnb, -50.4717998007807);

    casimir_lnab(casimir, 0.1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -106.459048710175);
    AssertAlmostEqual(&test, lnb, -108.884408615586);

    casimir_lnab(casimir, 0.1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1467.79140496873);
    AssertAlmostEqual(&test, lnb, -1474.38067358387);

    casimir_lnab(casimir, 0.1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8914.08310769686);
    AssertAlmostEqual(&test, lnb, -8923.86477975421);

    casimir_lnab(casimir, 0.1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19204.5587519404);
    AssertAlmostEqual(&test, lnb, -19215.7236390174);

    casimir_lnab(casimir, 0.1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -112073.625646231);
    AssertAlmostEqual(&test, lnb, -112088.006983658);

    casimir_lnab(casimir, 0.5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -4.58011874350317);
    AssertAlmostEqual(&test, lnb, -5.44356211096847);

    casimir_lnab(casimir, 0.5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -10.3287403421100);
    AssertAlmostEqual(&test, lnb, -11.0822410458808);

    casimir_lnab(casimir, 0.5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.7736822403988);
    AssertAlmostEqual(&test, lnb, -17.5545988424090);

    casimir_lnab(casimir, 0.5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -31.1635968795748);
    AssertAlmostEqual(&test, lnb, -32.1130907264185);

    casimir_lnab(casimir, 0.5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -72.6589684994314);
    AssertAlmostEqual(&test, lnb, -74.1370392739044);

    casimir_lnab(casimir, 0.5, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1144.29460635368);
    AssertAlmostEqual(&test, lnb, -1149.59136240170);

    casimir_lnab(casimir, 0.5, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -7303.03620998525);
    AssertAlmostEqual(&test, lnb, -7311.51838782470);

    casimir_lnab(casimir, 0.5, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15984.0739713144);
    AssertAlmostEqual(&test, lnb, -15993.9391388508);

    casimir_lnab(casimir, 0.5, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -95977.6375898558);
    AssertAlmostEqual(&test, lnb, -95990.7191354528);

    casimir_lnab(casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -2.51411997321138);
    AssertAlmostEqual(&test, lnb, -3.24502541928073);

    casimir_lnab(casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6.84713061301701);
    AssertAlmostEqual(&test, lnb, -7.51044079882530);

    casimir_lnab(casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -11.9053623981107);
    AssertAlmostEqual(&test, lnb, -12.5841718669796);

    casimir_lnab(casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.5268290312713);
    AssertAlmostEqual(&test, lnb, -24.3331165807176);

    casimir_lnab(casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -58.0962387617265);
    AssertAlmostEqual(&test, lnb, -59.3418302853163);

    casimir_lnab(casimir, 1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1004.97211251580);
    AssertAlmostEqual(&test, lnb, -1009.86736239736);

    casimir_lnab(casimir, 1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6609.19669557871);
    AssertAlmostEqual(&test, lnb, -6617.27261533794);

    casimir_lnab(casimir, 1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -14597.0873688981);
    AssertAlmostEqual(&test, lnb, -14606.5461238434);

    casimir_lnab(casimir, 1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -89045.4736173363);
    AssertAlmostEqual(&test, lnb, -89058.1487009062);

    casimir_lnab(casimir, 2, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -0.36743139891157);
    AssertAlmostEqual(&test, lnb, -0.87709660532933);

    casimir_lnab(casimir, 2, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -3.30249992908649);
    AssertAlmostEqual(&test, lnb, -3.85385129508083);

    casimir_lnab(casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -6.98355246226421);
    AssertAlmostEqual(&test, lnb, -7.57342791765374);

    casimir_lnab(casimir, 2, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -15.8516947627583);
    AssertAlmostEqual(&test, lnb, -16.5606199634532);

    casimir_lnab(casimir, 2, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -43.5124752122749);
    AssertAlmostEqual(&test, lnb, -44.6084884264939);

    casimir_lnab(casimir, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -865.649076887820);
    AssertAlmostEqual(&test, lnb, -870.258364290130);

    casimir_lnab(casimir, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5915.35881454682);
    AssertAlmostEqual(&test, lnb, -5923.14406611960);

    casimir_lnab(casimir, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13210.1026779173);
    AssertAlmostEqual(&test, lnb, -13219.2706113189);

    casimir_lnab(casimir, 2, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -82113.3117797070);
    AssertAlmostEqual(&test, lnb, -82125.6959931732);

    casimir_lnab(casimir, 10, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 8.517798831887984);
    AssertAlmostEqual(&test, lnb, 8.467113763998793);

    casimir_lnab(casimir, 10, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 7.775008941552725);
    AssertAlmostEqual(&test, lnb, 7.643492011204608);

    casimir_lnab(casimir, 10, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 6.676209861057171);
    AssertAlmostEqual(&test, lnb, 6.451806821284042);

    casimir_lnab(casimir, 10, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 3.512076374474462);
    AssertAlmostEqual(&test, lnb, 3.100045102492270);

    casimir_lnab(casimir, 10, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8.80139337111546);
    AssertAlmostEqual(&test, lnb, -9.63357351540353);

    casimir_lnab(casimir, 10, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -542.114864285450);
    AssertAlmostEqual(&test, lnb, -546.341193340921);

    casimir_lnab(casimir, 10, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -4304.36720893009);
    AssertAlmostEqual(&test, lnb, -4311.76298122633);

    casimir_lnab(casimir, 10, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9989.68503481926);
    AssertAlmostEqual(&test, lnb, -9998.46329907013);

    casimir_lnab(casimir, 10, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -66017.4003788555);
    AssertAlmostEqual(&test, lnb, -66029.3948802196);

    casimir_lnab(casimir, 100, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 96.37141733239700);
    AssertAlmostEqual(&test, lnb, 96.36994266274481);

    casimir_lnab(casimir, 100, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 96.29290781418510);
    AssertAlmostEqual(&test, lnb, 96.28849045293043);

    casimir_lnab(casimir, 100, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 96.17517113379122);
    AssertAlmostEqual(&test, lnb, 96.16635625460540);

    casimir_lnab(casimir, 100, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 95.82216025930016);
    AssertAlmostEqual(&test, lnb, 95.80027011815425);

    casimir_lnab(casimir, 100, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 94.25686446555718);
    AssertAlmostEqual(&test, lnb, 94.17888593299894);

    casimir_lnab(casimir, 100, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -69.5748732255245);
    AssertAlmostEqual(&test, lnb, -71.6909562647765);

    casimir_lnab(casimir, 100, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1999.13803884414);
    AssertAlmostEqual(&test, lnb, -2004.33040129920);

    casimir_lnab(casimir, 100, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5383.09926755156);
    AssertAlmostEqual(&test, lnb, -5389.67116522321);

    casimir_lnab(casimir, 100, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -42991.1212601541);
    AssertAlmostEqual(&test, lnb, -43000.9086331548);

    casimir_lnab(casimir, 1000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 991.9228548442002);
    AssertAlmostEqual(&test, lnb, 991.9228388161588);

    casimir_lnab(casimir, 1000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 991.9148708929320);
    AssertAlmostEqual(&test, lnb, 991.9148228095795);

    casimir_lnab(casimir, 1000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 991.9028950045822);
    AssertAlmostEqual(&test, lnb, 991.9027988401926);

    casimir_lnab(casimir, 1000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 991.8669676171037);
    AssertAlmostEqual(&test, lnb, 991.8667272234929);

    casimir_lnab(casimir, 1000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 991.7072953762883);
    AssertAlmostEqual(&test, lnb, 991.7064142158564);

    casimir_lnab(casimir, 1000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 971.8322010154538);
    AssertAlmostEqual(&test, lnb, 971.7543723379680);

    casimir_lnab(casimir, 1000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 524.2914511812079);
    AssertAlmostEqual(&test, lnb, 523.1909231925077);

    casimir_lnab(casimir, 1000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -660.131818104445);
    AssertAlmostEqual(&test, lnb, -662.330014458217);

    casimir_lnab(casimir, 1000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19942.7370611326);
    AssertAlmostEqual(&test, lnb, -19948.0396203752);

    casimir_lnab(casimir, 10000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9987.323411311616);
    AssertAlmostEqual(&test, lnb, 9987.323411151570);

    casimir_lnab(casimir, 10000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9987.322611471683);
    AssertAlmostEqual(&test, lnb, 9987.322610991545);

    casimir_lnab(casimir, 10000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9987.321411711824);
    AssertAlmostEqual(&test, lnb, 9987.321410751548);

    casimir_lnab(casimir, 10000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9987.317812432533);
    AssertAlmostEqual(&test, lnb, 9987.317810031846);

    casimir_lnab(casimir, 10000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9987.301815640893);
    AssertAlmostEqual(&test, lnb, 9987.301806838401);

    casimir_lnab(casimir, 10000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9985.304283101376);
    AssertAlmostEqual(&test, lnb, 9985.303475195475);

    casimir_lnab(casimir, 10000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9937.275394030872);
    AssertAlmostEqual(&test, lnb, 9937.255546567425);

/*
    casimir_lnab(casimir, 10000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9937.275394030872);
    AssertAlmostEqual(&test, lnb, 9937.255546567425);
*/

    casimir_lnab(casimir, 10000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 1.027363254709066);
    AssertAlmostEqual(&test, lnb, -2.74437493911335);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99982.71870341508);
    AssertAlmostEqual(&test, lnb, 99982.71870341348);

    casimir_lnab(casimir, 100000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99982.71862341668);
    AssertAlmostEqual(&test, lnb, 99982.71862341188);

    casimir_lnab(casimir, 100000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99982.71850341908);
    AssertAlmostEqual(&test, lnb, 99982.71850340948);

    casimir_lnab(casimir, 100000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99982.71814342628);
    AssertAlmostEqual(&test, lnb, 99982.71814340228);

    casimir_lnab(casimir, 100000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99982.71654345829);
    AssertAlmostEqual(&test, lnb, 99982.71654337029);

    casimir_lnab(casimir, 100000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99982.51674752239);
    AssertAlmostEqual(&test, lnb, 99982.51673944218);

    casimir_lnab(casimir, 100000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99977.70888543446);
    AssertAlmostEqual(&test, lnb, 99977.70868504855);

    casimir_lnab(casimir, 100000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99962.69981150634);
    AssertAlmostEqual(&test, lnb, 99962.69901100291);

/*
    casimir_lnab(casimir, 100000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99962.69981150634);
    AssertAlmostEqual(&test, lnb, 99962.69901100291);
*/


    userdata[0] = 1; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1162277522158);
    AssertAlmostEqual(&test, lnb, -44.1138622558291);

    casimir_lnab(casimir, 0.0001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -52.9188851380916);
    AssertAlmostEqual(&test, lnb, -67.4761841527582);

    casimir_lnab(casimir, 0.0001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -76.3989746758744);
    AssertAlmostEqual(&test, lnb, -91.4262935225587);

    casimir_lnab(casimir, 0.0001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -124.856526646546);
    AssertAlmostEqual(&test, lnb, -140.598207866872);

    casimir_lnab(casimir, 0.0001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -251.522125024187);
    AssertAlmostEqual(&test, lnb, -268.393976682016);

    casimir_lnab(casimir, 0.0001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2856.25034079854);
    AssertAlmostEqual(&test, lnb, -2877.47333574591);

    casimir_lnab(casimir, 0.0001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15828.7462565613);
    AssertAlmostEqual(&test, lnb, -15853.1642946594);

    casimir_lnab(casimir, 0.0001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -33026.9771784904);
    AssertAlmostEqual(&test, lnb, -33052.7785162891);

    casimir_lnab(casimir, 0.0001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -181158.086303596);
    AssertAlmostEqual(&test, lnb, -181187.104118978);

    casimir_lnab(casimir, 0.001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -23.2111711360542);
    AssertAlmostEqual(&test, lnb, -34.9044427585422);

    casimir_lnab(casimir, 0.001, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -41.4082090532443);
    AssertAlmostEqual(&test, lnb, -53.6615830615667);

    casimir_lnab(casimir, 0.001, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -60.2829786258421);
    AssertAlmostEqual(&test, lnb, -73.0065180991476);

    casimir_lnab(casimir, 0.001, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -99.5300703905526);
    AssertAlmostEqual(&test, lnb, -112.968088964047);

    casimir_lnab(casimir, 0.001, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -203.169727956525);
    AssertAlmostEqual(&test, lnb, -217.738004987780);

    casimir_lnab(casimir, 0.001, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2393.43254609425);
    AssertAlmostEqual(&test, lnb, -2412.35204647615);

    casimir_lnab(casimir, 0.001, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13523.8603802707);
    AssertAlmostEqual(&test, lnb, -13545.9749309849);

    casimir_lnab(casimir, 0.001, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -28419.5062083068);
    AssertAlmostEqual(&test, lnb, -28443.0040596203);

    casimir_lnab(casimir, 0.001, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -158129.934588741);
    AssertAlmostEqual(&test, lnb, -158156.648918357);

    casimir_lnab(casimir, 0.01, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -16.3302815601576);
    AssertAlmostEqual(&test, lnb, -25.7032616012340);

    casimir_lnab(casimir, 0.01, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -29.9177206115840);
    AssertAlmostEqual(&test, lnb, -39.8551194878818);

    casimir_lnab(casimir, 0.01, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -44.1858419741631);
    AssertAlmostEqual(&test, lnb, -54.5948437382827);

    casimir_lnab(casimir, 0.01, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -74.2214099144107);
    AssertAlmostEqual(&test, lnb, -85.3460440037785);

    casimir_lnab(casimir, 0.01, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -154.834328354609);
    AssertAlmostEqual(&test, lnb, -167.090091260455);

    casimir_lnab(casimir, 0.01, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1930.63102997493);
    AssertAlmostEqual(&test, lnb, -1947.23880850916);

    casimir_lnab(casimir, 0.01, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11218.9907186487);
    AssertAlmostEqual(&test, lnb, -11238.7936186152);

    casimir_lnab(casimir, 0.01, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -23812.0514448022);
    AssertAlmostEqual(&test, lnb, -23833.2376542655);

    casimir_lnab(casimir, 0.01, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -135101.799074173);
    AssertAlmostEqual(&test, lnb, -135126.201769059);

    casimir_lnab(casimir, 0.1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.67878208232218);
    AssertAlmostEqual(&test, lnb, -16.5798519455862);

    casimir_lnab(casimir, 0.1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -18.6227275360744);
    AssertAlmostEqual(&test, lnb, -26.1255930387970);

    casimir_lnab(casimir, 0.1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -28.2727571518085);
    AssertAlmostEqual(&test, lnb, -36.2598235873359);

    casimir_lnab(casimir, 0.1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -49.0875705350822);
    AssertAlmostEqual(&test, lnb, -57.8004631530715);

    casimir_lnab(casimir, 0.1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -106.666782541498);
    AssertAlmostEqual(&test, lnb, -116.518557207703);

    casimir_lnab(casimir, 0.1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1467.99106860221);
    AssertAlmostEqual(&test, lnb, -1482.20196848313);

    casimir_lnab(casimir, 0.1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8914.28205078219);
    AssertAlmostEqual(&test, lnb, -8931.68871286174);

    casimir_lnab(casimir, 0.1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19204.7576049186);
    AssertAlmostEqual(&test, lnb, -19223.5476567161);

    casimir_lnab(casimir, 0.1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -112073.824427117);
    AssertAlmostEqual(&test, lnb, -112095.831028534);

    casimir_lnab(casimir, 0.5, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -5.74294949798069);
    AssertAlmostEqual(&test, lnb, -10.4492855249682);

    casimir_lnab(casimir, 0.5, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -11.3828337554416);
    AssertAlmostEqual(&test, lnb, -16.7751195520322);

    casimir_lnab(casimir, 0.5, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -17.7843822441432);
    AssertAlmostEqual(&test, lnb, -23.6905863441216);

    casimir_lnab(casimir, 0.5, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -32.1374162071812);
    AssertAlmostEqual(&test, lnb, -38.7940961538267);

    casimir_lnab(casimir, 0.5, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -73.6041123377358);
    AssertAlmostEqual(&test, lnb, -81.4189088438124);

    casimir_lnab(casimir, 0.5, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1145.21328792964);
    AssertAlmostEqual(&test, lnb, -1157.40537535869);

    casimir_lnab(casimir, 0.5, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -7303.95250004258);
    AssertAlmostEqual(&test, lnb, -7319.34201995007);

    casimir_lnab(casimir, 0.5, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -15984.9899618627);
    AssertAlmostEqual(&test, lnb, -16001.7630810563);

    casimir_lnab(casimir, 0.5, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -95978.5533407051);
    AssertAlmostEqual(&test, lnb, -95998.5431773012);

    casimir_lnab(casimir, 1, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -4.39969882097945);
    AssertAlmostEqual(&test, lnb, -7.94097323503990);

    casimir_lnab(casimir, 1, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -8.62854457248949);
    AssertAlmostEqual(&test, lnb, -12.8848436031172);

    casimir_lnab(casimir, 1, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -13.6357394568173);
    AssertAlmostEqual(&test, lnb, -18.4171046262538);

    casimir_lnab(casimir, 1, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -25.2106780718991);
    AssertAlmostEqual(&test, lnb, -30.7520156030649);

    casimir_lnab(casimir, 1, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -59.7427994363058);
    AssertAlmostEqual(&test, lnb, -66.4500124609031);

    casimir_lnab(casimir, 1, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -1006.58392593180);
    AssertAlmostEqual(&test, lnb, -1017.67641434949);

    casimir_lnab(casimir, 1, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6610.80533245353);
    AssertAlmostEqual(&test, lnb, -6625.09604061716);

    casimir_lnab(casimir, 1, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -14598.6956071107);
    AssertAlmostEqual(&test, lnb, -14614.3700141523);

    casimir_lnab(casimir, 1, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -89047.0815363604);
    AssertAlmostEqual(&test, lnb, -89065.9727406733);

    casimir_lnab(casimir, 2, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -3.16898884259990);
    AssertAlmostEqual(&test, lnb, -5.46911331758430);

    casimir_lnab(casimir, 2, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -6.03847703487923);
    AssertAlmostEqual(&test, lnb, -9.04981331920652);

    casimir_lnab(casimir, 2, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -9.67511705100451);
    AssertAlmostEqual(&test, lnb, -13.2103513213590);

    casimir_lnab(casimir, 2, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -18.4957808872793);
    AssertAlmostEqual(&test, lnb, -22.7902853540574);

    casimir_lnab(casimir, 2, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -46.1156645868873);
    AssertAlmostEqual(&test, lnb, -51.5764198496199);

    casimir_lnab(casimir, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -868.213773640563);
    AssertAlmostEqual(&test, lnb, -878.062492523017);

    casimir_lnab(casimir, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5917.91989142655);
    AssertAlmostEqual(&test, lnb, -5930.96728461884);

    casimir_lnab(casimir, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -13212.6632973424);
    AssertAlmostEqual(&test, lnb, -13227.0944497354);

    casimir_lnab(casimir, 2, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -82115.8720323551);
    AssertAlmostEqual(&test, lnb, -82133.5200308591);

    casimir_lnab(casimir, 10, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 2.918113498974539);
    AssertAlmostEqual(&test, lnb, 2.706862240002826);

    casimir_lnab(casimir, 10, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 2.259513945914223);
    AssertAlmostEqual(&test, lnb, 1.751399424093723);

    casimir_lnab(casimir, 10, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 1.225692517811374);
    AssertAlmostEqual(&test, lnb, 0.409098726431880);

    casimir_lnab(casimir, 10, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -1.86534420143846);
    AssertAlmostEqual(&test, lnb, -3.24021545161202);

    casimir_lnab(casimir, 10, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -14.1294372134960);
    AssertAlmostEqual(&test, lnb, -16.5244490839536);

    casimir_lnab(casimir, 10, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -547.432138200871);
    AssertAlmostEqual(&test, lnb, -554.138282425717);

    casimir_lnab(casimir, 10, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -4309.68189891716);
    AssertAlmostEqual(&test, lnb, -4319.58589917357);

    casimir_lnab(casimir, 10, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -9994.99930211747);
    AssertAlmostEqual(&test, lnb, -10006.2870620205);

    casimir_lnab(casimir, 10, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -66022.7142914056);
    AssertAlmostEqual(&test, lnb, -66037.2189148783);

    casimir_lnab(casimir, 100, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 88.66104485692112);
    AssertAlmostEqual(&test, lnb, 88.65939629910903);

    casimir_lnab(casimir, 100, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 88.58270788698239);
    AssertAlmostEqual(&test, lnb, 88.57777047668709);

    casimir_lnab(casimir, 100, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 88.46522704665477);
    AssertAlmostEqual(&test, lnb, 88.45537687780671);

    casimir_lnab(casimir, 100, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 88.11296275996437);
    AssertAlmostEqual(&test, lnb, 88.08851979045763);

    casimir_lnab(casimir, 100, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 86.55063860793009);
    AssertAlmostEqual(&test, lnb, 86.46383549020551);

    casimir_lnab(casimir, 100, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -77.2732930843417);
    AssertAlmostEqual(&test, lnb, -79.4912280283300);

    casimir_lnab(casimir, 100, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -2006.84502099627);
    AssertAlmostEqual(&test, lnb, -2012.15323075469);

    casimir_lnab(casimir, 100, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -5390.80659878320);
    AssertAlmostEqual(&test, lnb, -5397.49490371415);

    casimir_lnab(casimir, 100, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -42998.8286788531);
    AssertAlmostEqual(&test, lnb, -43008.7326668031);

    casimir_lnab(casimir, 1000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 984.1000559271849);
    AssertAlmostEqual(&test, lnb, 984.1000398791289);

    casimir_lnab(casimir, 1000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 984.0920719959296);
    AssertAlmostEqual(&test, lnb, 984.0920238525354);

    casimir_lnab(casimir, 1000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 984.0800961375956);
    AssertAlmostEqual(&test, lnb, 984.0799998531284);

    casimir_lnab(casimir, 1000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 984.0441688401383);
    AssertAlmostEqual(&test, lnb, 984.0439281463771);

    casimir_lnab(casimir, 1000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 983.8844969989456);
    AssertAlmostEqual(&test, lnb, 983.8836147386687);

    casimir_lnab(casimir, 1000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 964.0094470129908);
    AssertAlmostEqual(&test, lnb, 963.9315248454191);

    casimir_lnab(casimir, 1000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 516.4688598564322);
    AssertAlmostEqual(&test, lnb, 515.3674995106236);

    casimir_lnab(casimir, 1000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -667.954506427466);
    AssertAlmostEqual(&test, lnb, -670.153811421426);

    casimir_lnab(casimir, 1000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, -19950.5598536550);
    AssertAlmostEqual(&test, lnb, -19955.8636540345);

    casimir_lnab(casimir, 10000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9979.499377794394);
    AssertAlmostEqual(&test, lnb, 9979.499377634346);

    casimir_lnab(casimir, 10000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9979.498577954463);
    AssertAlmostEqual(&test, lnb, 9979.498577474319);

    casimir_lnab(casimir, 10000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9979.497378194607);
    AssertAlmostEqual(&test, lnb, 9979.497377234320);

    casimir_lnab(casimir, 10000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 9979.493778915325);
    AssertAlmostEqual(&test, lnb, 9979.493776514608);

    casimir_lnab(casimir, 10000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9979.477782123726);
    AssertAlmostEqual(&test, lnb, 9979.477773321123);

    casimir_lnab(casimir, 10000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9977.480249589197);
    AssertAlmostEqual(&test, lnb, 9977.479441673205);

    casimir_lnab(casimir, 10000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9929.451360635204);
    AssertAlmostEqual(&test, lnb, 9929.431512926210);

/*
    casimir_lnab(casimir, 10000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 9929.451360635204);
    AssertAlmostEqual(&test, lnb, 9929.431512926210);
*/

    casimir_lnab(casimir, 10000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 1.593263824351985);
    AssertAlmostEqual(&test, lnb, -5.51459483493011);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99974.89465752917);
    AssertAlmostEqual(&test, lnb, 99974.89465752757);

    casimir_lnab(casimir, 100000, 2, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99974.89457753077);
    AssertAlmostEqual(&test, lnb, 99974.89457752597);

    casimir_lnab(casimir, 100000, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99974.89445753317);
    AssertAlmostEqual(&test, lnb, 99974.89445752357);

    casimir_lnab(casimir, 100000, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99974.89409754038);
    AssertAlmostEqual(&test, lnb, 99974.89409751637);

    casimir_lnab(casimir, 100000, 10, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99974.89249757238);
    AssertAlmostEqual(&test, lnb, 99974.89249748438);

    casimir_lnab(casimir, 100000, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99974.69270163648);
    AssertAlmostEqual(&test, lnb, 99974.69269355627);

    casimir_lnab(casimir, 100000, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99969.88483954856);
    AssertAlmostEqual(&test, lnb, 99969.88463916263);

    casimir_lnab(casimir, 100000, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99954.87576562048);
    AssertAlmostEqual(&test, lnb, 99954.87496511695);

/*
    casimir_lnab(casimir, 100000, 5000, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);
    AssertAlmostEqual(&test, lna, 99954.87576562048);
    AssertAlmostEqual(&test, lnb, 99954.87496511695);
*/


    casimir_free(casimir);

    return test_results(&test, stderr);
}

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

    userdata[0] = 500; userdata[1] = 1;
    casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);

    casimir_lnab(casimir, 0.0001, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -30.1159277683113);
    AssertAlmostEqual(&test, lnb, -32.1420093429172);

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
    AssertAlmostEqual(&test, lna, -1930.61093220782);
    AssertAlmostEqual(&test, lnb, -1926.42631568006);

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
    AssertAlmostEqual(&test, lna, -3.86857027565537);
    AssertAlmostEqual(&test, lnb, 2.709147136660740);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99987.32386122631);
    AssertAlmostEqual(&test, lnb, 99987.32386122471);

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
    AssertAlmostEqual(&test, lna, -1144.29384239071);
    AssertAlmostEqual(&test, lnb, -1140.33952106964);

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
    AssertAlmostEqual(&test, lna, -1004.97005968866);
    AssertAlmostEqual(&test, lnb, -1001.43925678769);

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
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, -865.643023940761);
    AssertAlmostEqual(&test, lnb, -861.796912691024);

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
    AssertAlmostEqual(&test, lna, -2.06377318741045);
    AssertAlmostEqual(&test, lnb, 6.314608715823475);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99984.10499740121);
    AssertAlmostEqual(&test, lnb, 99984.10499739961);

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
    AssertAlmostEqual(&test, lna, -541.941551538154);
    AssertAlmostEqual(&test, lnb, -537.732803116565);

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
    AssertAlmostEqual(&test, lna, -1.32018117501920);
    AssertAlmostEqual(&test, lnb, 4.476291671480886);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99982.71870341508);
    AssertAlmostEqual(&test, lnb, 99982.71870341348);

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
    AssertAlmostEqual(&test, lna, -5.62453292046332);
    AssertAlmostEqual(&test, lnb, 6.351547202347875);

    casimir_lnab(casimir, 100000, 1, &lna, &lnb, &sign_a, &sign_b);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, 1);
    AssertAlmostEqual(&test, lna, 99974.89465752917);
    AssertAlmostEqual(&test, lnb, 99974.89465752757);

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

#include <math.h>

#include "libcasimir.h"
#include "integration_perf.h"
#include "unittest.h"
#include "sfunc.h"
#include "utils.h"

#include "test_integration_perf.h"

static void _integrals(int l1, int l2, int m, double nT, casimir_integrals_t *cint);
static double A(int l1, int l2, int m, double nT, sign_t *sign);
static double B(int l1, int l2, int m, double nT, sign_t *sign);
static double C(int l1, int l2, int m, double nT, sign_t *sign);
static float80 _I(int nu, int m2, double tau);

static float80 _I(int nu, int m2, double tau)
{
    float80 v;
    const double nT = tau/2.0;
    integration_perf_t int_perf;

    casimir_integrate_perf_init(&int_perf, nT, m2/2, nu);
    v = casimir_integrate_perf_I(&int_perf, nu);
    casimir_integrate_perf_free(&int_perf);

    return v;
}

/* return integral of A (including Lambda) for l1,l2,m small */
static double A(int l1, int l2, int m, double nT, sign_t *sign)
{
    const double tau = 2*nT;

    if(m == 0)
    {
        *sign = 1;
        return -INFINITY;
    }
    if(m == 1)
    {
        if(l1 == 1 && l2 == 1)
        {
            *sign = +1;
            return -tau+log(3)-log(4)-log(tau);
        }
        if(l1 == 1 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(sqrt(5)/sqrt(3)*3./4)-2*log(tau)+log1p(tau);
        }
        if(l1 == 1 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(sqrt(21)/8)-3*log(tau)+log(2*pow_2(tau)+5*tau+5);
        }

        if(l1 == 2 && l2 == 1)
        {
            *sign = +1;
            return -tau+log(sqrt(5)/sqrt(3)*3./4)-2*log(tau)+log1p(tau);
        }
        if(l1 == 2 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(5./4)-3*log(tau)+log(pow_2(tau)+2*tau+2);
        }
        if(l1 == 2 && l2 == 3)
        {
            *sign = +1;
            return -tau+0.5*log(11340./20736)-4*log(tau)+log(2*pow_3(tau)+7*pow_2(tau)+15*tau+15);
        }

        if(l1 == 3 && l2 == 1)
        {
            *sign = +1;
            return -tau+log(sqrt(21)*6./48)-3*log(tau)+log(2*pow_2(tau)+5*tau+5);
        }
        if(l1 == 3 && l2 == 2)
        {
            *sign = -1;
            return -tau+0.5*log(11340./20736)-4*log(tau)+log(2*pow_3(tau)+7*pow_2(tau)+15*tau+15);
        }
        if(l1 == 3 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(7./576)-5*log(tau)+log(144*pow_4(tau)+720*pow_3(tau)+2520*pow_2(tau)+5400*tau+5400);
        }
    }
    if(m == 2)
    {
        if(l1 == 2 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(5/2.)-3*log(tau)+log1p(tau);
        }
        if(l1 == 2 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(sqrt(35./12960)*90.)-4*log(tau)+log(pow_2(tau)+3*tau+3);
        }

        if(l1 == 3 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(sqrt(35./12960)*90.)-4*log(tau)+log(pow_2(tau)+3*tau+3);
        }
        if(l1 == 3 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(7./360)-5*log(tau)+log(450*pow_3(tau)+2250*pow_2(tau)+5400*tau+5400);
        }
    }

    TERMINATE(0, "not implemented for l1=%d, l2=%d, m=%d", l1,l2,m);
    return 0;
}

/* return integral of A (including Lambda) for l1,l2,m small */
static double B(int l1, int l2, int m, double nT, sign_t *sign)
{
    const double tau = 2*nT;

    if(m == 0)
    {
        if(l1 == 1 && l2 == 1)
        {
            *sign = -1;
            return -tau+log(3)-3*log(tau)+log1p(tau);
        }
        if(l1 == 1 && l2 == 2)
        {
            *sign = +1;
            return -tau+log(3*sqrt(5))-4*log(tau)+log(pow_2(tau)+3*tau+3);
        }

        if(l1 == 2 && l2 == 1)
        {
            *sign = -1;
            return -tau+log(3*sqrt(5))-4*log(tau)+log(pow_2(tau)+3*tau+3);
        }
        if(l1 == 2 && l2 == 2)
        {
            *sign = +1;
            return -tau+log(5./6)-5*log(tau)+log(18*pow_3(tau)+90*pow_2(tau)+216*tau+216);
        }
    }

    if(m == 1)
    {
        if(l1 == 1 && l2 == 1)
        {
            *sign = -1;
            return -tau+log(3./4)-3*log(tau)+log(pow_2(tau)+2*tau+2);
        }
        if(l1 == 1 && l2 == 2)
        {
            *sign = +1;
            return -tau+log(sqrt(5./3)*3./4)-4*log(tau)+log(pow_3(tau)+5*pow_2(tau)+12*tau+12);
        }
        if(l1 == 1 && l2 == 3)
        {
            *sign = -1;
            return -tau+log(6.*sqrt(21)/48)-5*log(tau)+log(2*pow_4(tau)+19*pow_3(tau)+79*pow_2(tau)+180*tau+180);
        }

        if(l1 == 2 && l2 == 1)
        {
            *sign = -1;
            return -tau+log(sqrt(5./3)*3./4)-4*log(tau)+log(pow_3(tau)+5*pow_2(tau)+12*tau+12);
        }
        if(l1 == 2 && l2 == 2)
        {
            *sign = +1;
            return -tau+log(45./36)-5*log(tau)+log(pow_4(tau)+8*pow_3(tau)+40*pow_2(tau)+96*tau+96);
        }
        if(l1 == 2 && l2 == 3)
        {
            *sign = -1;
            return -tau+log(sqrt(35./20736))-6*log(tau)+log(36*pow_5(tau)+450*pow_4(tau)+3402*pow_3(tau)+14202*pow_2(tau)+32400*tau+32400);
        }

        if(l1 == 3 && l2 == 1)
        {
            *sign = -1;
            return -tau+log(6.*sqrt(21)/48)-5*log(tau)+log(2*pow_4(tau)+19*pow_3(tau)+79*pow_2(tau)+180*tau+180);
        }
        if(l1 == 3 && l2 == 2)
        {
            *sign = +1;
            return -tau+log(sqrt(35./20736))-6*log(tau)+log(36*pow_5(tau)+450*pow_4(tau)+3402*pow_3(tau)+14202*pow_2(tau)+32400*tau+32400);
        }
        if(l1 == 3 && l2 == 3)
        {
            *sign = -1;
            return -tau+log(7/576.)-7*log(tau)+log(144*pow_6(tau)+2448*pow_5(tau)+27288*pow_4(tau)+171720*pow_3(tau)+657720*pow_2(tau)+1458000*tau+1458000);
        }
    }
    if(m == 2)
    {
        if(l1 == 2 && l2 == 2)
        {
            *sign = +1;
            return -tau+log(5./144)-5*log(tau)+log(72*pow_3(tau)+360*pow_2(tau)+864*tau+864);
        }
        if(l1 == 2 && l2 == 3)
        {
            *sign = -1;
            return -tau+0.5*log(35./207360)-6*log(tau)+log(360*pow_4(tau)+3240*pow_3(tau)+14040*pow_2(tau)+32400*tau+32400);
        }

        if(l1 == 3 && l2 == 2)
        {
            *sign = +1;
            return -tau+0.5*log(35./207360)-6*log(tau)+log(360*pow_4(tau)+3240*pow_3(tau)+14040*pow_2(tau)+32400*tau+32400);
        }
        if(l1 == 3 && l2 == 3)
        {
            *sign = -1;
            return -tau+log(7./1440)-7*log(tau)+log(1800*pow_5(tau)+23400*pow_4(tau)+162000*pow_3(tau)+648000*pow_2(tau)+1458000*tau+1458000);
        }
    }

    TERMINATE(1, "not implemented for l1=%d, l2=%d, m=%d", l1,l2,m);
    return 0;
}

/* return integral of C (including Lambda) for l1,l2,m small */
static double C(int l1, int l2, int m, double nT, sign_t *sign)
{
    const double tau = 2*nT;

    if(m == 0)
    {
        *sign = 1;
        return -INFINITY;
    }
    if(m == 1)
    {
        if(l1 == 1 && l2 == 1)
        {
            *sign = +1;
            return -tau+log(3./4)-2*log(tau)+log1p(tau);
        }
        if(l1 == 1 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(sqrt(5./3)*3./4)-3*log(tau)+log(pow_2(tau)+4*tau+4);
        }
        if(l1 == 1 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(sqrt(21)/48)-4*log(tau)+log(12*pow_3(tau)+102*pow_2(tau)+270*tau+270);
        }

        if(l1 == 2 && l2 == 1)
        {
            *sign = +1;
            return -tau+log(sqrt(5./3)*3./4)-3*log(tau)+log(pow_2(tau)+2*tau+2);
        }
        if(l1 == 2 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(5./36)-4*log(tau)+log(9*pow_3(tau)+45*pow_2(tau)+108*tau+108);
        }
        if(l1 == 2 && l2 == 3)
        {
            *sign = +1;
            return -tau+0.5*log(35./20736)-5*log(tau)+log(36*pow_4(tau)+342*pow_3(tau)+1422*pow_2(tau)+3240*tau+3240);
        }

        if(l1 == 3 && l2 == 1)
        {
            *sign = +1;
            return -tau+log(sqrt(21)/48)-4*log(tau)+log(12*pow_3(tau)+42*pow_2(tau)+90*tau+90);
        }
        if(l1 == 3 && l2 == 2)
        {
            *sign = -1;
            return -tau+0.5*log(35./20736)-5*log(tau)+log(36*pow_4(tau)+234*pow_3(tau)+954*pow_2(tau)+2160*tau+2160);
        }
        if(l1 == 3 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(7./576)-6*log(tau)+log(144*pow_5(tau)+1584*pow_4(tau)+9720*pow_3(tau)+36720*pow_2(tau)+81000*tau+81000);
        }
    }
    if(m == 2)
    {
        if(l1 == 2 && l2 == 2)
        {
            *sign = -1;
            return -tau+log(5./72)-4*log(tau)+log(36*pow_2(tau)+108*tau+108);
        }
        if(l1 == 2 && l2 == 3)
        {
            *sign = +1;
            return -tau+0.5*log(35./51840)-5*log(tau)+log(180*pow_3(tau)+1260*pow_2(tau)+3240*tau+3240);
        }

        if(l1 == 3 && l2 == 2)
        {
            *sign = -1;
            return -tau+0.5*log(35./51840)-5*log(tau)+log(180*pow_3(tau)+900*pow_2(tau)+2160*tau+2160);
        }
        if(l1 == 3 && l2 == 3)
        {
            *sign = +1;
            return -tau+log(7./720)-6*log(tau)+log(900*pow_4(tau)+8100*pow_3(tau)+35100*pow_2(tau)+81000*tau+81000);
        }
    }

    TERMINATE(0, "not implemented for l1=%d, l2=%d, m=%d", l1,l2,m);
    return 0;
}

static void _integrals(int l1, int l2, int m, double nT, casimir_integrals_t *cint)
{
    integration_perf_t int_perf;
    int lmax = MAX(l1,l2);
    casimir_integrate_perf_init(&int_perf, nT, m, lmax);
    casimir_integrate_perf(&int_perf, l1, l2, cint);
    casimir_integrate_perf_free(&int_perf);
}

int test_integration_perf(void)
{
    double v;
    double list_nT[] = { 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1e0, 1.006, 2, 1e1, 1e2, 1e3, 1e4 };
    casimir_integrals_t cint;
    unittest_t test;
    unittest_init(&test, "Integration", "Test integration for various parameters");

    /* test vs analytical expressions for l1,l2,m small */
    for(size_t i = 0; i < sizeof(list_nT)/sizeof(list_nT[0]); i++)
    {
        double nT = list_nT[i];

        for(int l1 = 1; l1 <= 2; l1++)
            for(int l2 = 1; l2 <= 2; l2++)
            {
                sign_t signB = 0;
                double lnB = B(l1,l2,0,nT,&signB);

                _integrals(l1,l2,0,nT,&cint);

                AssertAlmostEqual(&test, cint.lnB_TM, lnB);
                AssertAlmostEqual(&test, cint.signB_TM, signB);
            }

        for(int m = 1; m <= 2; m++)
            for(int l1 = 2; l1 <= 3; l1++)
                for(int l2 = 2; l2 <= 3; l2++)
                {
                    sign_t signA = 0, signB = 0, signC = 0, signD = 0;
                    double lnA = A(l1,l2,m,nT,&signA);
                    double lnB = B(l1,l2,m,nT,&signB);
                    double lnC = C(l1,l2,m,nT,&signC);
                    double lnD = C(l2,l1,m,nT,&signD);

                    _integrals(l1,l2,m,nT,&cint);
                    AssertAlmostEqual(&test, cint.lnA_TM, lnA);
                    AssertAlmostEqual(&test, cint.signA_TM, signA);

                    AssertAlmostEqual(&test, cint.lnB_TM, lnB);
                    AssertAlmostEqual(&test, cint.signB_TM, signB);

                    AssertAlmostEqual(&test, cint.lnC_TM, lnC);
                    AssertAlmostEqual(&test, cint.signC_TM, signC);

                    AssertAlmostEqual(&test, cint.signD_TM, -signC);
                    AssertAlmostEqual(&test, cint.lnD_TM, lnD);
                    AssertAlmostEqual(&test, cint.signD_TM, MPOW(l1+l2+1)*signD);
                }
    }

    /* test integral I */
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

    v = _I(4000,2,2);
    AssertAlmostEqual(&test, v, 29179.23802073103788347465983473258110284);

    /* test for various combinations of l1,l2,m and nT */
    _integrals(1050,1050,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 13940.125756903571190096123124106829015888756201576644);
    AssertAlmostEqual(&test, cint.lnB_TM, 13967.9514623712062733110267893316781893739350523);
    AssertAlmostEqual(&test, cint.lnC_TM, 13954.03837148533464206);

    _integrals(1050,1,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 6244.485496882807753109);
    AssertAlmostEqual(&test, cint.lnB_TM, 6263.969792590335805);
    AssertAlmostEqual(&test, cint.lnC_TM, 6250.748896960319286622);

    _integrals(1,1050,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 6244.485496882807753109081179);
    AssertAlmostEqual(&test, cint.lnB_TM, 6263.969792590335805037405705);
    AssertAlmostEqual(&test, cint.lnC_TM, 6257.7054405868224491877);

    _integrals(500,1050,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 9813.28363230087848);
    AssertAlmostEqual(&test, cint.lnB_TM, 9839.7598665288218618873847);
    AssertAlmostEqual(&test, cint.lnC_TM, 9826.8923954024132364);

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
    AssertAlmostEqual(&test, cint.lnA_TM, 40.74560144887208);
    AssertEqual(&test, cint.signA_TM, -1);
    AssertAlmostEqual(&test, cint.lnB_TM, 56.75388164708835);
    AssertEqual(&test, cint.signB_TM, +1);
    AssertAlmostEqual(&test, cint.lnC_TM, 48.68297568585137);
    AssertEqual(&test, cint.signC_TM, -1);
    AssertAlmostEqual(&test, cint.lnD_TM, 48.68297568585137);
    AssertEqual(&test, cint.signD_TM, +1);

    _integrals(40,40,1,0.25,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 367.258716490769);
    AssertAlmostEqual(&test, cint.lnB_TM, 384.7742430107486);
    AssertAlmostEqual(&test, cint.lnC_TM, 376.010190217081);
    AssertAlmostEqual(&test, cint.lnD_TM, 376.010190217081);

    _integrals(40,40,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 257.7297499756845);
    AssertAlmostEqual(&test, cint.lnB_TM, 272.472669228606);
    AssertAlmostEqual(&test, cint.lnC_TM, 265.0949179301248);
    AssertAlmostEqual(&test, cint.lnD_TM, 265.0949179301248);

    _integrals(40,40,40,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 212.0868844486187);
    AssertAlmostEqual(&test, cint.lnB_TM, 219.4533537701274);
    AssertAlmostEqual(&test, cint.lnC_TM, 215.7638420201849);
    AssertAlmostEqual(&test, cint.lnD_TM, 215.7638420201849);

    _integrals(7,4,3,8,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -17.1928436859713);
    AssertAlmostEqual(&test, cint.lnB_TM, -16.0865392641165);
    AssertAlmostEqual(&test, cint.lnC_TM, -16.83090135860425);
    AssertAlmostEqual(&test, cint.lnD_TM, -16.48574215820564);

    _integrals(20, 20, 3, 13, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -5.202385074993125);
    AssertAlmostEqual(&test, cint.lnB_TM, -0.5773666089005467);
    AssertAlmostEqual(&test, cint.lnC_TM, -2.905312825257782);
    AssertAlmostEqual(&test, cint.lnD_TM, -2.905312825257782);

    _integrals(20, 20, 3, 0.001, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 368.3408666279195);
    AssertAlmostEqual(&test, cint.lnB_TM, 391.916763894729);
    AssertAlmostEqual(&test, cint.lnC_TM, 380.1161563573135);
    AssertAlmostEqual(&test, cint.lnD_TM, 380.1161563573135);

    _integrals(20, 20, 3, 0.01, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 278.540045473523);
    AssertAlmostEqual(&test, cint.lnB_TM, 297.5107725493697);
    AssertAlmostEqual(&test, cint.lnC_TM, 288.0127501055992);
    AssertAlmostEqual(&test, cint.lnD_TM, 288.0127501055992);

    _integrals(20, 20, 3, 0.1, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 188.7389740846992);
    AssertAlmostEqual(&test, cint.lnB_TM, 203.1045304770997);
    AssertAlmostEqual(&test, cint.lnC_TM, 195.9090931914071);
    AssertAlmostEqual(&test, cint.lnD_TM, 195.9090931914071);

    _integrals(20, 20, 3, 1, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 98.91288744548656);
    AssertAlmostEqual(&test, cint.lnB_TM, 108.6732240307656);
    AssertAlmostEqual(&test, cint.lnC_TM, 103.7803782982835);
    AssertAlmostEqual(&test, cint.lnD_TM, 103.7803782982835);

    _integrals(20, 20, 3, 10, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 6.660701378819416);
    AssertAlmostEqual(&test, cint.lnB_TM, 11.81202300232528);
    AssertAlmostEqual(&test, cint.lnC_TM, 9.221979692552173);
    AssertAlmostEqual(&test, cint.lnD_TM, 9.221979692552173);

    _integrals(20, 20, 3, 100, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -201.8065316849248);
    AssertAlmostEqual(&test, cint.lnB_TM, -200.7135958357346);
    AssertAlmostEqual(&test, cint.lnC_TM, -201.2774897833974);
    AssertAlmostEqual(&test, cint.lnD_TM, -201.2774897833974);

    _integrals(20, 10, 3, 0.1, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 131.0962726931826);
    AssertAlmostEqual(&test, cint.lnB_TM, 144.184735156638);
    AssertAlmostEqual(&test, cint.lnC_TM, 137.2769797021334);
    AssertAlmostEqual(&test, cint.lnD_TM, 137.9701257824486);

    _integrals(20, 15, 10, 5, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 22.41543637237407);
    AssertAlmostEqual(&test, cint.lnB_TM, 26.04398834917292);
    AssertAlmostEqual(&test, cint.lnC_TM, 24.07623857473362);
    AssertAlmostEqual(&test, cint.lnD_TM, 24.35571686776128);

    _integrals(50, 15, 10, 10, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 45.60713807155988);
    AssertAlmostEqual(&test, cint.lnB_TM, 49.99651015278684);
    AssertAlmostEqual(&test, cint.lnC_TM, 47.201688624453);
    AssertAlmostEqual(&test, cint.lnD_TM, 48.38666367876067);

    _integrals(100, 25, 20, 20, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 84.50837701530297);
    AssertAlmostEqual(&test, cint.lnB_TM, 88.66027006552913);
    AssertAlmostEqual(&test, cint.lnC_TM, 85.90224599747752);
    AssertAlmostEqual(&test, cint.lnD_TM, 87.2586869597479);

    _integrals(60, 55, 40, 0.11, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 645.1730223683922);
    AssertAlmostEqual(&test, cint.lnB_TM, 658.4063308419369);
    AssertAlmostEqual(&test, cint.lnC_TM, 651.7418041758159);
    AssertAlmostEqual(&test, cint.lnD_TM, 651.8288153934492);

    _integrals(40,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 185.27722707813169721211989855051);
    AssertAlmostEqual(&test, cint.lnB_TM, 198.1874611137788);
    AssertAlmostEqual(&test, cint.lnC_TM, 191.72604045861798912);

    _integrals(140,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 575.400220880156994701641252076629);
    AssertAlmostEqual(&test, cint.lnB_TM, 591.1921970497888542);
    AssertAlmostEqual(&test, cint.lnC_TM, 582.6670391327286);

    _integrals(240,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 1025.59108523802829059595981668750);
    AssertAlmostEqual(&test, cint.lnB_TM, 1042.807725206889176);
    AssertAlmostEqual(&test, cint.lnC_TM, 1033.30173547300948);

    _integrals(540,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 2561.62734676892999652846813001817);
    AssertAlmostEqual(&test, cint.lnB_TM, 2581.1132494456470771627);
    AssertAlmostEqual(&test, cint.lnC_TM, 2570.068090210887226712719);

    _integrals(1540,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 8589.22040500894307686493465726593);
    AssertAlmostEqual(&test, cint.lnB_TM, 8611.759673402576546371338);
    AssertAlmostEqual(&test, cint.lnC_TM, 8598.66439349824405764);

    _integrals(2000,40,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 11616.33666207935);
    AssertAlmostEqual(&test, cint.lnB_TM, 11639.648487984291931);
    AssertAlmostEqual(&test, cint.lnC_TM, 11626.03631835250693323);

    _integrals(2000,1,1,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 11358.3678862878775413115);
    AssertAlmostEqual(&test, cint.lnB_TM, 11377.9522208393254813978);
    AssertAlmostEqual(&test, cint.lnC_TM, 11364.359353960757194524);

    _integrals(200,1,0,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnB_TM, 682.0020149230455);

    _integrals(2000,1,0,2.5,&cint);
    AssertAlmostEqual(&test, cint.lnB_TM, 11378.29904124447803687331);

    _integrals(2000,40,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 13484.656037918688437);
    AssertAlmostEqual(&test, cint.lnB_TM, 13509.800445322035821193049198451);
    AssertAlmostEqual(&test, cint.lnC_TM, 13495.271984956561221915441695);


    _integrals(2000,2000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 29149.71533012066508345912);
    AssertAlmostEqual(&test, cint.lnB_TM, 29180.1186899274219145148);
    AssertAlmostEqual(&test, cint.lnC_TM, 29164.91688500840021851842);

    _integrals(2000,1000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 20993.7437127492067941699865154);
    AssertAlmostEqual(&test, cint.lnB_TM, 21022.8784778726214390197843);
    AssertAlmostEqual(&test, cint.lnC_TM, 21007.9643550261185783165);

    _integrals(2000,1,1,0.5,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 14577.246711903825292880294853452);
    AssertAlmostEqual(&test, cint.lnB_TM, 14600.0499192823994956727755285068);
    AssertAlmostEqual(&test, cint.lnC_TM, 14584.84761448839861741754014567215621);

    _integrals(2000,1,1,10,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 8585.7322779497324918261);
    AssertAlmostEqual(&test, cint.lnB_TM, 8602.54407061632099700814913677598);
    AssertAlmostEqual(&test, cint.lnC_TM, 8590.3374981457216867851698439);

    _integrals(4000,1,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 29164.30634417171784698325);
    AssertAlmostEqual(&test, cint.lnB_TM, 29187.80244882461237368);
    AssertAlmostEqual(&test, cint.lnC_TM, 29171.907246756275540666256);

    _integrals(4000,2000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 46169.30412583713741863288);
    AssertAlmostEqual(&test, cint.lnB_TM, 46201.211646391476308900991);
    AssertAlmostEqual(&test, cint.lnC_TM, 46184.911229183740231533516827);

    _integrals(4000,4000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 63868.6657467722265411811);
    AssertAlmostEqual(&test, cint.lnB_TM, 63901.841820324801967);
    AssertAlmostEqual(&test, cint.lnC_TM, 63885.2537210446057222743554823);

    _integrals(500, 500, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 5880.145418565579706653290311279924286563566945608969682962480558534581786);
    AssertAlmostEqual(&test, cint.lnB_TM, 5904.990886305616635574177976966139243575449897572646531806921581807454247);
    AssertAlmostEqual(&test, cint.lnC_TM, 5892.567652184402156243850461817790605864901360854261291722500861621079491);
    AssertAlmostEqual(&test, cint.lnD_TM, 5892.567652184402156243850461817790605864901360854261291722500861621079491);

    _integrals(499, 500, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 5873.247644850501700848005625848572206294394325680038496156282856954543478);
    AssertAlmostEqual(&test, cint.lnB_TM, 5898.089108585147827540904685136290828084996798178646217185451773987029079);
    AssertAlmostEqual(&test, cint.lnC_TM, 5885.668876966959763457112814499634871743366485689495287378961086374547958);
    AssertAlmostEqual(&test, cint.lnD_TM, 5885.666874964285038289757983818757214367337478685908564744676071164769226);

    _integrals(500, 499, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 5873.247644850501700848005625848572206294394325680038496156282856954543478);
    AssertAlmostEqual(&test, cint.lnB_TM, 5898.089108585147827540904685136290828084996798178646217185451773987029079);
    AssertAlmostEqual(&test, cint.lnC_TM, 5885.666874964285038289757983818757214367337478685908564744676071164769226);
    AssertAlmostEqual(&test, cint.lnD_TM, 5885.668876966959763457112814499634871743366485689495287378961086374547958);

    _integrals(499, 499, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 5866.350873639756594278938154316878502875775073387778089826355287796266225);
    AssertAlmostEqual(&test, cint.lnB_TM, 5891.188331362990374000578938922984773493968146030458032951120864889446638);
    AssertAlmostEqual(&test, cint.lnC_TM, 5878.769101247160551296372865772325533805896297034676801385673068702228523);
    AssertAlmostEqual(&test, cint.lnD_TM, 5878.769101247160551296372865772325533805896297034676801385673068702228523);

    _integrals(200, 200, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 1975.767096926493948364578286189828228822723485260602372040682771555386593);
    AssertAlmostEqual(&test, cint.lnB_TM, 1996.945898960585135363127239781841983882569826371795585966870453220538568);
    AssertAlmostEqual(&test, cint.lnC_TM, 1986.35524636210021458561037266323954886067698626041262792808131074279971);
    AssertAlmostEqual(&test, cint.lnD_TM, 1986.35524636210021458561037266323954886067698626041262792808131074279971);

    _integrals(201, 199, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 1975.76712170967901909630632990248917658766444079226758033464543591258041);
    AssertAlmostEqual(&test, cint.lnB_TM, 1996.945898743456163521589687237776467526903757736254664642401390499554195);
    AssertAlmostEqual(&test, cint.lnC_TM, 1986.350258603304040900412142964201303792493786158794697788127912783088005);
    AssertAlmostEqual(&test, cint.lnD_TM, 1986.360258686952479420258637033592832162969783950835808872636847179108004);

    _integrals(1, 2, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -2.339923455106124005554798782243353095693901455689575867856081014682324586);
    AssertAlmostEqual(&test, cint.lnB_TM, -0.6736624574273509838711334310684649842252665724200226758604799345224393973);
    AssertAlmostEqual(&test, cint.lnC_TM, -1.363077277081885731378935504250430542739743833222871115852865737072863219);
    AssertAlmostEqual(&test, cint.lnD_TM, -1.831883423580260237202605565635544138614431556772498548761012276038582097);

    _integrals(2, 1, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -2.339923455106124005554798782243353095693901455689575867856081014682324586);
    AssertAlmostEqual(&test, cint.lnB_TM, -0.6736624574273509838711334310684649842252665724200226758604799345224393973);
    AssertAlmostEqual(&test, cint.lnC_TM, -1.831883423580260237202605565635544138614431556772498548761012276038582097);
    AssertAlmostEqual(&test, cint.lnD_TM, -1.363077277081885731378935504250430542739743833222871115852865737072863219);

    _integrals(1, 1, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, -2.998811324689273716593897834310454621521738584722626358695413310706611303);
    AssertAlmostEqual(&test, cint.lnB_TM, -2.08729623546325568208182858341472166512674251862231861269070482888007109);
    AssertAlmostEqual(&test, cint.lnC_TM, -2.595336266989119347157555830395184063132956853912460002945057793600666948);
    AssertAlmostEqual(&test, cint.lnD_TM, -2.595336266989119347157555830395184063132956853912460002945057793600666948);

    _integrals(1, 500, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 2595.536799689232189132097458491885340008891921765359782557085782932445143);
    AssertAlmostEqual(&test, cint.lnB_TM, 2612.784371554692815301128903227217788398596035743007173804862451217505851);
    AssertAlmostEqual(&test, cint.lnC_TM, 2607.26688661763036080100424831203262202790938673585476359897048874019501);
    AssertAlmostEqual(&test, cint.lnD_TM, 2601.052286639743501207037424084674992544735252498822692375343512556374755);

    _integrals(500, 1, 1, 1.006, &cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 2595.536799689232189132097458491885340008891921765359782557085782932445143);
    AssertAlmostEqual(&test, cint.lnB_TM, 2612.784371554692815301128903227217788398596035743007173804862451217505851);
    AssertAlmostEqual(&test, cint.lnC_TM, 2601.052286639743501207037424084674992544735252498822692375343512556374755);
    AssertAlmostEqual(&test, cint.lnD_TM, 2607.26688661763036080100424831203262202790938673585476359897048874019501);

    _integrals(4000,4000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 63868.6657467722265411811000471834427);
    AssertAlmostEqual(&test, cint.lnB_TM, 63901.84182032480196706998988829691231125220392684855064349853042);

    _integrals(4000,2000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 46169.30412583713741863288035766593016712961550218654);
    AssertAlmostEqual(&test, cint.lnB_TM, 46201.2116463914763089009910622677286368913618141444483663);
    AssertAlmostEqual(&test, cint.lnC_TM, 46184.911229183740231533516827);

    _integrals(4000,1000,1,1,&cint);
    AssertAlmostEqual(&test, cint.lnA_TM, 37559.147788784669482290944857175749);
    AssertAlmostEqual(&test, cint.lnC_TM, 37573.8793900544332570436073011074230654516);

    return test_results(&test, stderr);
}

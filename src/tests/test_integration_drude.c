#include <math.h>

#include "libcasimir.h"
#include "integration_drude.h"
#include "unittest.h"
#include "sfunc.h"

#include "test_integration_drude.h"

#define Assert_Drude(test, x, y) _AssertAlmostEqual(__LINE__, (test), (x), (y), DRUDE_INTEG_ACCURACY)

static void drude_integrate(casimir_t* casimir, casimir_integrals_t* cint,
                            int l1, int l2, int m, int n, double T)
{
    integration_drude_t int_drude;
    casimir_integrate_drude_init(casimir, &int_drude, n * T, m, MAX(l1, l2));
    casimir_integrate_drude(&int_drude, l1, l2, cint);
    casimir_integrate_drude_free(&int_drude);
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
    casimir_init(&casimir, 1/0.5-1, 1);
    casimir_set_omegap_plane(&casimir, omegap);
    casimir_set_gamma_plane(&casimir, gamma_);

    {
        drude_integrate(&casimir, &cint, 250, 250, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 2584.8354863970);
        AssertEqual(&test, cint.signA_TE, +1);
        Assert_Drude(&test, cint.lnA_TM, 2588.0926913952);
        AssertEqual(&test, cint.signA_TM, -1);

    }

    {
        drude_integrate(&casimir, &cint, 3, 2, 1, 1, 1);

        Assert_Drude(&test, cint.lnA_TE, -0.62981145199252068602408);
        AssertEqual(&test, cint.signA_TE, +1);
        Assert_Drude(&test, cint.lnA_TM, -0.59589434712666196879313);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 2.7383266248198112347328);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 2.7945369735963442460566);
        AssertEqual(&test, cint.signB_TM, +1);

        Assert_Drude(&test, cint.lnC_TE, 0.74158885587407677484842);
        AssertEqual(&test, cint.signC_TE, +1);
        Assert_Drude(&test, cint.lnC_TM, 0.78654850992186874297884);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1.1256865694785886272321);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1.1711803659444285938208);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        drude_integrate(&casimir, &cint, 3, 2, 2, 1, 1);

        Assert_Drude(&test, cint.lnA_TE, -0.7132786835392315505014);
        AssertEqual(&test, cint.signA_TE, +1);
        Assert_Drude(&test, cint.lnA_TM, -0.67405454481613061169669);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1.5622997116874727152691);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1.6191154482067796685377);
        AssertEqual(&test, cint.signB_TM, +1);

        Assert_Drude(&test, cint.lnC_TE, 0.18064782184885164701636);
        AssertEqual(&test, cint.signC_TE, +1);
        Assert_Drude(&test, cint.lnC_TM, 0.22781988160375950931957);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 0.52072721231985447793985);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 0.56889321108023160164297);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        drude_integrate(&casimir, &cint, 3, 2, 1, 1, 2);

        Assert_Drude(&test, cint.lnA_TE, -4.1459191747317624709052);
        AssertEqual(&test, cint.signA_TE, +1);
        Assert_Drude(&test, cint.lnA_TM, -4.1197945729671869295841);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, -2.0359815567492122711267);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, -1.9903481981146134279244);
        AssertEqual(&test, cint.signB_TM, +1);

        Assert_Drude(&test, cint.lnC_TE, -3.356496453224571521287);
        AssertEqual(&test, cint.signC_TE, +1);
        Assert_Drude(&test, cint.lnC_TM, -3.3216912696961575288389);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, -3.0252864555803092122331);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, -2.9891216277980328174213);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        drude_integrate(&casimir, &cint, 4, 2, 1, 1, 2);

        Assert_Drude(&test, cint.lnA_TE, -3.4410543260111500276103);
        AssertEqual(&test, cint.signA_TE, +1);
        Assert_Drude(&test, cint.lnA_TM, -3.4078967115562010093212);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, -0.74990261955450098567219);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, -0.69545756637556400940693);
        AssertEqual(&test, cint.signB_TM, +1);

        Assert_Drude(&test, cint.lnC_TE, -2.5076062999924391950641);
        AssertEqual(&test, cint.signC_TE, +1);
        Assert_Drude(&test, cint.lnC_TM, -2.4646103122548521433085);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, -1.883875993266507818876);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, -1.8392724297362885494622);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        drude_integrate(&casimir, &cint, 4, 3, 1, 1, 2);

        Assert_Drude(&test, cint.lnA_TE, -2.6626325469377015011493);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, -2.6217256954918375062335);
        AssertEqual(&test, cint.signA_TM, +1);

        Assert_Drude(&test, cint.lnB_TE, 0.69662250934949833938181);
        AssertEqual(&test, cint.signB_TE, +1);
        Assert_Drude(&test, cint.lnB_TM, 0.76016712598231202383933);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, -1.2220766351698472482591);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, -1.1696039758614060731918);
        AssertEqual(&test, cint.signC_TM, +1);

        Assert_Drude(&test, cint.lnD_TE, -0.94454349314724581653695);
        AssertEqual(&test, cint.signD_TE, +1);
        Assert_Drude(&test, cint.lnD_TM, -0.89174332617421184874047);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 4, 3, 1, 1, 20);

        Assert_Drude(&test, cint.lnA_TE, -43.137947441949791493356);
        Assert_Drude(&test, cint.lnA_TM, -43.120372668298594245786);
        Assert_Drude(&test, cint.lnB_TE, -42.712836955286222313897);
        Assert_Drude(&test, cint.lnB_TM, -42.688853217438179634812);
        Assert_Drude(&test, cint.lnC_TE, -42.981679865536465578416);
        Assert_Drude(&test, cint.lnC_TM, -42.961747317757136470438);
        Assert_Drude(&test, cint.lnD_TE, -42.896313862496345284576);
        Assert_Drude(&test, cint.lnD_TM, -42.875331151064598797301);
    }

    {
        drude_integrate(&casimir, &cint, 4, 3, 1, 1, 0.01);

        Assert_Drude(&test, cint.lnA_TE, 29.361786303876121509058);
        Assert_Drude(&test, cint.lnA_TM, 29.727644399301097042065);
        Assert_Drude(&test, cint.lnB_TE, 43.288690290377248702015);
        Assert_Drude(&test, cint.lnB_TM, 43.774264506257246488789);
        Assert_Drude(&test, cint.lnC_TE, 36.104129446007204258962);
        Assert_Drude(&test, cint.lnC_TM, 36.53003707484319590354);
        Assert_Drude(&test, cint.lnD_TE, 36.39181148264103954446);
        Assert_Drude(&test, cint.lnD_TM, 36.817719115542167950645);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 0, 1, 2);

        Assert_Drude(&test, cint.lnB_TE, 80.261829579383622339087);
        Assert_Drude(&test, cint.lnB_TM, 80.616659994373914035408);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 1, 1, 2);

        Assert_Drude(&test, cint.lnA_TE, 69.659355648443184396324);
        Assert_Drude(&test, cint.lnA_TM, 69.996776884623239716192);
        Assert_Drude(&test, cint.lnB_TE, 80.212799262982493059119);
        Assert_Drude(&test, cint.lnB_TM, 80.567633810049211725342);
        Assert_Drude(&test, cint.lnC_TE, 74.92332975118821662198);
        Assert_Drude(&test, cint.lnC_TM, 75.269464691472079467208);
        Assert_Drude(&test, cint.lnD_TE, 74.92332975118821662198);
        Assert_Drude(&test, cint.lnD_TM, 75.269464691472079467208);
    }

    {
        drude_integrate(&casimir, &cint, 7, 7, 1, 1, 2);

        Assert_Drude(&test, cint.lnA_TE, 6.6971469912051709882689);
        Assert_Drude(&test, cint.lnA_TM, 6.8044855533996484279007);
        Assert_Drude(&test, cint.lnB_TE, 12.986970305176281775914);
        Assert_Drude(&test, cint.lnB_TM, 13.113563515622480422875);
        Assert_Drude(&test, cint.lnC_TE, 9.8017901880975337522829);
        Assert_Drude(&test, cint.lnC_TM, 9.9188754798338178698508);
        Assert_Drude(&test, cint.lnD_TE, 9.8017901880975337522829);
        Assert_Drude(&test, cint.lnD_TM, 9.9188754798338178698508);
    }

    {
        drude_integrate(&casimir, &cint, 60, 7, 1, 1, 4);

        Assert_Drude(&test, cint.lnA_TE, 111.142136572991446682);
        Assert_Drude(&test, cint.lnA_TM, 111.67030961586370531194);
        Assert_Drude(&test, cint.lnB_TE, 121.40133118469416742774);
        Assert_Drude(&test, cint.lnB_TM, 121.94548929011007049376);
        Assert_Drude(&test, cint.lnC_TE, 115.18974866771901373092);
        Assert_Drude(&test, cint.lnC_TM, 115.72592401117753046786);
        Assert_Drude(&test, cint.lnD_TE, 117.33853608403136877599);
        Assert_Drude(&test, cint.lnD_TM, 117.87470585363882880168);
    }

    {
        drude_integrate(&casimir, &cint, 60, 50, 1, 1, 2);

        Assert_Drude(&test, cint.lnA_TE, 316.6241201381509729508);
        Assert_Drude(&test, cint.lnA_TM, 317.54352315209682942627);
        Assert_Drude(&test, cint.lnB_TE, 331.23393499264542810725);
        Assert_Drude(&test, cint.lnB_TM, 332.16912114975629911808);
        Assert_Drude(&test, cint.lnC_TE, 323.83328989832645684183);
        Assert_Drude(&test, cint.lnC_TM, 324.76059133413377133818);
        Assert_Drude(&test, cint.lnD_TE, 324.01561254123456600265);
        Assert_Drude(&test, cint.lnD_TM, 324.94291396001144754842);
    }

    {
        drude_integrate(&casimir, &cint, 60, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 157.03473913126033297558);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 157.60288524555978642832);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 168.67977888414012100238);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 169.2647676475485558551);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 161.775456411727348781);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 162.35203030735193408451);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 163.92397950412386236859);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 164.50055191892345523067);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 80, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 229.1597844213887940905);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 229.89406457143829405924);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 241.61894305893617786564);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 242.36960285132202690367);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 234.16548444947210610469);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 234.90796111111401757985);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 236.60165519648842639912);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 237.34413097602283063093);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 100, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 306.06932889220031152131);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 306.96497345975102872487);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 319.16818541918637605961);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 320.07969458823978796966);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 311.28440340692612802115);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 312.18798734237934246671);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 313.94369995456691501291);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 314.84728331471311765311);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 120, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 386.83594323128082725092);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 387.8876642670721984086);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 400.46188992008071776);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 401.52892119280114443169);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 392.22415068164628194729);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 393.2835341808375725671);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 395.06575848369865149859);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 396.12514158411000949863);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 140, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 470.83515226908617063882);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 472.03731825287765654347);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 484.9094049849202701743);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 486.12630342758881220105);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 476.37098088090389735611);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 477.58052056713776968274);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 479.36673286632467557876);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 480.57627226379758413359);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 160, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 557.61732254018995416991);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 558.9641240774240227426);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 572.08164013464683280655);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 573.44258753277950544768);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 563.2818304865038356111);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 564.63571241269611007078);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 566.41110950973999424267);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 567.76499121982570015936);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 180, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 646.8430237940884399506);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 648.32861298950403491307);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 661.65258453933260639722);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 663.15173611443258024946);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 652.62158575898876955012);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 654.11396347151462447173);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 655.86864475862975153598);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 657.36102230517159947416);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 200, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 738.2468187234904583539);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 739.86541803510869960269);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 753.36605299727538560123);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 754.99764338175436799439);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 744.12779835366667173947);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 745.75290031062115515565);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 747.4802156387367917563);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 749.10531746548697028254);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 220, 7, 1, 1, 2.0);

        Assert_Drude(&test, cint.lnA_TE, 831.6153798455399456949);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 833.36136190640036566345);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 847.0153685285462429851);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 848.77378895540423099295);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 837.5892966203457351438);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 839.3415046949447144509);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 841.03702240982617910868);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 842.78923038047398507821);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 80, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 175.36573435960625587165);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 176.0661511230820160639);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 186.73638175840007575614);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 187.45227275413481610286);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 182.13786030088643692502);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 182.84601910010091710651);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 179.95289050868875764273);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 180.66105187994419171093);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 100, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 238.86268895832442365318);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 239.71519039592288738893);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 250.86454975778596677592);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 251.7320528124298849072);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 246.06302060626087693709);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 246.92302826176329044283);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 243.65496043495794765477);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 244.51496979577075616412);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 120, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 306.13652842451938959492);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 307.13647720621826995181);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 318.65961267082135116666);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 319.67406951388311410131);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 313.6893396860856971088);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 314.69654831785864390224);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 311.09898946274578258007);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 312.1061992912919397839);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 160, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 449.7790071372679505321);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 451.05896719346032146996);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 463.13286037231568336792);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 464.42630151666742434726);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 457.89195539777780881528);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 459.17866220189120513538);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 455.01395688191278753824);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 456.30066434664255071302);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 180, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 525.38273677174897945043);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 526.79503714729700591685);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 539.07922565254492455137);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 540.5044905871458019303);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 533.72620553512672882199);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 535.14499439463628212288);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 530.73043363580987294512);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 532.14922300642726468202);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 200, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 603.13928713058383832956);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 604.67882466837180389149);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 617.14333526719735288808);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 618.69532869795791447855);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 611.68946752333086844949);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 613.23523911909772368591);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 608.58834218633722848522);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 610.13411418560930004652);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 9, 220, 1, 1, 4.0);

        Assert_Drude(&test, cint.lnA_TE, 682.84005247565439576229);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 684.50181967017301397707);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 697.12310625261591077832);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 698.79683399531931474978);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 691.57759867858787484769);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 693.24535209716964874026);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 688.38116850409241018627);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 690.04892224665011146115);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 80, 80, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 1000.3265010697618353497);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 1002.9075071673729732209);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1022.4318818735723355838);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1025.0345142889238595537);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 1011.3760316918349128077);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 1013.9678759191470300465);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1011.3760316918349128077);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1013.9678759191470300465);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 100, 100, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 1298.8895358853346095728);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 1301.8685494640090580688);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1321.8922025201483994739);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1324.8893878580164147965);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 1310.3883443583759067664);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 1313.3764623888339949289);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1310.3883443583759067664);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1313.3764623888339949289);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 120, 120, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 1605.7090416488913082738);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 1609.0256882208207635538);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1629.444422099680087274);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1632.7766500102217970946);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 1617.5746301107227581445);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 1620.8990814290859244374);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1617.5746301107227581445);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1620.8990814290859244374);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 140, 140, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 1919.3701648903651043616);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 1922.9787381656604732805);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1943.7247312444170914015);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1947.3469002231727425414);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 1931.5456482310149757212);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 1935.1610302832843647346);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1931.5456482310149757212);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1935.1610302832843647346);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 160, 160, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 2238.8723967764613720782);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 2242.7374833377769533134);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 2263.7630948679067537412);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 2267.6402188935468735333);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 2251.3161722010048739484);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 2255.1872861688314099438);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 2251.3161722010048739484);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 2255.1872861688314099438);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 180, 180, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 2563.470686525182629642);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 2567.5642136168979948949);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 2588.8341146819640037798);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 2592.9384295687322122756);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 2576.1510027527435187158);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 2580.2499307710822658186);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 2576.1510027527435187158);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 2580.2499307710822658186);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 200, 200, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 2892.5886402561011254211);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 2896.8878913140468145497);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 2918.3748105776082988443);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 2922.6838278700639340846);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 2905.4804680514747237778);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 2909.784608025187174863);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 2905.4804680514747237778);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 2909.784608025187174863);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 220, 220, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 3225.7670587947173173257);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 3230.2533173674048736536);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 3251.9355465714229218854);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 3256.4307222979342672647);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 3238.8501601656328614046);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 3243.340882165396897821);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 3238.8501601656328614046);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 3243.340882165396897821);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 250, 250, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 2584.83548639709);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 2588.09269139522);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 2606.91190257296);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 2610.17653310336);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 2595.87269028866);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 2599.13361123991);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 2595.87269028866);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 2599.13361123991);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 300, 300, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 3214.70927052218);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 3218.30944850618);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 3237.51640639121);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 3241.12291037197);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 3226.11200205711);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 3229.71534540588);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 3226.11200205711);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 3229.71534540588);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 350, 350, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 3861.47392585065);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 3865.36917591569);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 3884.89873172236);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 3888.79947895669);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 3873.18561217467);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 3877.08361263686);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 3873.18561217467);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 3877.08361263686);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 400, 400, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 4522.68187003210);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 4526.83557371583);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 4546.64162364176);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 4550.80018113325);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 4534.66112001757);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 4538.81725203163);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 4534.66112001757);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 4538.81725203163);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 450, 450, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 5196.50632803109);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 5200.88971246079);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 5220.93786497312);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 5225.32559107076);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 5208.72153948751);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 5213.10709589997);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 5208.72153948751);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 5213.10709589997);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 500, 500, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 5881.53140000397);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 5886.12132226220);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 5906.38490670702);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 5910.97875416146);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 5893.95765216090);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 5898.54953796066);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 5893.95765216090);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 5898.54953796066);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 550, 550, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 6576.62729898170);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 6581.40477291657);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 6601.86248224733);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 6606.64353653098);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 6589.24443507388);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 6594.02369997087);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 6589.24443507388);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 6594.02369997087);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 600, 600, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 7280.87150394715);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 7285.82069017882);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 7306.45509851717);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 7311.40757513350);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 7293.66288372376);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 7298.61371581495);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 7293.66288372376);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 7298.61371581495);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 700, 700, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 8713.85357201649);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 8719.10788808928);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 8740.05434857256);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 8745.31149508620);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 8726.95360252240);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 8732.20933431268);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 8726.95360252240);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 8732.20933431268);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 600, 800, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 8713.87413933327);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 8719.12845540605);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 8740.05429660210);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 8745.31144311574);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 8727.10770123189);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 8732.36343302217);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 8726.82001915922);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 8732.07575094950);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 200, 50, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 1111.9532800424704341389);
    }

    {
        drude_integrate(&casimir, &cint, 113, 78, 1, 1, 0.001);

        Assert_Drude(&test, cint.lnA_TE, 2101.5903803939422367503);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 2108.8425535639998891074);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 2133.5819127500226934264);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 2140.8551639064542073329);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 2117.3981546718131620359);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 2124.6608945612724042261);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 2117.7688336638360201147);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 2125.0315735532952600116);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 116, 79, 1, 1, 0.005);

        Assert_Drude(&test, cint.lnA_TE, 1839.5275727640859964752);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 1845.2234996259595298362);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 1868.3812780355298324743);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 1874.0977357591566461407);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 1853.7597573714385472809);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 1859.4659758234835122571);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 1854.1438997100805539473);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 1859.8501181621254645226);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 119, 80, 1, 1, 0.01);

        Assert_Drude(&test, cint.lnA_TE, 1745.1807796059782470124);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 1750.2380158647398381132);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1772.7275553412358230707);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1777.8047838588761562616);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 1758.7530754981465532308);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 1763.8203325355563422724);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1759.1501723565945780551);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1764.2174293940041606889);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 122, 81, 1, 1, 0.05);

        Assert_Drude(&test, cint.lnA_TE, 1461.1105990660488679524);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 1464.6984095658989602131);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 1485.5169486434278696991);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 1489.1235236293070864008);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 1473.1064979948426688933);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 1476.7037115113266724543);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 1473.5160698851557277752);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 1477.1132834016350216383);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 125, 82, 1, 1, 0.1);

        Assert_Drude(&test, cint.lnA_TE, 1352.0773458174128010077);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 1355.1193197389409425546);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 1375.1741850141576720794);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 1378.2338235156462245001);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 1363.4125291926091971279);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 1366.4633530684125208715);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 1363.834123683627550851);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 1366.8849475594136416677);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 128, 83, 1, 1, 0.5);

        Assert_Drude(&test, cint.lnA_TE, 1045.5928923853170743104);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 1047.6514029204959882937);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 1065.5479984158034175226);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 1067.6212395198810053819);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 1055.3514647933877059448);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 1057.4173510517892183771);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 1055.7846544732613784604);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 1057.8505407313146176006);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 131, 84, 1, 1, 1.0);

        Assert_Drude(&test, cint.lnA_TE, 921.63791055013489151963);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 923.4688293258289086534);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 940.28072313475524458392);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 942.12516328619051845811);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 930.73478719144083256068);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 932.57247488844250189223);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 931.17916780797379593896);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 933.01685550373454945076);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 134, 85, 1, 1, 5.0);

        Assert_Drude(&test, cint.lnA_TE, 592.37447130269161042673);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 593.95480234791735271119);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 607.87124852432882248041);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 609.46369108406947938655);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 599.89296763881042236829);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 601.47936055457059204878);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 600.3481584154301980271);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 601.93455130426780397453);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 137, 86, 1, 1, 10.0);

        Assert_Drude(&test, cint.lnA_TE, 453.26912531298965194264);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, 454.82266488338732032815);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, 467.44997904216663783209);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, 469.01537997601448383347);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, 460.12446853376366857338);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, 461.68394486255553839154);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, 460.59011081773956695646);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, 462.14958704449861397233);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 141, 87, 1, 1, 50.0);

        Assert_Drude(&test, cint.lnA_TE, 93.056102414849808321233);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, 94.362026556706640857382);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, 104.10158197780293285514);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, 105.42020692508115500005);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, 98.33489616497447796846);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, 99.647183905364570631678);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, 98.817949323927396721432);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, 100.13023462853754440782);
        AssertEqual(&test, cint.signD_TM, -1);
    }

    {
        drude_integrate(&casimir, &cint, 144, 88, 1, 1, 100.0);

        Assert_Drude(&test, cint.lnA_TE, -90.430034469224345470401);
        AssertEqual(&test, cint.signA_TE, 1);
        Assert_Drude(&test, cint.lnA_TM, -89.508537297638933736979);
        AssertEqual(&test, cint.signA_TM, -1);

        Assert_Drude(&test, cint.lnB_TE, -80.707425151280069175637);
        AssertEqual(&test, cint.signB_TE, -1);
        Assert_Drude(&test, cint.lnB_TM, -79.772447064981016034979);
        AssertEqual(&test, cint.signB_TM, 1);

        Assert_Drude(&test, cint.lnC_TE, -85.81822802288151395291);
        AssertEqual(&test, cint.signC_TE, 1);
        Assert_Drude(&test, cint.lnC_TM, -84.889964585118132365301);
        AssertEqual(&test, cint.signC_TM, -1);

        Assert_Drude(&test, cint.lnD_TE, -85.325048992413004957189);
        AssertEqual(&test, cint.signD_TE, -1);
        Assert_Drude(&test, cint.lnD_TM, -84.396793791719501528202);
        AssertEqual(&test, cint.signD_TM, 1);
    }

    {
        drude_integrate(&casimir, &cint, 147, 89, 1, 1, 200.0);

        Assert_Drude(&test, cint.lnA_TE, -346.03971818453141707627);
        AssertEqual(&test, cint.signA_TE, -1);
        Assert_Drude(&test, cint.lnA_TM, -345.59165127117874662804);
        AssertEqual(&test, cint.signA_TM, 1);

        Assert_Drude(&test, cint.lnB_TE, -337.64896680622007463408);
        AssertEqual(&test, cint.signB_TE, 1);
        Assert_Drude(&test, cint.lnB_TM, -337.18850921765651037576);
        AssertEqual(&test, cint.signB_TM, -1);

        Assert_Drude(&test, cint.lnC_TE, -342.10067669190129502604);
        AssertEqual(&test, cint.signC_TE, -1);
        Assert_Drude(&test, cint.lnC_TM, -341.6463707980743521324);
        AssertEqual(&test, cint.signC_TM, 1);

        Assert_Drude(&test, cint.lnD_TE, -341.59675984241556707119);
        AssertEqual(&test, cint.signD_TE, 1);
        Assert_Drude(&test, cint.lnD_TM, -341.14247335882537258888);
        AssertEqual(&test, cint.signD_TM, -1);
    }


    /* test drude integration for perfect reflectors */
    casimir_set_omegap_plane(&casimir, INFINITY);
    casimir_set_gamma_plane(&casimir, 0.0);

    {
        drude_integrate(&casimir, &cint, 1050, 1050, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 13940.125756903571190096123124106829015888756201576644);
        Assert_Drude(&test, cint.lnB_TM, 13967.9514623712062733110267893316781893739350523);
        Assert_Drude(&test, cint.lnC_TM, 13954.03837148533464206);
    }

    {
        drude_integrate(&casimir, &cint, 1050, 1, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 6244.485496882807753109);
        Assert_Drude(&test, cint.lnB_TM, 6263.969792590335805);
        Assert_Drude(&test, cint.lnC_TM, 6250.748896960319286622);
    }

    {
        drude_integrate(&casimir, &cint, 1, 1050, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 6244.485496882807753109081179);
        Assert_Drude(&test, cint.lnB_TM, 6263.969792590335805037405705);
        Assert_Drude(&test, cint.lnC_TM, 6257.7054405868224491877);
    }

    {
        drude_integrate(&casimir, &cint, 500, 1050, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 9813.28363230087848);
        Assert_Drude(&test, cint.lnB_TM, 9839.7598665288218618873847);
        Assert_Drude(&test, cint.lnC_TM, 9826.8923954024132364);
    }

    {
        drude_integrate(&casimir, &cint, 1, 1, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, -2.980829253011726);
        Assert_Drude(&test, cint.lnB_TM, -2.0645385211375711716721);
        Assert_Drude(&test, cint.lnC_TM, -2.5753641449035618548786);
        Assert_Drude(&test, cint.lnD_TM, -2.5753641449035618548786);
    }

    {
        drude_integrate(&casimir, &cint, 241, 73, 1, 1, 30);
        Assert_Drude(&test, cint.lnA_TM, 406.63047665158294437064);
        Assert_Drude(&test, cint.lnB_TM, 419.71230683599700819362);
        Assert_Drude(&test, cint.lnC_TM, 412.57255550309976814896);
        Assert_Drude(&test, cint.lnD_TM, 413.766977852356385781);
    }

    {
        drude_integrate(&casimir, &cint, 241, 1, 1, 1, 30);
        Assert_Drude(&test, cint.lnA_TM, 249.75276347175786475423);
        Assert_Drude(&test, cint.lnB_TM, 258.05248402595679167552);
        Assert_Drude(&test, cint.lnC_TM, 251.17334248392289626321);
        Assert_Drude(&test, cint.lnD_TM, 256.62788585958419530558);
    }

    {
        drude_integrate(&casimir, &cint, 241, 241, 1, 1, 30);
        Assert_Drude(&test, cint.lnA_TM, 838.84852861683729524124);
        Assert_Drude(&test, cint.lnB_TM, 853.98316452183914507246);
        Assert_Drude(&test, cint.lnC_TM, 846.41479992430049881808);
    }

    {
        drude_integrate(&casimir, &cint, 3, 2, 1, 1, 2);
        Assert_Drude(&test, cint.lnA_TM, -4.094372316589062);
        Assert_Drude(&test, cint.lnB_TM, -1.970116759119433);
        Assert_Drude(&test, cint.lnC_TM, -3.298725852652321);
    }

    {
        drude_integrate(&casimir, &cint, 4, 4, 0, 1, 0.005);
        Assert_Drude(&test, cint.lnB_TM, 56.977025325953406);
    }


    {
        drude_integrate(&casimir, &cint, 4, 4, 1, 1, 0.005);
        Assert_Drude(&test, cint.lnA_TM, 40.74560144887208);
        AssertEqual(&test, cint.signA_TM, -1);
        Assert_Drude(&test, cint.lnB_TM, 56.75388164708835);
        AssertEqual(&test, cint.signB_TM, +1);
        Assert_Drude(&test, cint.lnC_TM, 48.68297568585137);
        AssertEqual(&test, cint.signC_TM, -1);
        Assert_Drude(&test, cint.lnD_TM, 48.68297568585137);
        AssertEqual(&test, cint.signD_TM, +1);
    }

    {
        drude_integrate(&casimir, &cint, 40, 40, 1, 1, 0.25);
        Assert_Drude(&test, cint.lnA_TM, 367.258716490769);
        Assert_Drude(&test, cint.lnB_TM, 384.7742430107486);
        Assert_Drude(&test, cint.lnC_TM, 376.010190217081);
        Assert_Drude(&test, cint.lnD_TM, 376.010190217081);
    }

    {
        drude_integrate(&casimir, &cint, 40, 40, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 257.7297499756845);
        Assert_Drude(&test, cint.lnB_TM, 272.472669228606);
        Assert_Drude(&test, cint.lnC_TM, 265.0949179301248);
        Assert_Drude(&test, cint.lnD_TM, 265.0949179301248);
    }

    {
        drude_integrate(&casimir, &cint, 40, 40, 40, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 212.0868844486187);
        Assert_Drude(&test, cint.lnB_TM, 219.4533537701274);
        Assert_Drude(&test, cint.lnC_TM, 215.7638420201849);
        Assert_Drude(&test, cint.lnD_TM, 215.7638420201849);
    }

    {
        drude_integrate(&casimir, &cint, 7, 4, 3, 1, 8);
        Assert_Drude(&test, cint.lnA_TM, -17.1928436859713);
        Assert_Drude(&test, cint.lnB_TM, -16.0865392641165);
        Assert_Drude(&test, cint.lnC_TM, -16.83090135860425);
        Assert_Drude(&test, cint.lnD_TM, -16.48574215820564);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 13);
        Assert_Drude(&test, cint.lnA_TM, -5.202385074993125);
        Assert_Drude(&test, cint.lnB_TM, -0.5773666089005467);
        Assert_Drude(&test, cint.lnC_TM, -2.905312825257782);
        Assert_Drude(&test, cint.lnD_TM, -2.905312825257782);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 0.001);
        Assert_Drude(&test, cint.lnA_TM, 368.3408666279195);
        Assert_Drude(&test, cint.lnB_TM, 391.916763894729);
        Assert_Drude(&test, cint.lnC_TM, 380.1161563573135);
        Assert_Drude(&test, cint.lnD_TM, 380.1161563573135);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 0.01);
        Assert_Drude(&test, cint.lnA_TM, 278.540045473523);
        Assert_Drude(&test, cint.lnB_TM, 297.5107725493697);
        Assert_Drude(&test, cint.lnC_TM, 288.0127501055992);
        Assert_Drude(&test, cint.lnD_TM, 288.0127501055992);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 0.1);
        Assert_Drude(&test, cint.lnA_TM, 188.7389740846992);
        Assert_Drude(&test, cint.lnB_TM, 203.1045304770997);
        Assert_Drude(&test, cint.lnC_TM, 195.9090931914071);
        Assert_Drude(&test, cint.lnD_TM, 195.9090931914071);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 98.91288744548656);
        Assert_Drude(&test, cint.lnB_TM, 108.6732240307656);
        Assert_Drude(&test, cint.lnC_TM, 103.7803782982835);
        Assert_Drude(&test, cint.lnD_TM, 103.7803782982835);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 10);
        Assert_Drude(&test, cint.lnA_TM, 6.660701378819416);
        Assert_Drude(&test, cint.lnB_TM, 11.81202300232528);
        Assert_Drude(&test, cint.lnC_TM, 9.221979692552173);
        Assert_Drude(&test, cint.lnD_TM, 9.221979692552173);
    }

    {
        drude_integrate(&casimir, &cint, 20, 20, 3, 1, 100);
        Assert_Drude(&test, cint.lnA_TM, -201.8065316849248);
        Assert_Drude(&test, cint.lnB_TM, -200.7135958357346);
        Assert_Drude(&test, cint.lnC_TM, -201.2774897833974);
        Assert_Drude(&test, cint.lnD_TM, -201.2774897833974);
    }

    {
        drude_integrate(&casimir, &cint, 20, 10, 3, 1, 0.1);
        Assert_Drude(&test, cint.lnA_TM, 131.0962726931826);
        Assert_Drude(&test, cint.lnB_TM, 144.184735156638);
        Assert_Drude(&test, cint.lnC_TM, 137.2769797021334);
        Assert_Drude(&test, cint.lnD_TM, 137.9701257824486);
    }

    {
        drude_integrate(&casimir, &cint, 20, 15, 10, 1, 5);
        Assert_Drude(&test, cint.lnA_TM, 22.41543637237407);
        Assert_Drude(&test, cint.lnB_TM, 26.04398834917292);
        Assert_Drude(&test, cint.lnC_TM, 24.07623857473362);
        Assert_Drude(&test, cint.lnD_TM, 24.35571686776128);
    }

    {
        drude_integrate(&casimir, &cint, 50, 15, 10, 1, 10);
        Assert_Drude(&test, cint.lnA_TM, 45.60713807155988);
        Assert_Drude(&test, cint.lnB_TM, 49.99651015278684);
        Assert_Drude(&test, cint.lnC_TM, 47.201688624453);
        Assert_Drude(&test, cint.lnD_TM, 48.38666367876067);
    }

    {
        drude_integrate(&casimir, &cint, 100, 25, 20, 1, 20);
        Assert_Drude(&test, cint.lnA_TM, 84.50837701530297);
        Assert_Drude(&test, cint.lnB_TM, 88.66027006552913);
        Assert_Drude(&test, cint.lnC_TM, 85.90224599747752);
        Assert_Drude(&test, cint.lnD_TM, 87.2586869597479);
    }

    {
        drude_integrate(&casimir, &cint, 60, 55, 40, 1, 0.11);
        Assert_Drude(&test, cint.lnA_TM, 645.1730223683922);
        Assert_Drude(&test, cint.lnB_TM, 658.4063308419369);
        Assert_Drude(&test, cint.lnC_TM, 651.7418041758159);
        Assert_Drude(&test, cint.lnD_TM, 651.8288153934492);
    }

    {
        drude_integrate(&casimir, &cint, 40, 40, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 185.27722707813169721211989855051);
        Assert_Drude(&test, cint.lnB_TM, 198.1874611137788);
        Assert_Drude(&test, cint.lnC_TM, 191.72604045861798912);
    }

    {
        drude_integrate(&casimir, &cint, 140, 40, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 575.400220880156994701641252076629);
        Assert_Drude(&test, cint.lnB_TM, 591.1921970497888542);
        Assert_Drude(&test, cint.lnC_TM, 582.6670391327286);
    }

    {
        drude_integrate(&casimir, &cint, 240, 40, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 1025.59108523802829059595981668750);
        Assert_Drude(&test, cint.lnB_TM, 1042.807725206889176);
        Assert_Drude(&test, cint.lnC_TM, 1033.30173547300948);
    }

    {
        drude_integrate(&casimir, &cint, 540, 40, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 2561.62734676892999652846813001817);
        Assert_Drude(&test, cint.lnB_TM, 2581.1132494456470771627);
        Assert_Drude(&test, cint.lnC_TM, 2570.068090210887226712719);
    }

    {
        drude_integrate(&casimir, &cint, 1540, 40, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 8589.22040500894307686493465726593);
        Assert_Drude(&test, cint.lnB_TM, 8611.759673402576546371338);
        Assert_Drude(&test, cint.lnC_TM, 8598.66439349824405764);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 40, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 11616.33666207935);
        Assert_Drude(&test, cint.lnB_TM, 11639.648487984291931);
        Assert_Drude(&test, cint.lnC_TM, 11626.03631835250693323);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 1, 1, 1, 2.5);
        Assert_Drude(&test, cint.lnA_TM, 11358.3678862878775413115);
        Assert_Drude(&test, cint.lnB_TM, 11377.9522208393254813978);
        Assert_Drude(&test, cint.lnC_TM, 11364.359353960757194524);
    }

    {
        drude_integrate(&casimir, &cint, 200, 1, 0, 1, 2.5);
        Assert_Drude(&test, cint.lnB_TM, 682.0020149230455);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 1, 0, 1, 2.5);
        Assert_Drude(&test, cint.lnB_TM, 11378.29904124447803687331);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 40, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 13484.656037918688437);
        Assert_Drude(&test, cint.lnB_TM, 13509.800445322035821193049198451);
        Assert_Drude(&test, cint.lnC_TM, 13495.271984956561221915441695);
    }


    {
        drude_integrate(&casimir, &cint, 2000, 2000, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 29149.71533012066508345912);
        Assert_Drude(&test, cint.lnB_TM, 29180.1186899274219145148);
        Assert_Drude(&test, cint.lnC_TM, 29164.91688500840021851842);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 1000, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 20993.7437127492067941699865154);
        Assert_Drude(&test, cint.lnB_TM, 21022.8784778726214390197843);
        Assert_Drude(&test, cint.lnC_TM, 21007.9643550261185783165);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 1, 1, 1, 0.5);
        Assert_Drude(&test, cint.lnA_TM, 14577.246711903825292880294853452);
        Assert_Drude(&test, cint.lnB_TM, 14600.0499192823994956727755285068);
        Assert_Drude(&test, cint.lnC_TM, 14584.84761448839861741754014567215621);
    }

    {
        drude_integrate(&casimir, &cint, 2000, 1, 1, 1, 10);
        Assert_Drude(&test, cint.lnA_TM, 8585.7322779497324918261);
        Assert_Drude(&test, cint.lnB_TM, 8602.54407061632099700814913677598);
        Assert_Drude(&test, cint.lnC_TM, 8590.3374981457216867851698439);
    }

    {
        drude_integrate(&casimir, &cint, 4000, 1, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 29164.30634417171784698325);
        Assert_Drude(&test, cint.lnB_TM, 29187.80244882461237368);
        Assert_Drude(&test, cint.lnC_TM, 29171.907246756275540666256);
    }

    {
        drude_integrate(&casimir, &cint, 4000, 2000, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 46169.30412583713741863288);
        Assert_Drude(&test, cint.lnB_TM, 46201.211646391476308900991);
        Assert_Drude(&test, cint.lnC_TM, 46184.911229183740231533516827);
    }

    {
        drude_integrate(&casimir, &cint, 4000, 4000, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 63868.6657467722265411811);
        Assert_Drude(&test, cint.lnB_TM, 63901.841820324801967);
        Assert_Drude(&test, cint.lnC_TM, 63885.2537210446057222743554823);
    }

    {
        drude_integrate(&casimir, &cint, 500, 500, 1, 1, 1.006);
        Assert_Drude(&test, cint.lnA_TM, 5880.145418565579706653290311279924286563566945608969682962480558534581786);
        Assert_Drude(&test, cint.lnB_TM, 5904.990886305616635574177976966139243575449897572646531806921581807454247);
        Assert_Drude(&test, cint.lnC_TM, 5892.567652184402156243850461817790605864901360854261291722500861621079491);
        Assert_Drude(&test, cint.lnD_TM, 5892.567652184402156243850461817790605864901360854261291722500861621079491);
    }

    {
        drude_integrate(&casimir, &cint, 499, 500, 1, 1, 1.006);
        Assert_Drude(&test, cint.lnA_TM, 5873.247644850501700848005625848572206294394325680038496156282856954543478);
        Assert_Drude(&test, cint.lnB_TM, 5898.089108585147827540904685136290828084996798178646217185451773987029079);
        Assert_Drude(&test, cint.lnC_TM, 5885.668876966959763457112814499634871743366485689495287378961086374547958);
        Assert_Drude(&test, cint.lnD_TM, 5885.666874964285038289757983818757214367337478685908564744676071164769226);
    }

    {
        drude_integrate(&casimir, &cint, 500, 499, 1, 1, 1.006);
        Assert_Drude(&test, cint.lnA_TM, 5873.247644850501700848005625848572206294394325680038496156282856954543478);
        Assert_Drude(&test, cint.lnB_TM, 5898.089108585147827540904685136290828084996798178646217185451773987029079);
        Assert_Drude(&test, cint.lnC_TM, 5885.666874964285038289757983818757214367337478685908564744676071164769226);
        Assert_Drude(&test, cint.lnD_TM, 5885.668876966959763457112814499634871743366485689495287378961086374547958);
    }

    {
        drude_integrate(&casimir, &cint, 499, 499, 1, 1,  1.006);
        Assert_Drude(&test, cint.lnA_TM, 5866.350873639756594278938154316878502875775073387778089826355287796266225);
        Assert_Drude(&test, cint.lnB_TM, 5891.188331362990374000578938922984773493968146030458032951120864889446638);
        Assert_Drude(&test, cint.lnC_TM, 5878.769101247160551296372865772325533805896297034676801385673068702228523);
        Assert_Drude(&test, cint.lnD_TM, 5878.769101247160551296372865772325533805896297034676801385673068702228523);
    }

    {
        drude_integrate(&casimir, &cint, 200, 200, 1, 1, 1.006);
        Assert_Drude(&test, cint.lnA_TM, 1975.767096926493948364578286189828228822723485260602372040682771555386593);
        Assert_Drude(&test, cint.lnB_TM, 1996.945898960585135363127239781841983882569826371795585966870453220538568);
        Assert_Drude(&test, cint.lnC_TM, 1986.35524636210021458561037266323954886067698626041262792808131074279971);
        Assert_Drude(&test, cint.lnD_TM, 1986.35524636210021458561037266323954886067698626041262792808131074279971);
    }

    {
        drude_integrate(&casimir, &cint, 201, 199, 1, 1,  1.006);
        Assert_Drude(&test, cint.lnA_TM, 1975.76712170967901909630632990248917658766444079226758033464543591258041);
        Assert_Drude(&test, cint.lnB_TM, 1996.945898743456163521589687237776467526903757736254664642401390499554195);
        Assert_Drude(&test, cint.lnC_TM, 1986.350258603304040900412142964201303792493786158794697788127912783088005);
        Assert_Drude(&test, cint.lnD_TM, 1986.360258686952479420258637033592832162969783950835808872636847179108004);
    }

    {
        drude_integrate(&casimir, &cint, 1, 2, 1, 1,  1.006);
        Assert_Drude(&test, cint.lnA_TM, -2.339923455106124005554798782243353095693901455689575867856081014682324586);
        Assert_Drude(&test, cint.lnB_TM, -0.6736624574273509838711334310684649842252665724200226758604799345224393973);
        Assert_Drude(&test, cint.lnC_TM, -1.363077277081885731378935504250430542739743833222871115852865737072863219);
        Assert_Drude(&test, cint.lnD_TM, -1.831883423580260237202605565635544138614431556772498548761012276038582097);
    }

    {
        drude_integrate(&casimir, &cint, 2, 1, 1, 1,  1.006);
        Assert_Drude(&test, cint.lnA_TM, -2.339923455106124005554798782243353095693901455689575867856081014682324586);
        Assert_Drude(&test, cint.lnB_TM, -0.6736624574273509838711334310684649842252665724200226758604799345224393973);
        Assert_Drude(&test, cint.lnC_TM, -1.831883423580260237202605565635544138614431556772498548761012276038582097);
        Assert_Drude(&test, cint.lnD_TM, -1.363077277081885731378935504250430542739743833222871115852865737072863219);
    }

    {
        drude_integrate(&casimir, &cint, 1, 1, 1, 1, 1.006);
        Assert_Drude(&test, cint.lnA_TM, -2.998811324689273716593897834310454621521738584722626358695413310706611303);
        Assert_Drude(&test, cint.lnB_TM, -2.08729623546325568208182858341472166512674251862231861269070482888007109);
        Assert_Drude(&test, cint.lnC_TM, -2.595336266989119347157555830395184063132956853912460002945057793600666948);
        Assert_Drude(&test, cint.lnD_TM, -2.595336266989119347157555830395184063132956853912460002945057793600666948);
    }

    {
        drude_integrate(&casimir, &cint, 1, 500, 1, 1,  1.006);
        Assert_Drude(&test, cint.lnA_TM, 2595.536799689232189132097458491885340008891921765359782557085782932445143);
        Assert_Drude(&test, cint.lnB_TM, 2612.784371554692815301128903227217788398596035743007173804862451217505851);
        Assert_Drude(&test, cint.lnC_TM, 2607.26688661763036080100424831203262202790938673585476359897048874019501);
        Assert_Drude(&test, cint.lnD_TM, 2601.052286639743501207037424084674992544735252498822692375343512556374755);
    }

    {
        drude_integrate(&casimir, &cint, 500, 1, 1, 1,  1.006);
        Assert_Drude(&test, cint.lnA_TM, 2595.536799689232189132097458491885340008891921765359782557085782932445143);
        Assert_Drude(&test, cint.lnB_TM, 2612.784371554692815301128903227217788398596035743007173804862451217505851);
        Assert_Drude(&test, cint.lnC_TM, 2601.052286639743501207037424084674992544735252498822692375343512556374755);
        Assert_Drude(&test, cint.lnD_TM, 2607.26688661763036080100424831203262202790938673585476359897048874019501);
    }

    {
        drude_integrate(&casimir, &cint, 4000 ,4000 ,1 ,1 , 1);
        Assert_Drude(&test, cint.lnA_TM, 63868.6657467722265411811000471834427);
        Assert_Drude(&test, cint.lnB_TM, 63901.84182032480196706998988829691231125220392684855064349853042);
    }

    {
        drude_integrate(&casimir, &cint, 4000 ,2000 ,1 ,1, 1);
        Assert_Drude(&test, cint.lnA_TM, 46169.30412583713741863288035766593016712961550218654);
        Assert_Drude(&test, cint.lnB_TM, 46201.2116463914763089009910622677286368913618141444483663);
        Assert_Drude(&test, cint.lnC_TM, 46184.911229183740231533516827);
    }

    {
        drude_integrate(&casimir, &cint, 4000, 1000, 1, 1, 1);
        Assert_Drude(&test, cint.lnA_TM, 37559.147788784669482290944857175749);
        Assert_Drude(&test, cint.lnC_TM, 37573.8793900544332570436073011074230654516);
    }


    casimir_free(&casimir);

    return test_results(&test, stderr);
}


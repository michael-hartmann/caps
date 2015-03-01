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
#include "test_lbinom.h"
#include "test_logdet.h"
#include "test_mie.h"
#include "test_mie_drude.h"
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
    test_lbinom();
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
    test_logdet(CORES);

    test_casimirF();
    
    return 0;
}

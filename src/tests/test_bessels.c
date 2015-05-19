#include "sfunc.h"
#include "unittest.h"

#include "test_bessels.h"

int test_besselI()
{
    unittest_t test;
    unittest_init(&test, "Bessel function I_nu", "Test modified Bessel function I_nu");

    AssertAlmostEqual(&test, bessel_lnInu(0,1e-6), -7.133546631626697);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-5), -5.982254085113174);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-4), -4.830961536966152);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-3), -3.679668825469134);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-2), -2.528359779027661);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e-1), -1.375417787678169);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e0),  -0.064351991073531);
    AssertAlmostEqual(&test, bessel_lnInu(0,5e0),   3.276297109617906);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e1),   7.929768918237150);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e2),   96.77847637380128);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e3),   995.6271838273042);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e4),   9994.475891280807);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e5),   99993.32459873431);
    AssertAlmostEqual(&test, bessel_lnInu(0,1e6),   999992.1733061878);

    AssertAlmostEqual(&test, bessel_lnInu(1,1e-6), -22.04766947825915);
    AssertAlmostEqual(&test, bessel_lnInu(1,1),    -1.225791352644727);
    AssertAlmostEqual(&test, bessel_lnInu(1,3),     1.131235470744604);
    AssertAlmostEqual(&test, bessel_lnInu(1,1e6),   999992.1733051878);

    AssertAlmostEqual(&test, bessel_lnInu(2,1e-6), -37.47261794865755);
    AssertAlmostEqual(&test, bessel_lnInu(2,1),    -2.862970265776753);
    AssertAlmostEqual(&test, bessel_lnInu(2,5),     2.622265862896675);
    AssertAlmostEqual(&test, bessel_lnInu(2,1e6),   999992.1733031878);

    AssertAlmostEqual(&test, bessel_lnInu(3,1), -4.824473578629219);

    AssertAlmostEqual(&test, bessel_lnInu(23,1e-6), -394.1439513814884);
    AssertAlmostEqual(&test, bessel_lnInu(23,5),    -31.40382021014728);
    AssertAlmostEqual(&test, bessel_lnInu(23,1e6),   999992.1730301876);

    AssertAlmostEqual(&test, bessel_lnInu(119,1e-6), -2189.202200199878);
    AssertAlmostEqual(&test, bessel_lnInu(119,0.5),  -621.0792579289692);
    AssertAlmostEqual(&test, bessel_lnInu(119,3),    -406.9458492626251);
    AssertAlmostEqual(&test, bessel_lnInu(119,30),   -129.9524456900199);
    AssertAlmostEqual(&test, bessel_lnInu(119,300),   272.6929318295042);
    AssertAlmostEqual(&test, bessel_lnInu(119,1e6),   999992.1661661842);

    AssertAlmostEqual(&test, bessel_lnInu(702,1e-6),  -14098.666835519577122094);
    AssertAlmostEqual(&test, bessel_lnInu(702,1e-4),  -10863.534779862939382744);
    AssertAlmostEqual(&test, bessel_lnInu(702,3),     -3621.4923374733442116413);
    AssertAlmostEqual(&test, bessel_lnInu(702,1234),   1034.4300403851143436433);
    AssertAlmostEqual(&test, bessel_lnInu(702,12345),  12319.387046237228462572);

    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-6), -20431.4944983961827997);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-3), -13520.2853417743050538);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-2), -11216.5489562090494164);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e-1), -8912.81256819721365224);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1),    -6609.07593552739598009);
    AssertAlmostEqual(&test, bessel_lnInu(1000,3),    -5509.91234371294526732);
    AssertAlmostEqual(&test, bessel_lnInu(1000,1e3),   527.852986878681152219);

    AssertAlmostEqual(&test, bessel_lnInu(2000,7), -10704.166550337010374);

    return test_results(&test, stderr);
}

int test_besselK()
{
    unittest_t test;
    unittest_init(&test, "Bessel function K_nu", "Test modified Bessel function K_nu");

    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-6), 7.133545631626864);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-5), 5.9822440851298415);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-4), 4.830861538632819);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-3), 3.6786689921357962);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-2), 2.5183764456387734);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e-1), 1.2770838991417504);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e0), -0.7742086473552725);
    AssertAlmostEqual(&test, bessel_lnKnu(0,5e0), -5.5789276035723227);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e1), -10.925501193852295);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e2), -102.07679374034932);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e3), -1003.2280862868463);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e4), -10004.379378833343);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e5), -100005.53067137984);
    AssertAlmostEqual(&test, bessel_lnKnu(0,1e6), -1.0000066819639263e6);

    AssertAlmostEqual(&test, bessel_lnKnu(1,1e-6), 20.949057189590638);
    AssertAlmostEqual(&test, bessel_lnKnu(1,1), -0.08106146679532716);
    AssertAlmostEqual(&test, bessel_lnKnu(1,3),   -3.0358327192375464858953059975);
    AssertAlmostEqual(&test, bessel_lnKnu(1,1e6), -1.00000668196292633e6);

    AssertAlmostEqual(&test, bessel_lnKnu(2,1e-6), 35.863180036223355);
    AssertAlmostEqual(&test, bessel_lnKnu(2,1), 1.17170150170004);
    AssertAlmostEqual(&test, bessel_lnKnu(2,5), -5.036603312746961080665958204772);
    AssertAlmostEqual(&test, bessel_lnKnu(2,1e6), -1.0000066819609263389096196908e6);

    AssertAlmostEqual(&test, bessel_lnKnu(3,1), 2.8367092652889516);

    AssertAlmostEqual(&test, bessel_lnKnu(4,1e15), -1.00000000000001704359684e15);

    AssertAlmostEqual(&test, bessel_lnKnu(23,1e-6), 390.29380377977833);
    AssertAlmostEqual(&test, bessel_lnKnu(23,5), 27.5314997887589672718741222750056);
    AssertAlmostEqual(&test, bessel_lnKnu(23,1e6), -1.00000668168792647542217765299e6);

    AssertAlmostEqual(&test, bessel_lnKnu(119,1e-6), 2183.7257366479472175742539693253862993069);
    AssertAlmostEqual(&test, bessel_lnKnu(119,0.5), 615.6027856231534);
    AssertAlmostEqual(&test, bessel_lnKnu(119,3), 401.4690706673959);
    AssertAlmostEqual(&test, bessel_lnKnu(119,30), 124.44542144141829);
    AssertAlmostEqual(&test, bessel_lnKnu(119,300), -279.16349731660983);

    AssertAlmostEqual(&test, bessel_lnKnu(702,1e-4), 10856.28698728117152647293);
    AssertAlmostEqual(&test, bessel_lnKnu(702,3), 3614.24453577321548255948381274);
    AssertAlmostEqual(&test, bessel_lnKnu(702,1234), -1042.3815655681729711090061175312483747);
    AssertAlmostEqual(&test, bessel_lnKnu(702,12345), -12329.50281632819683895772331427);

    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-6),  20423.89309606159906635627);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-3),  13512.68393943972082096386);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-2),  11208.94755387441573291358);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e-1),  8905.211165857634910126568);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,3),     5502.310936882873879713131);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e3),  -535.8007129753599475405978);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e6),  -1.0000061814642183370632e6);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e10), -1.000000001128708406232e10);
    AssertAlmostEqual(&test, bessel_lnKnu(1000,1e15), -1.000000000000017043596e15);

    return test_results(&test, stderr);
}

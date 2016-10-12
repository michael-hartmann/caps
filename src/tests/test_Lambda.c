#include "libcasimir.h"
#include "unittest.h"

#include "test_Lambda.h"

int test_Lambda(void)
{
    unittest_t test;
    unittest_init(&test, "Lambda", "Test Lambda function for various parameters");

    AssertAlmostEqual(&test, casimir_lnLambda(1,1,0),           0.40546510810816438197801311546434913657);
    AssertAlmostEqual(&test, casimir_lnLambda(1,1,1),          -0.2876820724517809274392190059938274315);

    AssertAlmostEqual(&test, casimir_lnLambda(2,1,1),          -1.1308815492368952772317071947645221687);
    AssertAlmostEqual(&test, casimir_lnLambda(4,5,3),          -10.119213444166137830580776265774926608);
    AssertAlmostEqual(&test, casimir_lnLambda(5,7,3),          -12.079209031704799296645550329264124926);
    AssertAlmostEqual(&test, casimir_lnLambda(16,6,4),         -20.020841061629258759138303870262239174);

    AssertAlmostEqual(&test, casimir_lnLambda(50,50,0),        -3.2287281213112123793793323757932149304);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,1),        -11.072576759463684209642863018499419017);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,50),       -366.96810367687470252345932574544885295);
    AssertAlmostEqual(&test, casimir_lnLambda(50,20,10),       -71.708384125276581706972600369949945466);

    AssertAlmostEqual(&test, casimir_lnLambda(100,1,0),        -1.7557603433310553429384254137655599901);
    AssertAlmostEqual(&test, casimir_lnLambda(100,1,1),        -6.7124792850257034071071320626355070602);
    AssertAlmostEqual(&test, casimir_lnLambda(100,6,4),        -28.190076579450425813590378183217220139);
    AssertAlmostEqual(&test, casimir_lnLambda(100,33,17),      -140.56151146632312845247188373793277075);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,0),      -3.9169857947702750678548639429954691168);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,50),     -460.45932469242092685829036419158558034);
    AssertAlmostEqual(&test, casimir_lnLambda(100,201,10),     -103.40367557445323498326245865777798370);

    AssertAlmostEqual(&test, casimir_lnLambda(200,200,0),      -4.6076608473005432426678531767638037076);
    AssertAlmostEqual(&test, casimir_lnLambda(200,100,70),     -689.79951714054617706762753098998776201);
    AssertAlmostEqual(&test, casimir_lnLambda(500,500,0),      -5.5224594201918359560709544932885934909);

    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,0),    -6.2151077237136242278894880259320868284);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1),    -20.031617782010981865164246152958807085);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,2),    -33.848125842304345491790305669883313001);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,3),    -47.664629906577744997358213920918706752);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,10),   -144.37987862654270719988337448491110499);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,20),   -282.54265122891145618177303437612939551);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,50),   -696.99897104306439045438129052648809049);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,100),  -1387.5321439096580157639308668460135449);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,499),  -6855.7610076931390806110616109722327157);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,500),  -6869.2908341813142468437485944654174625);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,501),  -6882.8193291113699005405673096813816430);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,999),  -13205.138555777978298175816431569216187);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1000), -13212.739458237520380537287638054727456);

    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,0),    -6.6204063732833962735715537724344127507);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,1),    -21.247513592007159661785514225321118969);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,2),    -35.874619922433837419434454683947723733);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,500),  -7301.0165783989268037048388170757359000);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,1000), -14460.235969925823579321977240581831393);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,1500), -21030.645259418831097631857903852345427);

    AssertAlmostEqual(&test, casimir_lnLambda(1500,1,0),    -3.1074706325876159457967703284850318070);

    AssertAlmostEqual(&test, casimir_lnLambda(5000,1,1), -12.573207215572773561556083406674521745);
    AssertAlmostEqual(&test, casimir_lnLambda(5000,500,1), -21.406202989196329346517820278453015337);
    AssertAlmostEqual(&test, casimir_lnLambda(5000,2000,1), -23.484521169044081354765678610825503880);
    AssertAlmostEqual(&test, casimir_lnLambda(5000,5000,1), -24.858732358693766195672947807355216058);

    return test_results(&test, stderr);
}

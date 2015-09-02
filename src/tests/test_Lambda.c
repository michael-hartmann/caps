#include "libcasimir.h"
#include "unittest.h"

#include "test_Lambda.h"

int test_Lambda(void)
{
    unittest_t test;
    unittest_init(&test, "Lambda", "Test Lambda function for various parameters");

    AssertAlmostEqual(&test, casimir_lnLambda(1,1,0,NULL),           0.40546510810816438197801311546434913657);
    AssertAlmostEqual(&test, casimir_lnLambda(1,1,1,NULL),          -0.2876820724517809274392190059938274315);

    AssertAlmostEqual(&test, casimir_lnLambda(2,1,1,NULL),          -1.1308815492368952772317071947645221687);
    AssertAlmostEqual(&test, casimir_lnLambda(4,5,3,NULL),          -10.119213444166137830580776265774926608);
    AssertAlmostEqual(&test, casimir_lnLambda(5,7,3,NULL),          -12.079209031704799296645550329264124926);
    AssertAlmostEqual(&test, casimir_lnLambda(16,6,4,NULL),         -20.020841061629258759138303870262239174);

    AssertAlmostEqual(&test, casimir_lnLambda(50,50,0,NULL),        -3.2287281213112123793793323757932149304);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,1,NULL),        -11.072576759463684209642863018499419017);
    AssertAlmostEqual(&test, casimir_lnLambda(50,50,50,NULL),       -366.96810367687470252345932574544885295);
    AssertAlmostEqual(&test, casimir_lnLambda(50,20,10,NULL),       -71.708384125276581706972600369949945466);

    AssertAlmostEqual(&test, casimir_lnLambda(100,1,0,NULL),        -1.7557603433310553429384254137655599901);
    AssertAlmostEqual(&test, casimir_lnLambda(100,1,1,NULL),        -6.7124792850257034071071320626355070602);
    AssertAlmostEqual(&test, casimir_lnLambda(100,6,4,NULL),        -28.190076579450425813590378183217220139);
    AssertAlmostEqual(&test, casimir_lnLambda(100,33,17,NULL),      -140.56151146632312845247188373793277075);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,0,NULL),      -3.9169857947702750678548639429954691168);
    AssertAlmostEqual(&test, casimir_lnLambda(100,100,50,NULL),     -460.45932469242092685829036419158558034);
    AssertAlmostEqual(&test, casimir_lnLambda(100,201,10,NULL),     -103.40367557445323498326245865777798370);

    AssertAlmostEqual(&test, casimir_lnLambda(200,200,0,NULL),      -4.6076608473005432426678531767638037076);
    AssertAlmostEqual(&test, casimir_lnLambda(200,100,70,NULL),     -689.79951714054617706762753098998776201);
    AssertAlmostEqual(&test, casimir_lnLambda(500,500,0,NULL),      -5.5224594201918359560709544932885934909);

    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,0,NULL),    -6.2151077237136242278894880259320868284);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1,NULL),    -20.031617782010981865164246152958807085);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,2,NULL),    -33.848125842304345491790305669883313001);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,3,NULL),    -47.664629906577744997358213920918706752);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,10,NULL),   -144.37987862654270719988337448491110499);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,20,NULL),   -282.54265122891145618177303437612939551);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,50,NULL),   -696.99897104306439045438129052648809049);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,100,NULL),  -1387.5321439096580157639308668460135449);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,499,NULL),  -6855.7610076931390806110616109722327157);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,500,NULL),  -6869.2908341813142468437485944654174625);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,501,NULL),  -6882.8193291113699005405673096813816430);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,999,NULL),  -13205.138555777978298175816431569216187);
    AssertAlmostEqual(&test, casimir_lnLambda(1000,1000,1000,NULL), -13212.739458237520380537287638054727456);

    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,0,NULL),    -6.6204063732833962735715537724344127507);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,1,NULL),    -21.247513592007159661785514225321118969);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,2,NULL),    -35.874619922433837419434454683947723733);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,500,NULL),  -7301.0165783989268037048388170757359000);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,1000,NULL), -14460.235969925823579321977240581831393);
    AssertAlmostEqual(&test, casimir_lnLambda(1500,1500,1500,NULL), -21030.645259418831097631857903852345427);

    AssertAlmostEqual(&test, casimir_lnLambda(1500,1,0,NULL),    -3.1074706325876159457967703284850318070);

    return test_results(&test, stderr);
}

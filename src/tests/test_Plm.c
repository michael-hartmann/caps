#include <math.h>
#include "unittest.h"
#include "sfunc.h"
#include "libcasimir.h"

#include "test_Plm.h"

int test_Plm()
{
    unittest_t test;

    unittest_init(&test, "Plm", "Test associated Legendre polynomials");

    AssertAlmostEqual(&test, Plm(1,0,1.1,1,1), 0.09531017980432493);
    AssertAlmostEqual(&test, Plm(1,0,2,1,1), 0.6931471805599453);
    AssertAlmostEqual(&test, Plm(1,0,5,1,1), 1.6094379124341003);
    AssertAlmostEqual(&test, Plm(1,0,10,1,1), 2.3025850929940455);
    AssertAlmostEqual(&test, Plm(1,0,100,1,1), 4.605170185988091);
    AssertAlmostEqual(&test, Plm(1,0,1000,1,1), 6.907755278982137);
    AssertAlmostEqual(&test, Plm(1,0,1e+06,1,1), 13.815510557964274);
    AssertAlmostEqual(&test, Plm(1,1,1.1,1,1), -0.7803238741323337);
    AssertAlmostEqual(&test, Plm(1,1,2,1,1), 0.5493061443340548);
    AssertAlmostEqual(&test, Plm(1,1,5,1,1), 1.5890269151739727);
    AssertAlmostEqual(&test, Plm(1,1,10,1,1), 2.297559925067295);
    AssertAlmostEqual(&test, Plm(1,1,100,1,1), 4.605120183487924);
    AssertAlmostEqual(&test, Plm(1,1,1000,1,1), 6.907754778981887);
    AssertAlmostEqual(&test, Plm(1,1,1e+06,1,1), 13.815510557963773);
    AssertAlmostEqual(&test, Plm(5,0,1.1,1,1), 1.1310847224188405);
    AssertAlmostEqual(&test, Plm(5,0,2,1,1), 5.224401683597868);
    AssertAlmostEqual(&test, Plm(5,0,5,1,1), 10.06581896445358);
    AssertAlmostEqual(&test, Plm(5,0,10,1,1), 13.5654694258405);
    AssertAlmostEqual(&test, Plm(5,0,100,1,1), 25.08943299974896);
    AssertAlmostEqual(&test, Plm(5,0,1000,1,1), 36.602468468510885);
    AssertAlmostEqual(&test, Plm(5,0,1e+06,1,1), 71.14124597453196);
    AssertAlmostEqual(&test, Plm(5,1,1.1,1,1), 2.54332404330678);
    AssertAlmostEqual(&test, Plm(5,1,2,1,1), 6.816269473090174);
    AssertAlmostEqual(&test, Plm(5,1,5,1,1), 11.672959264491194);
    AssertAlmostEqual(&test, Plm(5,1,10,1,1), 15.17434719989084);
    AssertAlmostEqual(&test, Plm(5,1,100,1,1), 26.69886535617332);
    AssertAlmostEqual(&test, Plm(5,1,1000,1,1), 38.211906325389386);
    AssertAlmostEqual(&test, Plm(5,1,1e+06,1,1), 72.75068388696599);
    AssertAlmostEqual(&test, Plm(5,5,1.1,1,1), 2.949565556832074);
    AssertAlmostEqual(&test, Plm(5,5,2,1,1), 9.597715649164016);
    AssertAlmostEqual(&test, Plm(5,5,5,1,1), 14.796319503363605);
    AssertAlmostEqual(&test, Plm(5,5,10,1,1), 18.338984552830215);
    AssertAlmostEqual(&test, Plm(5,5,100,1,1), 29.876785844933366);
    AssertAlmostEqual(&test, Plm(5,5,1000,1,1), 41.38995882240317);
    AssertAlmostEqual(&test, Plm(5,5,1e+06,1,1), 75.9287377173126);
    AssertAlmostEqual(&test, Plm(10,0,1.1,1,1), 2.9853524385054198);
    AssertAlmostEqual(&test, Plm(10,0,2,1,1), 11.47273368360779);
    AssertAlmostEqual(&test, Plm(10,0,5,1,1), 21.193564845688325);
    AssertAlmostEqual(&test, Plm(10,0,10,1,1), 28.19740034187603);
    AssertAlmostEqual(&test, Plm(10,0,100,1,1), 51.24678451823574);
    AssertAlmostEqual(&test, Plm(10,0,1000,1,1), 74.27286993040246);
    AssertAlmostEqual(&test, Plm(10,0,1e+06,1,1), 143.35042508864336);
    AssertAlmostEqual(&test, Plm(10,1,1.1,1,1), 5.206849812025071);
    AssertAlmostEqual(&test, Plm(10,1,2,1,1), 13.767103131424927);
    AssertAlmostEqual(&test, Plm(10,1,5,1,1), 23.495063352349288);
    AssertAlmostEqual(&test, Plm(10,1,10,1,1), 30.499720210145995);
    AssertAlmostEqual(&test, Plm(10,1,100,1,1), 53.54936697944592);
    AssertAlmostEqual(&test, Plm(10,1,1000,1,1), 76.5754549970807);
    AssertAlmostEqual(&test, Plm(10,1,1e+06,1,1), 145.65301018163737);
    AssertAlmostEqual(&test, Plm(10,5,1.1,1,1), 11.570910101939365);
    AssertAlmostEqual(&test, Plm(10,5,2,1,1), 21.586765038605602);
    AssertAlmostEqual(&test, Plm(10,5,5,1,1), 31.483362977337464);
    AssertAlmostEqual(&test, Plm(10,5,10,1,1), 38.50769304102551);
    AssertAlmostEqual(&test, Plm(10,5,100,1,1), 61.56363955417693);
    AssertAlmostEqual(&test, Plm(10,5,1000,1,1), 84.5897901028007);
    AssertAlmostEqual(&test, Plm(10,5,1e+06,1,1), 153.66734591893618);
    AssertAlmostEqual(&test, Plm(10,10,1.1,1,1), 12.496493340755178);
    AssertAlmostEqual(&test, Plm(10,10,2,1,1), 25.792793525419064);
    AssertAlmostEqual(&test, Plm(10,10,5,1,1), 36.19000123381824);
    AssertAlmostEqual(&test, Plm(10,10,10,1,1), 43.27533133275146);
    AssertAlmostEqual(&test, Plm(10,10,100,1,1), 66.35093391695776);
    AssertAlmostEqual(&test, Plm(10,10,1000,1,1), 89.37727987189739);
    AssertAlmostEqual(&test, Plm(10,10,1e+06,1,1), 158.45483766171625);
    AssertAlmostEqual(&test, Plm(50,0,1.1,1,1), 19.916497983164188);
    AssertAlmostEqual(&test, Plm(50,0,2,1,1), 63.3546620129539);
    AssertAlmostEqual(&test, Plm(50,0,5,1,1), 112.09588794463903);
    AssertAlmostEqual(&test, Plm(50,0,10,1,1), 147.13153649500967);
    AssertAlmostEqual(&test, Plm(50,0,100,1,1), 262.3837545036006);
    AssertAlmostEqual(&test, Plm(50,0,1000,1,1), 377.514234199385);
    AssertAlmostEqual(&test, Plm(50,0,1e+06,1,1), 722.9020105222214);
    AssertAlmostEqual(&test, Plm(50,1,1.1,1,1), 23.81416668339634);
    AssertAlmostEqual(&test, Plm(50,1,2,1,1), 67.26511989971114);
    AssertAlmostEqual(&test, Plm(50,1,5,1,1), 116.00770261583082);
    AssertAlmostEqual(&test, Plm(50,1,10,1,1), 151.0435086107852);
    AssertAlmostEqual(&test, Plm(50,1,100,1,1), 266.2957770039401);
    AssertAlmostEqual(&test, Plm(50,1,1000,1,1), 381.42625719976263);
    AssertAlmostEqual(&test, Plm(50,1,1e+06,1,1), 726.8140335276496);
    AssertAlmostEqual(&test, Plm(50,5,1.1,1,1), 38.91341611149594);
    AssertAlmostEqual(&test, Plm(50,5,2,1,1), 82.66938336265711);
    AssertAlmostEqual(&test, Plm(50,5,5,1,1), 131.44451317185312);
    AssertAlmostEqual(&test, Plm(50,5,10,1,1), 166.4840975823886);
    AssertAlmostEqual(&test, Plm(50,5,100,1,1), 281.7375751890318);
    AssertAlmostEqual(&test, Plm(50,5,1000,1,1), 396.8680673857681);
    AssertAlmostEqual(&test, Plm(50,5,1e+06,1,1), 742.2558438348672);
    AssertAlmostEqual(&test, Plm(50,10,1.1,1,1), 56.6687112091108);
    AssertAlmostEqual(&test, Plm(50,10,2,1,1), 101.35554051734911);
    AssertAlmostEqual(&test, Plm(50,10,5,1,1), 150.23218621891817);
    AssertAlmostEqual(&test, Plm(50,10,10,1,1), 185.28357503153336);
    AssertAlmostEqual(&test, Plm(50,10,100,1,1), 300.54083123177764);
    AssertAlmostEqual(&test, Plm(50,10,1000,1,1), 415.67136093134974);
    AssertAlmostEqual(&test, Plm(50,10,1e+06,1,1), 761.0591377592366);
    AssertAlmostEqual(&test, Plm(50,50,1.1,1,1), 141.5880558691765);
    AssertAlmostEqual(&test, Plm(50,50,2,1,1), 208.06955679249592);
    AssertAlmostEqual(&test, Plm(50,50,5,1,1), 260.0555953344918);
    AssertAlmostEqual(&test, Plm(50,50,10,1,1), 295.48224582915793);
    AssertAlmostEqual(&test, Plm(50,50,100,1,1), 410.8602587501894);
    AssertAlmostEqual(&test, Plm(50,50,1000,1,1), 525.9919885248875);
    AssertAlmostEqual(&test, Plm(50,50,1e+06,1,1), 871.3797774739818);
    AssertAlmostEqual(&test, Plm(100,0,1.1,1,1), 41.74776988801881);
    AssertAlmostEqual(&test, Plm(100,0,2,1,1), 128.85703644825924);
    AssertAlmostEqual(&test, Plm(100,0,5,1,1), 226.37212162047132);
    AssertAlmostEqual(&test, Plm(100,0,10,1,1), 296.44734878125274);
    AssertAlmostEqual(&test, Plm(100,0,100,1,1), 526.9530490939379);
    AssertAlmostEqual(&test, Plm(100,0,1000,1,1), 757.2140210491165);
    AssertAlmostEqual(&test, Plm(100,0,1e+06,1,1), 1447.9895738216865);
    AssertAlmostEqual(&test, Plm(100,1,1.1,1,1), 46.3458523370313);
    AssertAlmostEqual(&test, Plm(100,1,2,1,1), 133.46142863478607);
    AssertAlmostEqual(&test, Plm(100,1,5,1,1), 230.97718817389858);
    AssertAlmostEqual(&test, Plm(100,1,10,1,1), 301.0524936509405);
    AssertAlmostEqual(&test, Plm(100,1,100,1,1), 531.5582190286508);
    AssertAlmostEqual(&test, Plm(100,1,1000,1,1), 761.819191232592);
    AssertAlmostEqual(&test, Plm(100,1,1e+06,1,1), 1452.5947440076745);
    AssertAlmostEqual(&test, Plm(100,5,1.1,1,1), 64.49512261401553);
    AssertAlmostEqual(&test, Plm(100,5,2,1,1), 151.76190508973625);
    AssertAlmostEqual(&test, Plm(100,5,5,1,1), 249.29384752416067);
    AssertAlmostEqual(&test, Plm(100,5,10,1,1), 319.3710325604637);
    AssertAlmostEqual(&test, Plm(100,5,100,1,1), 549.8773594968228);
    AssertAlmostEqual(&test, Plm(100,5,1000,1,1), 780.1383376710668);
    AssertAlmostEqual(&test, Plm(100,5,1e+06,1,1), 1470.9138905064508);
    AssertAlmostEqual(&test, Plm(100,10,1.1,1,1), 86.62947342742189);
    AssertAlmostEqual(&test, Plm(100,10,2,1,1), 174.36600447369366);
    AssertAlmostEqual(&test, Plm(100,10,5,1,1), 271.94849484328086);
    AssertAlmostEqual(&test, Plm(100,10,10,1,1), 342.03155311883495);
    AssertAlmostEqual(&test, Plm(100,10,100,1,1), 572.5397599017762);
    AssertAlmostEqual(&test, Plm(100,10,1000,1,1), 802.8007567332141);
    AssertAlmostEqual(&test, Plm(100,10,1e+06,1,1), 1493.5763097570402);
    AssertAlmostEqual(&test, Plm(100,50,1.1,1,1), 241.06115080063847);
    AssertAlmostEqual(&test, Plm(100,50,2,1,1), 342.19333262236023);
    AssertAlmostEqual(&test, Plm(100,50,5,1,1), 441.3749906552409);
    AssertAlmostEqual(&test, Plm(100,50,10,1,1), 511.6456869778622);
    AssertAlmostEqual(&test, Plm(100,50,100,1,1), 742.2140295117537);
    AssertAlmostEqual(&test, Plm(100,50,1000,1,1), 972.4756233714953);
    AssertAlmostEqual(&test, Plm(100,50,1e+06,1,1), 1663.2511824254705);
    AssertAlmostEqual(&test, Plm(100,100,1.1,1,1), 352.14550616761403);
    AssertAlmostEqual(&test, Plm(100,100,2,1,1), 485.1085080142529);
    AssertAlmostEqual(&test, Plm(100,100,5,1,1), 589.0805850982447);
    AssertAlmostEqual(&test, Plm(100,100,10,1,1), 659.9338860875769);
    AssertAlmostEqual(&test, Plm(100,100,100,1,1), 890.6899119296398);
    AssertAlmostEqual(&test, Plm(100,100,1000,1,1), 1120.9533714790362);
    AssertAlmostEqual(&test, Plm(100,100,1e+06,1,1), 1811.7289493772248);
    AssertAlmostEqual(&test, Plm(500,0,1.1,1,1), 218.369931617909);
    AssertAlmostEqual(&test, Plm(500,0,2,1,1), 654.8363204761886);
    AssertAlmostEqual(&test, Plm(500,0,5,1,1), 1142.5410497358569);
    AssertAlmostEqual(&test, Plm(500,0,10,1,1), 1492.9327632020668);
    AssertAlmostEqual(&test, Plm(500,0,100,1,1), 2645.4662763264428);
    AssertAlmostEqual(&test, Plm(500,0,1000,1,1), 3796.771185904025);
    AssertAlmostEqual(&test, Plm(500,0,1e+06,1,1), 7250.6489502698905);
    AssertAlmostEqual(&test, Plm(500,1,1.1,1,1), 224.58313594464965);
    AssertAlmostEqual(&test, Plm(500,1,2,1,1), 661.0507736951938);
    AssertAlmostEqual(&test, Plm(500,1,5,1,1), 1148.7556371924848);
    AssertAlmostEqual(&test, Plm(500,1,10,1,1), 1499.1473662576052);
    AssertAlmostEqual(&test, Plm(500,1,100,1,1), 2651.6808843748113);
    AssertAlmostEqual(&test, Plm(500,1,1000,1,1), 3802.985794001947);
    AssertAlmostEqual(&test, Plm(500,1,1e+06,1,1), 7256.863558368313);
    AssertAlmostEqual(&test, Plm(500,5,1.1,1,1), 249.38781929785907);
    AssertAlmostEqual(&test, Plm(500,5,2,1,1), 685.8854287299744);
    AssertAlmostEqual(&test, Plm(500,5,5,1,1), 1173.5935139152787);
    AssertAlmostEqual(&test, Plm(500,5,10,1,1), 1523.9856173540122);
    AssertAlmostEqual(&test, Plm(500,5,100,1,1), 2676.5192552991175);
    AssertAlmostEqual(&test, Plm(500,5,1000,1,1), 3827.8241661155325);
    AssertAlmostEqual(&test, Plm(500,5,1e+06,1,1), 7281.70193049391);
    AssertAlmostEqual(&test, Plm(500,10,1.1,1,1), 280.2850888002344);
    AssertAlmostEqual(&test, Plm(500,10,2,1,1), 716.876338307201);
    AssertAlmostEqual(&test, Plm(500,10,5,1,1), 1204.5944910827982);
    AssertAlmostEqual(&test, Plm(500,10,10,1,1), 1554.987764436067);
    AssertAlmostEqual(&test, Plm(500,10,100,1,1), 2707.521776843168);
    AssertAlmostEqual(&test, Plm(500,10,1000,1,1), 3858.8266913760813);
    AssertAlmostEqual(&test, Plm(500,10,1e+06,1,1), 7312.704455791996);
    AssertAlmostEqual(&test, Plm(500,50,1.1,1,1), 523.073734528164);
    AssertAlmostEqual(&test, Plm(500,50,2,1,1), 962.6445786135785);
    AssertAlmostEqual(&test, Plm(500,50,5,1,1), 1450.6847466448119);
    AssertAlmostEqual(&test, Plm(500,50,10,1,1), 1801.1154548591724);
    AssertAlmostEqual(&test, Plm(500,50,100,1,1), 2953.6614498983918);
    AssertAlmostEqual(&test, Plm(500,50,1000,1,1), 4104.966483359222);
    AssertAlmostEqual(&test, Plm(500,50,1e+06,1,1), 7558.844248976339);
    AssertAlmostEqual(&test, Plm(500,100,1.1,1,1), 815.4370579952206);
    AssertAlmostEqual(&test, Plm(500,100,2,1,1), 1264.1198067814141);
    AssertAlmostEqual(&test, Plm(500,100,5,1,1), 1753.1644352752814);
    AssertAlmostEqual(&test, Plm(500,100,10,1,1), 2103.7120973975634);
    AssertAlmostEqual(&test, Plm(500,100,100,1,1), 3256.2955362655425);
    AssertAlmostEqual(&test, Plm(500,100,1000,1,1), 4407.600941375931);
    AssertAlmostEqual(&test, Plm(500,100,1e+06,1,1), 7861.4787107468);
    // AssertAlmostEqual(&test, Plm(500,500,1.1,1,1), ???);
    AssertAlmostEqual(&test, Plm(500,500,2,1,1), 3228.877201915062);
    AssertAlmostEqual(&test, Plm(500,500,5,1,1), 3748.7375873350206);
    AssertAlmostEqual(&test, Plm(500,500,10,1,1), 4103.004092281682);
    AssertAlmostEqual(&test, Plm(500,500,100,1,1), 5256.7842214919965);
    AssertAlmostEqual(&test, Plm(500,500,1000,1,1), 6408.101519238978);
    AssertAlmostEqual(&test, Plm(500,500,1e+06,1,1), 9861.979408729921);
    AssertAlmostEqual(&test, Plm(1000,0,1.1,1,1), 439.8074345368616);
    AssertAlmostEqual(&test, Plm(1000,0,2,1,1), 1312.9688009771617);
    AssertAlmostEqual(&test, Plm(1000,0,5,1,1), 2288.4104333445875);
    AssertAlmostEqual(&test, Plm(1000,0,10,1,1), 2989.1977370442632);
    AssertAlmostEqual(&test, Plm(1000,0,100,1,1), 5294.266010535109);
    AssertAlmostEqual(&test, Plm(1000,0,1000,1,1), 7596.87584208464);
    AssertAlmostEqual(&test, Plm(1000,0,1e+06,1,1), 14504.631370941557);
    AssertAlmostEqual(&test, Plm(1000,1,1.1,1,1), 446.71448877540297);
    AssertAlmostEqual(&test, Plm(1000,1,2,1,1), 1319.8764788611857);
    AssertAlmostEqual(&test, Plm(1000,1,5,1,1), 2295.318178307942);
    AssertAlmostEqual(&test, Plm(1000,1,10,1,1), 2996.105489803071);
    AssertAlmostEqual(&test, Plm(1000,1,100,1,1), 5301.173765789076);
    AssertAlmostEqual(&test, Plm(1000,1,1000,1,1), 7603.783597363372);
    AssertAlmostEqual(&test, Plm(1000,1,1e+06,1,1), 14511.53912622054);
    AssertAlmostEqual(&test, Plm(1000,5,1.1,1,1), 474.31867010444955);
    AssertAlmostEqual(&test, Plm(1000,5,2,1,1), 1347.4956274665913);
    AssertAlmostEqual(&test, Plm(1000,5,5,1,1), 2322.9389368154248);
    AssertAlmostEqual(&test, Plm(1000,5,10,1,1), 3023.726435401402);
    AssertAlmostEqual(&test, Plm(1000,5,100,1,1), 5328.794771271238);
    AssertAlmostEqual(&test, Plm(1000,5,1000,1,1), 7631.404603439875);
    AssertAlmostEqual(&test, Plm(1000,5,1e+06,1,1), 14539.160132303046);
    AssertAlmostEqual(&test, Plm(1000,10,1.1,1,1), 508.76974368580466);
    AssertAlmostEqual(&test, Plm(1000,10,2,1,1), 1381.993471123552);
    AssertAlmostEqual(&test, Plm(1000,10,5,1,1), 2357.4418113933493);
    AssertAlmostEqual(&test, Plm(1000,10,10,1,1), 3058.2298946378514);
    AssertAlmostEqual(&test, Plm(1000,10,100,1,1), 5363.298417644634);
    AssertAlmostEqual(&test, Plm(1000,10,1000,1,1), 7665.908251670591);
    AssertAlmostEqual(&test, Plm(1000,10,1e+06,1,1), 14573.663780552522);
    AssertAlmostEqual(&test, Plm(1000,50,1.1,1,1), 782.1991216934603);
    AssertAlmostEqual(&test, Plm(1000,50,2,1,1), 1656.9173691658418);
    AssertAlmostEqual(&test, Plm(1000,50,5,1,1), 2632.526680489959);
    AssertAlmostEqual(&test, Plm(1000,50,10,1,1), 3333.3334725075015);
    AssertAlmostEqual(&test, Plm(1000,50,100,1,1), 5638.407983877662);
    AssertAlmostEqual(&test, Plm(1000,50,1000,1,1), 7941.017877337837);
    AssertAlmostEqual(&test, Plm(1000,50,1e+06,1,1), 14848.773406820066);
    AssertAlmostEqual(&test, Plm(1000,100,1.1,1,1), 1118.4853446625614);
    AssertAlmostEqual(&test, Plm(1000,100,2,1,1), 1997.8478295053453);
    AssertAlmostEqual(&test, Plm(1000,100,5,1,1), 2973.9599454242243);
    AssertAlmostEqual(&test, Plm(1000,100,10,1,1), 3674.8251986109312);
    AssertAlmostEqual(&test, Plm(1000,100,100,1,1), 5979.918423380138);
    AssertAlmostEqual(&test, Plm(1000,100,1000,1,1), 8282.528502572215);
    AssertAlmostEqual(&test, Plm(1000,100,1e+06,1,1), 15190.284033930384);
    AssertAlmostEqual(&test, Plm(1000,500,1.1,1,1), 3582.309249564578);
    AssertAlmostEqual(&test, Plm(1000,500,2,1,1), 4594.6097374355295);
    AssertAlmostEqual(&test, Plm(1000,500,5,1,1), 5586.632588306393);
    AssertAlmostEqual(&test, Plm(1000,500,10,1,1), 6289.365612497314);
    AssertAlmostEqual(&test, Plm(1000,500,100,1,1), 8595.057476986954);
    AssertAlmostEqual(&test, Plm(1000,500,1000,1,1), 10897.673499581337);
    AssertAlmostEqual(&test, Plm(1000,500,1e+06,1,1), 17805.4290909695);
    // AssertAlmostEqual(&test, Plm(1000,1000,1.1,1,1), ???);
    AssertAlmostEqual(&test, Plm(1000,1000,2,1,1), 7150.555135799753);
    AssertAlmostEqual(&test, Plm(1000,1000,5,1,1), 8190.27590663967);
    AssertAlmostEqual(&test, Plm(1000,1000,10,1,1), 8898.808916532993);
    AssertAlmostEqual(&test, Plm(1000,1000,100,1,1), 11206.369174953621);
    AssertAlmostEqual(&test, Plm(1000,1000,1000,1,1), 13509.003770447584);
    AssertAlmostEqual(&test, Plm(1000,1000,1e+06,1,1), 20416.75954942947);
    AssertAlmostEqual(&test, Plm(2000,0,1.1,1,1), 883.0290901488979);
    AssertAlmostEqual(&test, Plm(2000,0,2,1,1), 2629.5801771345295);
    AssertAlmostEqual(&test, Plm(2000,0,5,1,1), 4580.4955905256975);
    AssertAlmostEqual(&test, Plm(2000,0,10,1,1), 5982.074071765259);
    AssertAlmostEqual(&test, Plm(2000,0,100,1,1), 10592.21186505218);
    AssertAlmostEqual(&test, Plm(2000,0,1000,1,1), 15197.431540536312);
    AssertAlmostEqual(&test, Plm(2000,0,1e+06,1,1), 29012.94259837524);
    AssertAlmostEqual(&test, Plm(2000,1,1.1,1,1), 890.6296422989395);
    AssertAlmostEqual(&test, Plm(2000,1,2,1,1), 2637.1810409077684);
    AssertAlmostEqual(&test, Plm(2000,1,5,1,1), 4588.096487828742);
    AssertAlmostEqual(&test, Plm(2000,1,10,1,1), 5989.67497296503);
    AssertAlmostEqual(&test, Plm(2000,1,100,1,1), 10599.81276749922);
    AssertAlmostEqual(&test, Plm(2000,1,1000,1,1), 15205.032442995729);
    AssertAlmostEqual(&test, Plm(2000,1,1e+06,1,1), 29020.543500834785);
    AssertAlmostEqual(&test, Plm(2000,5,1.1,1,1), 921.019840981977);
    AssertAlmostEqual(&test, Plm(2000,5,2,1,1), 2667.578718520729);
    AssertAlmostEqual(&test, Plm(2000,5,5,1,1), 4618.494970156808);
    AssertAlmostEqual(&test, Plm(2000,5,10,1,1), 6020.073548814538);
    AssertAlmostEqual(&test, Plm(2000,5,100,1,1), 10630.211373283117);
    AssertAlmostEqual(&test, Plm(2000,5,1000,1,1), 15235.431049076726);
    AssertAlmostEqual(&test, Plm(2000,5,1e+06,1,1), 29050.94210691878);
    AssertAlmostEqual(&test, Plm(2000,10,1.1,1,1), 958.9805485310268);
    AssertAlmostEqual(&test, Plm(2000,10,2,1,1), 2705.5627973939354);
    AssertAlmostEqual(&test, Plm(2000,10,5,1,1), 4656.481563761844);
    AssertAlmostEqual(&test, Plm(2000,10,10,1,1), 6058.060434674034);
    AssertAlmostEqual(&test, Plm(2000,10,100,1,1), 10668.19835268758);
    AssertAlmostEqual(&test, Plm(2000,10,1000,1,1), 15273.418029409615);
    AssertAlmostEqual(&test, Plm(2000,10,1e+06,1,1), 29088.929087261047);
    AssertAlmostEqual(&test, Plm(2000,50,1.1,1,1), 1261.5811042164298);
    AssertAlmostEqual(&test, Plm(2000,50,2,1,1), 3008.9109702192545);
    AssertAlmostEqual(&test, Plm(2000,50,5,1,1), 4959.910205706487);
    AssertAlmostEqual(&test, Plm(2000,50,10,1,1), 6361.498428723967);
    AssertAlmostEqual(&test, Plm(2000,50,100,1,1), 10971.639340174106);
    AssertAlmostEqual(&test, Plm(2000,50,1000,1,1), 15576.85904660582);
    AssertAlmostEqual(&test, Plm(2000,50,1e+06,1,1), 29392.370104757327);
    AssertAlmostEqual(&test, Plm(2000,100,1.1,1,1), 1637.1036366485025);
    AssertAlmostEqual(&test, Plm(2000,100,2,1,1), 3386.7665032145783);
    AssertAlmostEqual(&test, Plm(2000,100,5,1,1), 5338.017175967623);
    AssertAlmostEqual(&test, Plm(2000,100,10,1,1), 6739.634623846346);
    AssertAlmostEqual(&test, Plm(2000,100,100,1,1), 11349.784889756313);
    AssertAlmostEqual(&test, Plm(2000,100,1000,1,1), 15955.004689030764);
    AssertAlmostEqual(&test, Plm(2000,100,1e+06,1,1), 29770.515748120008);
    AssertAlmostEqual(&test, Plm(2000,500,1.1,1,1), 4530.154954240339);
    AssertAlmostEqual(&test, Plm(2000,500,2,1,1), 6352.051211695172);
    AssertAlmostEqual(&test, Plm(2000,500,5,1,1), 8311.325050080186);
    AssertAlmostEqual(&test, Plm(2000,500,10,1,1), 9713.877319627418);
    AssertAlmostEqual(&test, Plm(2000,500,100,1,1), 14324.326904634738);
    AssertAlmostEqual(&test, Plm(2000,500,1000,1,1), 18929.549674874524);
    AssertAlmostEqual(&test, Plm(2000,500,1e+06,1,1), 32745.06076397126);
    AssertAlmostEqual(&test, Plm(2000,1000,1.1,1,1), 7860.963699943079);
    AssertAlmostEqual(&test, Plm(2000,1000,2,1,1), 9885.673416917001);
    AssertAlmostEqual(&test, Plm(2000,1000,5,1,1), 11869.74194034297);
    AssertAlmostEqual(&test, Plm(2000,1000,10,1,1), 13275.210870832698);
    AssertAlmostEqual(&test, Plm(2000,1000,100,1,1), 17886.59553305319);
    AssertAlmostEqual(&test, Plm(2000,1000,1000,1,1), 22491.827587530606);
    AssertAlmostEqual(&test, Plm(2000,1000,1e+06,1,1), 36307.33877040076);
    AssertAlmostEqual(&test, Plm(5000,0,1.1,1,1), 2213.27569287925);
    AssertAlmostEqual(&test, Plm(5000,0,2,1,1), 6579.995754239426);
    AssertAlmostEqual(&test, Plm(5000,0,5,1,1), 11457.33249056974);
    AssertAlmostEqual(&test, Plm(5000,0,10,1,1), 14961.284502089478);
    AssertAlmostEqual(&test, Plm(5000,0,100,1,1), 26486.63085401582);
    AssertAlmostEqual(&test, Plm(5000,0,1000,1,1), 37999.68006129632);
    AssertAlmostEqual(&test, Plm(5000,0,1e+06,1,1), 72538.4577060812);
    AssertAlmostEqual(&test, Plm(5000,1,1.1,1,1), 2221.792745997359);
    AssertAlmostEqual(&test, Plm(5000,1,2,1,1), 6588.512931959001);
    AssertAlmostEqual(&test, Plm(5000,1,5,1,1), 11465.849681698872);
    AssertAlmostEqual(&test, Plm(5000,1,10,1,1), 14969.801694777063);
    AssertAlmostEqual(&test, Plm(5000,1,100,1,1), 26495.148047202234);
    AssertAlmostEqual(&test, Plm(5000,1,1000,1,1), 38008.197254487684);
    AssertAlmostEqual(&test, Plm(5000,1,1e+06,1,1), 72546.97489927262);
    AssertAlmostEqual(&test, Plm(5000,5,1.1,1,1), 2255.8561564051015);
    AssertAlmostEqual(&test, Plm(5000,5,2,1,1), 6622.579332800245);
    AssertAlmostEqual(&test, Plm(5000,5,5,1,1), 11499.916404369476);
    AssertAlmostEqual(&test, Plm(5000,5,10,1,1), 15003.868454850488);
    AssertAlmostEqual(&test, Plm(5000,5,100,1,1), 26529.214819247612);
    AssertAlmostEqual(&test, Plm(5000,5,1000,1,1), 38042.26402665189);
    AssertAlmostEqual(&test, Plm(5000,5,1e+06,1,1), 72581.04167143803);
    AssertAlmostEqual(&test, Plm(5000,10,1.1,1,1), 2298.4246117857556);
    AssertAlmostEqual(&test, Plm(5000,10,2,1,1), 6665.15713326439);
    AssertAlmostEqual(&test, Plm(5000,10,5,1,1), 11542.49521055019);
    AssertAlmostEqual(&test, Plm(5000,10,10,1,1), 15046.44737791502);
    AssertAlmostEqual(&test, Plm(5000,10,100,1,1), 26571.79377972449);
    AssertAlmostEqual(&test, Plm(5000,10,1000,1,1), 38084.84298750007);
    AssertAlmostEqual(&test, Plm(5000,10,1e+06,1,1), 72623.62063228997);
    AssertAlmostEqual(&test, Plm(5000,50,1.1,1,1), 2638.539374646917);
    AssertAlmostEqual(&test, Plm(5000,50,2,1,1), 7005.570921840936);
    AssertAlmostEqual(&test, Plm(5000,50,5,1,1), 11882.941181909891);
    AssertAlmostEqual(&test, Plm(5000,50,10,1,1), 15386.89708955445);
    AssertAlmostEqual(&test, Plm(5000,50,100,1,1), 26912.244688558854);
    AssertAlmostEqual(&test, Plm(5000,50,1000,1,1), 38425.293908216525);
    AssertAlmostEqual(&test, Plm(5000,50,1e+06,1,1), 72964.07155312643);
    AssertAlmostEqual(&test, Plm(5000,100,1.1,1,1), 3062.597933636185);
    AssertAlmostEqual(&test, Plm(5000,100,2,1,1), 7430.563724364048);
    AssertAlmostEqual(&test, Plm(5000,100,5,1,1), 12308.034553792588);
    AssertAlmostEqual(&test, Plm(5000,100,10,1,1), 15812.002149781383);
    AssertAlmostEqual(&test, Plm(5000,100,100,1,1), 27337.35349001806);
    AssertAlmostEqual(&test, Plm(5000,100,1000,1,1), 38850.40274680726);
    AssertAlmostEqual(&test, Plm(5000,100,1e+06,1,1), 73389.1803920922);
    AssertAlmostEqual(&test, Plm(5000,500,1.1,1,1), 6411.206340644299);
    AssertAlmostEqual(&test, Plm(5000,500,2,1,1), 10808.900960793368);
    AssertAlmostEqual(&test, Plm(5000,500,5,1,1), 15689.588541138106);
    AssertAlmostEqual(&test, Plm(5000,500,10,1,1), 19193.93014021267);
    AssertAlmostEqual(&test, Plm(5000,500,100,1,1), 30719.401198371503);
    AssertAlmostEqual(&test, Plm(5000,500,1000,1,1), 42232.45164336937);
    AssertAlmostEqual(&test, Plm(5000,500,1e+06,1,1), 76771.22930065551);
    AssertAlmostEqual(&test, Plm(5000,1000,1.1,1,1), 10485.8028979215);
    AssertAlmostEqual(&test, Plm(5000,1000,2,1,1), 14974.427949133653);
    AssertAlmostEqual(&test, Plm(5000,1000,5,1,1), 19865.1495998723);
    AssertAlmostEqual(&test, Plm(5000,1000,10,1,1), 23370.65965964008);
    AssertAlmostEqual(&test, Plm(5000,1000,100,1,1), 34896.504817425244);
    AssertAlmostEqual(&test, Plm(5000,1000,1000,1,1), 46409.558975573345);
    AssertAlmostEqual(&test, Plm(5000,1000,1e+06,1,1), 80948.33667036322);
    AssertAlmostEqual(&test, Plm(10000,0,1.1,1,1), 4430.77038620328);
    AssertAlmostEqual(&test, Plm(10000,0,2,1,1), 13164.438675839137);
    AssertAlmostEqual(&test, Plm(10000,0,5,1,1), 22919.14427702755);
    AssertAlmostEqual(&test, Plm(10000,0,10,1,1), 29927.05217156812);
    AssertAlmostEqual(&test, Plm(10000,0,100,1,1), 52977.74612097734);
    AssertAlmostEqual(&test, Plm(10000,0,1000,1,1), 76003.84454791597);
    AssertAlmostEqual(&test, Plm(10000,0,1e+06,1,1), 145081.39983761078);
    AssertAlmostEqual(&test, Plm(10000,1,1.1,1,1), 4439.980656547012);
    AssertAlmostEqual(&test, Plm(10000,1,2,1,1), 13173.64900847564);
    AssertAlmostEqual(&test, Plm(10000,1,5,1,1), 22928.354616368437);
    AssertAlmostEqual(&test, Plm(10000,1,10,1,1), 29936.262511688194);
    AssertAlmostEqual(&test, Plm(10000,1,100,1,1), 52986.95646134681);
    AssertAlmostEqual(&test, Plm(10000,1,1000,1,1), 76013.05488828792);
    AssertAlmostEqual(&test, Plm(10000,1,1e+06,1,1), 145090.61017798274);
    AssertAlmostEqual(&test, Plm(10000,5,1.1,1,1), 4476.819337207217);
    AssertAlmostEqual(&test, Plm(10000,5,2,1,1), 13210.489184162148);
    AssertAlmostEqual(&test, Plm(10000,5,5,1,1), 22965.194952960173);
    AssertAlmostEqual(&test, Plm(10000,5,10,1,1), 29973.10286698038);
    AssertAlmostEqual(&test, Plm(10000,5,100,1,1), 53023.796822624674);
    AssertAlmostEqual(&test, Plm(10000,5,1000,1,1), 76049.89524962519);
    AssertAlmostEqual(&test, Plm(10000,5,1e+06,1,1), 145127.45053932062);
    AssertAlmostEqual(&test, Plm(10000,10,1.1,1,1), 4522.862285676423);
    AssertAlmostEqual(&test, Plm(10000,10,2,1,1), 13256.536804585901);
    AssertAlmostEqual(&test, Plm(10000,10,5,1,1), 23011.243076212744);
    AssertAlmostEqual(&test, Plm(10000,10,10,1,1), 30019.151048671865);
    AssertAlmostEqual(&test, Plm(10000,10,100,1,1), 53069.84502302139);
    AssertAlmostEqual(&test, Plm(10000,10,1000,1,1), 76095.94345020756);
    AssertAlmostEqual(&test, Plm(10000,10,1e+06,1,1), 145173.49873990484);
    AssertAlmostEqual(&test, Plm(10000,50,1.1,1,1), 4890.989633810495);
    AssertAlmostEqual(&test, Plm(10000,50,2,1,1), 13624.813653147077);
    AssertAlmostEqual(&test, Plm(10000,50,5,1,1), 23379.536015277678);
    AssertAlmostEqual(&test, Plm(10000,50,10,1,1), 30387.445857781695);
    AssertAlmostEqual(&test, Plm(10000,50,100,1,1), 53438.140430698666);
    AssertAlmostEqual(&test, Plm(10000,50,1000,1,1), 76464.23886382558);
    AssertAlmostEqual(&test, Plm(10000,50,1e+06,1,1), 145541.7941535829);
    AssertAlmostEqual(&test, Plm(10000,100,1.1,1,1), 5350.607526954631);
    AssertAlmostEqual(&test, Plm(10000,100,2,1,1), 14084.89870864949);
    AssertAlmostEqual(&test, Plm(10000,100,5,1,1), 23839.67135337468);
    AssertAlmostEqual(&test, Plm(10000,100,10,1,1), 30847.587039765254);
    AssertAlmostEqual(&test, Plm(10000,100,100,1,1), 53898.28348320525);
    AssertAlmostEqual(&test, Plm(10000,100,1000,1,1), 76924.381934897);
    AssertAlmostEqual(&test, Plm(10000,100,1e+06,1,1), 146001.9372248418);
    AssertAlmostEqual(&test, Plm(10000,500,1.1,1,1), 9005.767857665567);
    AssertAlmostEqual(&test, Plm(10000,500,2,1,1), 17754.987132949453);
    AssertAlmostEqual(&test, Plm(10000,500,5,1,1), 27511.368637020114);
    AssertAlmostEqual(&test, Plm(10000,500,10,1,1), 34519.471324789796);
    AssertAlmostEqual(&test, Plm(10000,500,100,1,1), 57570.22762477764);
    AssertAlmostEqual(&test, Plm(10000,500,1000,1,1), 80596.32667054408);
    AssertAlmostEqual(&test, Plm(10000,500,1e+06,1,1), 149673.8819664892);
    AssertAlmostEqual(&test, Plm(10000,1000,1.1,1,1), 13519.734372125722);
    AssertAlmostEqual(&test, Plm(10000,1000,2,1,1), 22315.34403779489);
    AssertAlmostEqual(&test, Plm(10000,1000,5,1,1), 32076.750933853808);
    AssertAlmostEqual(&test, Plm(10000,1000,10,1,1), 39085.437963553515);
    AssertAlmostEqual(&test, Plm(10000,1000,100,1,1), 62136.38131289363);
    AssertAlmostEqual(&test, Plm(10000,1000,1000,1,1), 85162.48221514323);
    AssertAlmostEqual(&test, Plm(10000,1000,1e+06,1,1), 154240.03752983926);
    AssertAlmostEqual(&test, Plm(15000,0,1.1,1,1), 6648.408923905316);
    AssertAlmostEqual(&test, Plm(15000,0,2,1,1), 19749.025431431186);
    AssertAlmostEqual(&test, Plm(15000,0,5,1,1), 34381.09989636012);
    AssertAlmostEqual(&test, Plm(15000,0,10,1,1), 44892.963673791644);
    AssertAlmostEqual(&test, Plm(15000,0,100,1,1), 79469.00522064215);
    AssertAlmostEqual(&test, Plm(15000,0,1000,1,1), 114008.15286723853);
    AssertAlmostEqual(&test, Plm(15000,0,1e+06,1,1), 217624.48580184323);
    AssertAlmostEqual(&test, Plm(15000,1,1.1,1,1), 6658.024682701772);
    AssertAlmostEqual(&test, Plm(15000,1,2,1,1), 19758.641231754387);
    // AssertAlmostEqual(&test, Plm(15000,1,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,1,10,1,1), 44902.57947910379);
    AssertAlmostEqual(&test, Plm(15000,1,100,1,1), 79478.62102612057);
    AssertAlmostEqual(&test, Plm(15000,1,1000,1,1), 114017.7686727186);
    AssertAlmostEqual(&test, Plm(15000,1,1e+06,1,1), 217634.1016073233);
    AssertAlmostEqual(&test, Plm(15000,5,1.1,1,1), 6696.486117481743);
    AssertAlmostEqual(&test, Plm(15000,5,2,1,1), 19797.103663176185);
    // AssertAlmostEqual(&test, Plm(15000,5,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,5,10,1,1), 44941.0420302604);
    AssertAlmostEqual(&test, Plm(15000,5,100,1,1), 79517.08358126757);
    AssertAlmostEqual(&test, Plm(15000,5,1000,1,1), 114056.23122790518);
    AssertAlmostEqual(&test, Plm(15000,5,1e+06,1,1), 217672.5641625103);
    AssertAlmostEqual(&test, Plm(15000,10,1.1,1,1), 6744.559309710822);
    AssertAlmostEqual(&test, Plm(15000,10,2,1,1), 19845.179969910194);
    // AssertAlmostEqual(&test, Plm(15000,10,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,10,10,1,1), 44989.11871116567);
    AssertAlmostEqual(&test, Plm(15000,10,100,1,1), 79565.16027464278);
    AssertAlmostEqual(&test, Plm(15000,10,1000,1,1), 114104.30792140416);
    AssertAlmostEqual(&test, Plm(15000,10,1e+06,1,1), 217720.64085601055);
    AssertAlmostEqual(&test, Plm(15000,50,1.1,1,1), 7129.000732855435);
    AssertAlmostEqual(&test, Plm(15000,50,2,1,1), 20229.721056584967);
    // AssertAlmostEqual(&test, Plm(15000,50,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,50,10,1,1), 45373.671771315356);
    AssertAlmostEqual(&test, Plm(15000,50,100,1,1), 79949.71373383075);
    AssertAlmostEqual(&test, Plm(15000,50,1000,1,1), 114488.86138455257);
    AssertAlmostEqual(&test, Plm(15000,50,1e+06,1,1), 218105.19431919893);
    AssertAlmostEqual(&test, Plm(15000,100,1.1,1,1), 7609.191914188202);
    AssertAlmostEqual(&test, Plm(15000,100,2,1,1), 20710.22367860547);
    // AssertAlmostEqual(&test, Plm(15000,100,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,100,10,1,1), 45854.21181037582);
    AssertAlmostEqual(&test, Plm(15000,100,100,1,1), 80430.25501988578);
    AssertAlmostEqual(&test, Plm(15000,100,1000,1,1), 114969.40268298394);
    AssertAlmostEqual(&test, Plm(15000,100,1e+06,1,1), 218585.73561775533);
    AssertAlmostEqual(&test, Plm(15000,500,1.1,1,1), 11436.23685715974);
    AssertAlmostEqual(&test, Plm(15000,500,2,1,1), 24547.228458838814);
    // AssertAlmostEqual(&test, Plm(15000,500,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,500,10,1,1), 49692.413880574);
    AssertAlmostEqual(&test, Plm(15000,500,100,1,1), 84268.49699385418);
    AssertAlmostEqual(&test, Plm(15000,500,1000,1,1), 118807.64505299553);
    AssertAlmostEqual(&test, Plm(15000,500,1e+06,1,1), 222423.97799176705);
    AssertAlmostEqual(&test, Plm(15000,1000,1.1,1,1), 16183.571356075016);
    AssertAlmostEqual(&test, Plm(15000,1000,2,1,1), 29325.609657082314);
    // AssertAlmostEqual(&test, Plm(15000,1000,5,1,1), ???);
    AssertAlmostEqual(&test, Plm(15000,1000,10,1,1), 54474.535918872665);
    AssertAlmostEqual(&test, Plm(15000,1000,100,1,1), 89050.7437307356);
    AssertAlmostEqual(&test, Plm(15000,1000,1000,1,1), 123589.89302751188);
    AssertAlmostEqual(&test, Plm(15000,1000,1e+06,1,1), 227206.2259787838);

    return test_results(&test, stderr);
}

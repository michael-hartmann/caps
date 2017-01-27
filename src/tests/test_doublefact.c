#include "sfunc.h"
#include "unittest.h"

#include "test_doublefact.h"

int test_doublefact()
{
    unittest_t test;

    unittest_init(&test, "doublefact", "Test double factorial");

    AssertAlmostEqual(&test, ln_factorial2(0), 0);
    AssertAlmostEqual(&test, ln_factorial2(1), 0);
    AssertAlmostEqual(&test, ln_factorial2(2), 0.6931471805599453);
    AssertAlmostEqual(&test, ln_factorial2(3), 1.09861228866811);
    AssertAlmostEqual(&test, ln_factorial2(4), 2.079441541679836);
    AssertAlmostEqual(&test, ln_factorial2(5), 2.70805020110221);
    AssertAlmostEqual(&test, ln_factorial2(6), 3.871201010907891);
    AssertAlmostEqual(&test, ln_factorial2(7), 4.653960350157523);
    AssertAlmostEqual(&test, ln_factorial2(8), 5.950642552587727);
    AssertAlmostEqual(&test, ln_factorial2(9), 6.851184927493743);
    AssertAlmostEqual(&test, ln_factorial2(10), 8.253227645581772);
    AssertAlmostEqual(&test, ln_factorial2(11), 9.249080200292113);
    AssertAlmostEqual(&test, ln_factorial2(12), 10.73813429536977);
    AssertAlmostEqual(&test, ln_factorial2(13), 11.81402955775365);
    AssertAlmostEqual(&test, ln_factorial2(14), 13.37719162498503);
    AssertAlmostEqual(&test, ln_factorial2(15), 14.52207975885586);
    AssertAlmostEqual(&test, ln_factorial2(16), 16.14978034722481);
    AssertAlmostEqual(&test, ln_factorial2(17), 17.35529310291207);
    AssertAlmostEqual(&test, ln_factorial2(18), 19.04015210512098);
    AssertAlmostEqual(&test, ln_factorial2(19), 20.29973208207852);
    AssertAlmostEqual(&test, ln_factorial2(20), 22.03588437867497);
    AssertAlmostEqual(&test, ln_factorial2(21), 23.34425451980194);
    AssertAlmostEqual(&test, ln_factorial2(22), 25.12692683203328);
    AssertAlmostEqual(&test, ln_factorial2(23), 26.47974873573109);
    AssertAlmostEqual(&test, ln_factorial2(24), 28.30498066238123);
    AssertAlmostEqual(&test, ln_factorial2(25), 29.69862456059929);
    AssertAlmostEqual(&test, ln_factorial2(26), 31.56307720040271);
    AssertAlmostEqual(&test, ln_factorial2(27), 32.99446142660362);
    AssertAlmostEqual(&test, ln_factorial2(28), 34.89528171057792);
    AssertAlmostEqual(&test, ln_factorial2(29), 36.3617572565901);
    AssertAlmostEqual(&test, ln_factorial2(30), 38.29647909224007);
    AssertAlmostEqual(&test, ln_factorial2(31), 39.79574446107524);
    AssertAlmostEqual(&test, ln_factorial2(32), 41.7622149950398);
    AssertAlmostEqual(&test, ln_factorial2(33), 43.29225202254172);
    AssertAlmostEqual(&test, ln_factorial2(34), 45.28857551965596);
    AssertAlmostEqual(&test, ln_factorial2(35), 46.84760008403114);
    AssertAlmostEqual(&test, ln_factorial2(36), 48.87209445811207);
    AssertAlmostEqual(&test, ln_factorial2(37), 50.45851799667535);
    AssertAlmostEqual(&test, ln_factorial2(38), 52.50968061783846);
    AssertAlmostEqual(&test, ln_factorial2(39), 54.12207964280501);
    AssertAlmostEqual(&test, ln_factorial2(40), 56.19856007195239);
    AssertAlmostEqual(&test, ln_factorial2(41), 57.83565170950931);
    AssertAlmostEqual(&test, ln_factorial2(42), 59.93622969023576);
    AssertAlmostEqual(&test, ln_factorial2(43), 61.59685182520288);
    AssertAlmostEqual(&test, ln_factorial2(44), 63.72041932415402);
    AssertAlmostEqual(&test, ln_factorial2(45), 65.40351431497319);
    AssertAlmostEqual(&test, ln_factorial2(46), 67.54906072064311);
    AssertAlmostEqual(&test, ln_factorial2(47), 69.25366191668326);
    AssertAlmostEqual(&test, ln_factorial2(48), 71.42026173155101);
    AssertAlmostEqual(&test, ln_factorial2(49), 73.14548221479387);
    AssertAlmostEqual(&test, ln_factorial2(50), 75.33228473697915);
    AssertAlmostEqual(&test, ln_factorial2(51), 77.0773078475182);
    AssertAlmostEqual(&test, ln_factorial2(52), 79.28352845556059);
    AssertAlmostEqual(&test, ln_factorial2(53), 81.04759976107033);
    AssertAlmostEqual(&test, ln_factorial2(54), 83.27251250212485);
    AssertAlmostEqual(&test, ln_factorial2(55), 85.0549329463028);
    AssertAlmostEqual(&test, ln_factorial2(56), 87.29786419286);
    AssertAlmostEqual(&test, ln_factorial2(57), 89.09798421413734);
    AssertAlmostEqual(&test, ln_factorial2(58), 91.35830720340643);
    AssertAlmostEqual(&test, ln_factorial2(59), 93.17552165804307);
    AssertAlmostEqual(&test, ln_factorial2(60), 95.45265176562853);
    AssertAlmostEqual(&test, ln_factorial2(61), 97.28639552221638);
    AssertAlmostEqual(&test, ln_factorial2(62), 99.57978615067361);
    AssertAlmostEqual(&test, ln_factorial2(63), 101.4295302486079);
    AssertAlmostEqual(&test, ln_factorial2(64), 103.7386692340333);
    AssertAlmostEqual(&test, ln_factorial2(65), 105.6039175185036);
    AssertAlmostEqual(&test, ln_factorial2(66), 107.9283239760597);
    AssertAlmostEqual(&test, ln_factorial2(67), 109.8086101378945);
    AssertAlmostEqual(&test, ln_factorial2(68), 112.1478316812358);
    AssertAlmostEqual(&test, ln_factorial2(69), 114.0427166424918);
    AssertAlmostEqual(&test, ln_factorial2(70), 116.3963269232852);
    AssertAlmostEqual(&test, ln_factorial2(71), 118.3053965195331);
    AssertAlmostEqual(&test, ln_factorial2(72), 120.6729930423012);
    AssertAlmostEqual(&test, ln_factorial2(73), 122.5958559606815);
    AssertAlmostEqual(&test, ln_factorial2(74), 124.9770581355054);
    AssertAlmostEqual(&test, ln_factorial2(75), 126.9133440742178);
    AssertAlmostEqual(&test, ln_factorial2(76), 129.3077914757917);
    AssertAlmostEqual(&test, ln_factorial2(77), 131.2571494960715);
    AssertAlmostEqual(&test, ln_factorial2(78), 133.6645003024813);
    AssertAlmostEqual(&test, ln_factorial2(79), 135.6265973485385);
    AssertAlmostEqual(&test, ln_factorial2(80), 138.0465269371552);
    AssertAlmostEqual(&test, ln_factorial2(81), 140.0210465032109);
    AssertAlmostEqual(&test, ln_factorial2(82), 142.4532461844195);
    AssertAlmostEqual(&test, ln_factorial2(83), 144.4398871110075);
    AssertAlmostEqual(&test, ln_factorial2(84), 146.8840629832628);
    AssertAlmostEqual(&test, ln_factorial2(85), 148.8825383674979);
    AssertAlmostEqual(&test, ln_factorial2(86), 151.3384102795163);
    AssertAlmostEqual(&test, ln_factorial2(87), 153.3484464861524);
    AssertAlmostEqual(&test, ln_factorial2(88), 155.8157470939945);
    AssertAlmostEqual(&test, ln_factorial2(89), 157.8370828558846);
    AssertAlmostEqual(&test, ln_factorial2(90), 160.3155567643248);
    AssertAlmostEqual(&test, ln_factorial2(91), 162.3479423624014);
    AssertAlmostEqual(&test, ln_factorial2(92), 164.8373453413738);
    AssertAlmostEqual(&test, ln_factorial2(93), 166.8805418555547);
    AssertAlmostEqual(&test, ln_factorial2(94), 169.3806401236438);
    AssertAlmostEqual(&test, ln_factorial2(95), 171.4344187471552);
    AssertAlmostEqual(&test, ln_factorial2(96), 173.9449883151116);
    AssertAlmostEqual(&test, ln_factorial2(97), 176.0091297256586);
    AssertAlmostEqual(&test, ln_factorial2(98), 178.5299557937822);
    AssertAlmostEqual(&test, ln_factorial2(99), 180.6042495757932);
    AssertAlmostEqual(&test, ln_factorial2(100), 183.1351259797703);
    
    AssertAlmostEqual(&test, ln_factorial2(333), 803.80686991691279488);
    AssertAlmostEqual(&test, ln_factorial2(334), 806.93898026792161962);

    AssertAlmostEqual(&test, ln_factorial2(499), 1303.9984315293167964);
    AssertAlmostEqual(&test, ln_factorial2(500), 1307.3320269308392879);

    AssertAlmostEqual(&test, ln_factorial2(1000), 2957.904048740129);
    AssertAlmostEqual(&test, ln_factorial2(2000), 6605.275359048108);
    AssertAlmostEqual(&test, ln_factorial2(2000), 6605.275359048108);
    AssertAlmostEqual(&test, ln_factorial2(3000), 10514.12695575767);
    AssertAlmostEqual(&test, ln_factorial2(4000), 14592.8187116337);
    AssertAlmostEqual(&test, ln_factorial2(5000), 18797.81397341256);
    AssertAlmostEqual(&test, ln_factorial2(6000), 23103.46639472538);
    AssertAlmostEqual(&test, ln_factorial2(7000), 27492.82821959754);
    AssertAlmostEqual(&test, ln_factorial2(8000), 31953.85326683448);
    AssertAlmostEqual(&test, ln_factorial2(9000), 36477.5342268222);
    AssertAlmostEqual(&test, ln_factorial2(10000), 41056.87941167649);
    AssertAlmostEqual(&test, ln_factorial2(20000), 89040.3996424138);
    AssertAlmostEqual(&test, ln_factorial2(30000), 139640.0167564932);
    AssertAlmostEqual(&test, ln_factorial2(40000), 191938.5653483976);
    AssertAlmostEqual(&test, ln_factorial2(50000), 245500.4393676755);

    return test_results(&test, stderr);
}

#include "sfunc.h"
#include "unittest.h"

#include "test_doublefact.h"

int test_doublefact()
{
    unittest_t test;

    unittest_init(&test, "doublefact", "Test double factorial");

    AssertEqual(&test, ln_factorial2(0), 0);
    AssertEqual(&test, ln_factorial2(1), 0);

    AssertAlmostEqual(&test, ln_factorial2(2),   0.6931471805599453094);
    AssertAlmostEqual(&test, ln_factorial2(3),   1.0986122886681096913);
    AssertAlmostEqual(&test, ln_factorial2(4),   2.0794415416798359282);
    AssertAlmostEqual(&test, ln_factorial2(5),   2.7080502011022100659);

    AssertAlmostEqual(&test, ln_factorial2(51),  77.077307847518205164);
    AssertAlmostEqual(&test, ln_factorial2(52),  79.283528455560580029);
    
    AssertAlmostEqual(&test, ln_factorial2(100), 183.13512597977029753);
    AssertAlmostEqual(&test, ln_factorial2(101), 185.21937009263445205);
    
    AssertAlmostEqual(&test, ln_factorial2(333), 803.80686991691279488);
    AssertAlmostEqual(&test, ln_factorial2(334), 806.93898026792161962);

    AssertAlmostEqual(&test, ln_factorial2(499), 1303.9984315293167964);
    AssertAlmostEqual(&test, ln_factorial2(500), 1307.3320269308392879);

    return test_results(&test, stderr);
}

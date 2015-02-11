#include "sfunc.h"
#include "unittest.h"

#include "test_doublefact.h"

int test_doublefact()
{
    unittest_t test;

    unittest_init(&test, "doublefact", "Test double factorial");

    AssertEqual(&test, ln_doublefact(0),   0);
    AssertEqual(&test, ln_doublefact(1),   0);

    AssertAlmostEqual(&test, ln_doublefact(2),   0.693147180559945309417232121458176568075500134360255254120680);
    AssertAlmostEqual(&test, ln_doublefact(3),   1.098612288668109691395245236922525704647490557822749451734694);
    AssertAlmostEqual(&test, ln_doublefact(4),   2.079441541679835928251696364374529704226500403080765762362040);
    AssertAlmostEqual(&test, ln_doublefact(5),   2.708050201102210065996004570148713344173091912091267173647342);

    AssertAlmostEqual(&test, ln_doublefact(51),  77.07730784751820516445630791916220149454265081179451116106491);
    AssertAlmostEqual(&test, ln_doublefact(52),  79.28352845556058002961361747099464590876118666235724979722338);
    
    AssertAlmostEqual(&test, ln_doublefact(100), 183.1351259797702975383987999237883518258868824082753118006299);
    AssertAlmostEqual(&test, ln_doublefact(101), 185.2193700926344520565653917127802753588279212418095676769807);
    
    AssertAlmostEqual(&test, ln_doublefact(333), 803.8068699169127948868947011937225159269950980762322648014825);
    AssertAlmostEqual(&test, ln_doublefact(334), 806.9389802679216196222732317043854384535329497966731766081127);

    AssertAlmostEqual(&test, ln_doublefact(499), 1303.998431529316796415491154267014961848434758393410615755047);
    AssertAlmostEqual(&test, ln_doublefact(500), 1307.332026930839287985819861539370764480277316393340117718467);

    return test_results(&test, stderr);
}

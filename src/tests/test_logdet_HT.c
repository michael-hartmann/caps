#include "libcasimir.h"
#include "unittest.h"

#include "test_logdet_HT.h"

double data[][5] =
{
    /* L/R, m, lmax, EE, MM */
    { 1, 0, 30, -0.039359509611808, -0.020914358073834 },
    { 1, 1, 30, -0.020914358073834, -0.011383372269677 },
    { 1, 2, 30, -0.0014890154006822, -0.0010439636682519 },
    { 1, 3, 30, -0.00010684253505352, -8.3068160092048e-05 },
    { 1, 4, 30, -7.6706194831114e-06, -6.3114778894253e-06 },
    { 1, 5, 30, -5.5072400278604e-07, -4.6967304805243e-07 },

    { 0.5, 0, 30, -0.11933466080803, -0.066894561496325 },
    { 0.5, 1, 30, -0.066894561496325, -0.039189251736552 },
    { 0.5, 2, 30, -0.0095534474374675, -0.0070091590402225 },
    { 0.5, 3, 30, -0.0013895949379127, -0.0011164478044288 },
    { 0.5, 4, 30, -0.00020264950726963, -0.00017108137032858 },
    { 0.5, 5, 30, -2.956425727734e-05, -2.5751203034666e-05 },

    { 0.1, 0, 70, -0.75828376979107, -0.49961137357808 },
    { 0.1, 1, 70, -0.49961137357808, -0.3573053616443 },
    { 0.1, 2, 70, -0.19269085420424, -0.15998670819642 },
    { 0.1, 3, 70, -0.077455287593816, -0.068141090725342 },
    { 0.1, 4, 70, -0.031594057315748, -0.028649692906629 },
    { 0.1, 5, 70, -0.012960911806234, -0.0119739786861 },

    { -1, -1, -1, -1, -1}
};

static void _logdet(double LbyR, int m, int lmax, double *EE, double *MM)
{
    casimir_t casimir;

    casimir_init(&casimir, LbyR, 1);
    casimir_set_lmax(&casimir, lmax);
    casimir_logdetD0(&casimir, m, EE, MM);
    casimir_free(&casimir);
}

int test_logdet_HT(void)
{
    unittest_t test;
    unittest_init(&test, "logdetD (HT)", "calculate logdetD for T→∞");

    for(int i = 0; ; i++)
    {
        double LbyR, EE_correct, MM_correct, EE_computed = 0, MM_computed = 0;
        int lmax,m;
        
        LbyR = data[i][0];
        m    = data[i][1];
        lmax = data[i][2];
        EE_correct = data[i][3];
        MM_correct = data[i][4];

        if(LbyR <= 0)
            break;

        _logdet(LbyR, m, lmax, &EE_computed, &MM_computed);

        AssertAlmostEqual(&test, EE_correct, EE_computed);
        AssertAlmostEqual(&test, MM_correct, MM_computed);
    }

    return test_results(&test, stderr);
}

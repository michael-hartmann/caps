#include "libcasimir.h"
#include "unittest.h"
#include "integration_drude.h"
#include "integration_perf.h"

#include "test_logdet.h"

static double data[][5] =
{
    /* L/R, nT, m, lmax, logdetD */
    {1, 0.1, 1, 50, -0.0326493072428 },
    {1, 0.5, 1, 50, -0.0345873881519 },
    {1, 1,   1, 50, -0.0306060278148 },
    {1, 2,   1, 50, -0.0158943867684 },
    {1, 5,   1, 50, -0.00102643843796 },
    {1, 10,  1, 50, -7.56809680963e-06 },

    {1, 0.1, 5, 50, -1.01943442088e-06 },
    {1, 0.5, 5, 50, -9.96556270262e-07 },
    {1, 1,   5, 50, -9.2800603016e-07 },
    {1, 2,   5, 50, -6.99048536402e-07 },
    {1, 5,   5, 50, -1.26756284443e-07 },
    {1, 10,  5, 50, -2.04558168797e-09 },

    {0.1, 0.1, 1, 75, -0.859291428533 },
    {0.1, 0.2, 1, 75, -0.864187288762 },
    {0.1, 0.5, 1, 75, -0.875454187407 },
    {0.1, 1,   1, 75, -0.857971481783 },
    {0.1, 2,   1, 75, -0.749247839073 },
    {0.1, 5,   1, 75, -0.428752561189 },
    {0.1, 10,  1, 75, -0.165907753049 },

    {0.1, 0.1, 5, 75, -0.0249311229898 },
    {0.1, 0.2, 5, 75, -0.024919229991 },
    {0.1, 0.5, 5, 75, -0.0248362338789 },
    {0.1, 1,   5, 75, -0.0245435935019 },
    {0.1, 2,   5, 75, -0.0234342078822 },
    {0.1, 5,   5, 75, -0.017746834386 },
    {0.1, 10,  5, 75, -0.00884026743848 },

    {0.02, 0.1, 1, 300, -3.78882984277 },
    {0.02, 0.2, 1, 300, -3.79591725861 },
    {0.02, 0.5, 1, 300, -3.81575540076 },
    {0.02, 1,   1, 300, -3.80753561696 },
    {0.02, 2,   1, 300, -3.68632708782 },
    {0.02, 5,   1, 300, -3.17058502525 },
    {0.02, 10,  1, 300, -2.44322158172 },

    {0.02, 0.1, 5, 300, -0.685917710837 },
    {0.02, 0.2, 5, 300, -0.68585848661 },
    {0.02, 0.5, 5, 300, -0.685444536019 },
    {0.02, 1,   5, 300, -0.68397532647 },
    {0.02, 2,   5, 300, -0.678249499702 },
    {0.02, 5,   5, 300, -0.643790097746 },
    {0.02, 10,  5, 300, -0.560158367159 },

    {0.006, 1.006, 1, 500, -9.27424711542347 },

    { -1, -1, -1, -1, -1 }
};

static double _logdet(double LbyR, double nT, int m, int lmax);

static double _logdet(double LbyR, double nT, int m, int lmax)
{
    casimir_t casimir;
    double logdet;

    casimir_init(&casimir, LbyR, nT);
    casimir_set_lmax(&casimir, lmax);
    logdet = casimir_logdetD(&casimir, 1, m);
    casimir_free(&casimir);

    return logdet;
}

int test_logdet(void)
{
    unittest_t test;

    unittest_init(&test, "logdet M", "calculate logdet");

    for(int i = 0; ; i++)
    {
        double LbyR = data[i][0];

        if(LbyR < 0)
            break;

        double nT = data[i][1];
        int m     = data[i][2];
        int lmax  = data[i][3];

        double logdetD_expected = data[i][4];
        double logdetD_computed = _logdet(LbyR, nT, m, lmax);

        AssertAlmostEqual(&test, logdetD_expected, logdetD_computed);
    }

    return test_results(&test, stderr);
}

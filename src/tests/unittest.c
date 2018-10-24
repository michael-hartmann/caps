#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "unittest.h"

void unittest_init(unittest_t *test, const char *func, const char *desc, double eps)
{
    test->passed = test->failed = 0;
    test->func   = func;
    test->desc   = desc;
    test->eps    = eps;
    test->start  = now();
}

int test_results(unittest_t *test, FILE *stream)
{
    double t = now()-test->start;

    fprintf(stream, "[%9d/%9d]\t%-16s\t%-32s", test->passed, test->passed+test->failed, test->func, test->desc);

    if(t < 0.1)
        fprintf(stream, " (%.3fms)\t", t*1000);
    else
        fprintf(stream, " (%.3fs)\t", t);

    if(test->failed == 0)
        fprintf(stream, " [" KGRN "OK" KNRM "]\n");
    else
        fprintf(stream, " [" KRED "FAIL" KNRM "]\n");

    return test->failed;
}

int _AssertEqual(int line, unittest_t *test, int x, int y)
{
    if(x == y)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, KRED "FAILED" KNRM ": %d != %d on line %d\n", x, y, line);
        return 1;
    }
}

int _Assert(int line, unittest_t *test, int boolean)
{
    if(boolean)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, KRED "FAILED" KNRM ": on line %d\n", line);
        return 1;
    }
}

int _AssertAlmostEqual(int line, unittest_t *test, double x, double y, double eps)
{
    if(x == y || fabs(1-x/y) < eps)
    {
        test->passed++;
        return 0;
    }
    else
    {
        test->failed++;
        fprintf(stderr, KRED "FAILED" KNRM ": %.20g != %.20g (%g) on line %d\n", x, y, fabs(1-x/y), line);
        return 1;
    }
}

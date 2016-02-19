/**
 * @file   utils.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2015
 * @brief  helper functions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include "utils.h"

void (*error_handler)(const char *) = default_error_handler;

void set_defaut_error_handler(void (*f)(const char *))
{
    error_handler = f;
}

void default_error_handler(const char *reason)
{
    fprintf(stderr, "%s\n", reason);
    exit(1);
}

void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if(p == NULL)
        error_handler("malloc failed.");

    return p;
}

void *xrealloc(void *p, size_t size)
{
    void *p2 = realloc(p, size);
    if(p2 == NULL)
    {
        free(p);
        error_handler("realloc failed.");
    }

    return p2;
}

void _xfree(void *p)
{
    if(p == NULL)
        error_handler("free on NULL.");

    free(p);
}

/* This function returns the seconds since 1st Jan 1970 in µs precision */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

/* This function counts how many times the character c in string str appears */
int cinstr(const char *str, char c)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c)
            i++;

    return i;
}

/* This function returns a pointer to the n-th occurrence of the character c in
 * the string s. If the character c occures less than n times, NULL is
 * returned. 
 */
const char *indexn(const char *str, char c, int n)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c && ++i == n)
            return str;

    return NULL;
}

void swap(double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

void sec2human(double t, int *h, int *m, int *s)
{
    int time = round(t);

    *s = time % 60;
    *m = (time % 3600)/60;
    *h = time/3600;
}

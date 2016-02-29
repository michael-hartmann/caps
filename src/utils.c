/**
 * @file   utils.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   February, 2016
 * @brief  helper functions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include "utils.h"

void (*error_handler)(const char *) = default_error_handler;

/** @brief Set error handler
 *
 * Set callback to a function f that is called when an error occured during
 * calls to malloc, realloc or free. The function is given a string with an
 * description of the error.
 *
 * The default error handler prints the message to stderr and exits
 *
 * @param f callback to error handler
 */
void set_defaut_error_handler(void (*f)(const char *))
{
    error_handler = f;
}


/** @brief Default error handler
 *
 * This function prints the string reason to stderr and calls exit(1).
 *
 * @param reason string with description of error
 */
void default_error_handler(const char *reason)
{
    fprintf(stderr, "%s\n", reason);
    exit(1);
}


/** @brief Wrapper for malloc
 *
 * This function is a wrapper for malloc. If an error occures, the error handler will be called.
 *
 * @param size size of bytes to allocate
 * @retval ptr pointer to memory
 */
void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if(p == NULL)
        error_handler("malloc failed.");

    return p;
}


/** @brief Wrapper for realloc
 *
 * This function is a wrapper for realloc. If an error occures, the error
 * handler will be called.
 *
 * @param oldptr ptr to old memory
 * @param size size
 * @retval newptr pointer to new memory
 */
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

/** @brief Wrapper for xfree
 *
 * This function is a wrapper for free. If an error occures, the error handler
 * will be called.
 *
 * You usually want to use the macro xfree that also sets the pointer p to
 * NULL.
 *
 * @param ptr ptr to free
 */
void _xfree(void *ptr)
{
    if(ptr == NULL)
        error_handler("free on NULL.");

    free(ptr);
}

/** @brief Get time
 *
 * This function returns the seconds since 1st Jan 1970 in µs precision.
 *
 * @retval time seconds sind 1st Jan 1970
 */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

/** @brief Find character in string
 *
 * This function counts how many times the character c in string str appears.
 *
 * @param str string
 * @param c character
 * retval how many times c is in str
 */
int cinstr(const char *str, char c)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c)
            i++;

    return i;
}

/** @brief Find n-th occurence of character in string
 *
 * This function returns a pointer to the n-th occurrence of the character c in
 * the string s. If the character c occures less than n times, NULL is
 * returned.
 *
 * @param str string
 * @param c character
 * @param n occurence
 * @retval NULL if c occures less than n times in str
 * @retval ptr pointer to n-th occurence of c
 */
const char *indexn(const char *str, char c, int n)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c && ++i == n)
            return str;

    return NULL;
}

/** @brief Swap values a and b
 *
 * After the call: a=b and b=a.
 *
 * @param a first double
 * @param b second double
 */
void swap(double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}


/** @brief Convert seconds to HH:MM:SS
 *
 * Convert seconds to hours, minutes and seconds.
 *
 * @param t seconds
 * @param h hours
 * @param m minutes
 * @param s seconds
 */
void sec2human(double t, int *h, int *m, int *s)
{
    int time = round(t);

    *s = time % 60;
    *m = (time % 3600)/60;
    *h = time/3600;
}

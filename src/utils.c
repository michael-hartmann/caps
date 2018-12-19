/**
 * @file   utils.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2017
 * @brief  wrappers for malloc, calloc realloc, and a few more useful functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "utils.h"

/** @brief Wrapper for malloc
 *
 * This function is a wrapper for malloc. If malloc fails \ref TERMINATE is
 * called.
 *
 * @param size size of bytes to allocate
 * @retval ptr pointer to memory
 */
void *xmalloc(size_t size)
{
    void *p = malloc(size);

    TERMINATE(p == NULL, "malloc failed: size=%zu", size);

    return p;
}


/** @brief Wrapper for calloc
 *
 * This function is a wrapper for calloc. If calloc fails \ref TERMINATE is
 * called.
 *
 * @param mnemb
 * @param size
 * @retval ptr pointer to memory
 */
void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);

    TERMINATE(p == NULL, "calloc failed: nmemb=%zu, size=%zu", nmemb, size);

    return p;
}


/** @brief Wrapper for realloc
 *
 * This function is a wrapper for realloc. If realloc fails \ref TERMINATE is
 * called.
 *
 * @param oldptr ptr to old memory
 * @param size size
 * @retval newptr pointer to new memory
 */
void *xrealloc(void *p, size_t size)
{
    void *p2 = realloc(p, size);

    TERMINATE(p2 == NULL, "realloc failed: p=%p, size=%zu", p, size);

    return p2;
}

/** @brief Seconds since 01/01/1970
 *
 * This function returns the seconds since 1st Jan 1970 in µs precision.
 *
 * @retval time seconds since 1st Jan 1970
 */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

/** @brief Write time into string
 *
 * Write current time in a human readable format into string s. The output is
 * similar to "Aug 30 2018 14:37:35".
 *
 * @param str array
 * @param len length of array
 */
void time_as_string(char *s, size_t len)
{
    time_t rawtime;
    struct tm *info;

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(s, len, "%c", info);
}

/** @brief Disable buffering to stderr and stdout
 *
 */
void disable_buffering(void)
{
    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
}

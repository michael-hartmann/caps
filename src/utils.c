/**
 * @file   utils.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   January, 2017
 * @brief  wrappers for malloc, calloc realloc, and a few more useful functions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include "utils.h"

/** @brief Wrapper for malloc
 *
 * This function is a wrapper for malloc. If malloc fails TERMINATE is called.
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
 * This function is a wrapper for calloc. If calloc fails TERMINATE is called.
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
 * This function is a wrapper for realloc. If realloc fails TERMINATE is
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

/** @brief Wrapper for xfree
 *
 * This function is a wrapper for free.
 *
 * You usually want to use the macro xfree that also sets the pointer p to
 * NULL.
 *
 * @param ptr ptr to free
 */
void _xfree(void *ptr)
{
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

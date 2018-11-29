/**
 * @file   utils.h
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   July, 2017
 * @brief  wrappers for malloc, calloc and realloc, assert-like macros, now()-function
 */

#ifndef __UTILS_H
#define __UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

/*! name of compile */
#if defined(__clang__)
    #define COMPILER "clang/llvm"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
    #define COMPILER "icc"
#elif defined(__GNUC__) || defined(__GNUG__)
    #define COMPILER "gcc"
#else
    #define COMPILER "unknown"
#endif

/*! macro similar to assert that prints a warning to stderr and aborts */
#define TERMINATE(cond, ...) if(cond) { fprintf(stderr, "Fatal error: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); abort(); }

/*! macro similar to assert that prints a warning to stderr */
#define WARN(cond, ...)      if(cond) { fprintf(stderr, "Warning: ");     fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); }

/*! macro for free that sets pointer p to NULL after freeing memory */
#define xfree(p) do { _xfree(p); p = NULL; } while(0)

double now(void);
void time_as_string(char *s, size_t len);

void *xmalloc(size_t size);
void *xrealloc(void *p, size_t size);
void *xcalloc(size_t nmemb, size_t size);
void _xfree(void *p);

void disable_buffering(void);

#ifdef __cplusplus
}
#endif

#endif

/**
 * @file   utils.h
 * @author Michael Hartmann <caps@speicherleck.de>
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
    #define COMPILER "clang"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
    #define COMPILER "icc"
#elif defined(__GNUC__) || defined(__GNUG__)
    #define COMPILER "gcc"
#elif defined(_MSC_VER)
    #define COMPILER "msvc"
#else
    #define COMPILER "unknown"
#endif

#ifndef TERMINATE
/*! Macro similar to assert that prints a warning to stderr and aborts */
#define TERMINATE(cond, ...) if(cond) { fprintf(stderr, "Fatal error: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); abort(); }
#endif

#ifndef WARN
/*! macro similar to assert that prints a warning to stderr */
#define WARN(cond, ...)      if(cond) { fprintf(stderr, "Warning: ");     fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); }
#endif

/*! macro for free that sets pointer p to NULL after freeing memory */
#define xfree(p) do { free(p); p = NULL; } while(0)

double now(void);
void time_as_string(char *s, size_t len);

void *xmalloc(size_t size);
void *xrealloc(void *p, size_t size);
void *xcalloc(size_t nmemb, size_t size);

void disable_buffering(void);

void strrep(char *s, const char a, const char b);
void strim(char *str);

#ifdef __cplusplus
}
#endif

#endif

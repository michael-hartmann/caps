#ifndef __UTILS_H
#define __UTILS_H

#include <stdlib.h>

#if defined(__clang__)
    #define COMPILER "clang/llvm"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
    #define COMPILER "icc"
#elif defined(__GNUC__) || defined(__GNUG__)
    #define COMPILER "gcc"
#else
    #define COMPILER "unknown"
#endif

#define TERMINATE(cond, ...) if(cond) { fprintf(stderr, "Fatal error: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); exit(1); }
#define WARN(cond, ...)      if(cond) { fprintf(stderr, "Warning: ");     fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); }


#define xfree(p) do { _xfree(p); p = NULL; } while (0)

double now(void);

void *xmalloc(size_t size);
void *xrealloc(void *p, size_t size);
void _xfree(void *p);

int cinstr(const char *str, char c);
const char *indexn(const char *str, char c, int n);

void swap(int *a, int *b);

#endif

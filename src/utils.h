#ifndef __UTILS_H
#define __UTILS_H

#ifndef NDEBUG
#define TERMINATE(cond, ...) if(cond) { fprintf(stderr, "fatal error: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); exit(1); }
#define WARN(cond, ...) if(cond) { fprintf(stderr, "warning: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, " (in %s, %s:%d)\n", __func__, __FILE__, __LINE__); }
#else
#define TERMINATE(cond, ...)
#define WARN(cond, ...)
#endif


#define xfree(p) do { _xfree(p); p = NULL; } while (0)

double now(void);

void set_defaut_error_handler(void (*f)(const char *));
void default_error_handler(const char *str);

void *xmalloc(size_t size);
void *xrealloc(void *p, size_t size);
void _xfree(void *p);

int cinstr(const char *str, char c);
const char *indexn(const char *str, char c, int n);

void swap(double *a, double *b);

#endif

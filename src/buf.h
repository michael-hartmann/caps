/**
 * @file buf.h
 * @author Christopher Wellons
 * @date   September, 2018
 * @brief growable memory buffers for C99
 *
 * This is free and unencumbered software released into the public domain.
 *
 * Original version from https://github.com/skeeto/growable-buf.
 *
 * Note: buf_push(), buf_grow(), buf_trunc(), and buf_free() may change
 * the buffer pointer, and any previously-taken pointers should be
 * considered invalidated.
 *
 * Example usage:
 *
 *     float *values = 0;
 *     for (size_t i = 0; i < 25; i++)
 *         buf_push(values, rand() / (float)RAND_MAX);
 *     for (size_t i = 0; i < buf_size(values); i++)
 *         printf("values[%zu] = %f\n", i, values[i]);
 *     buf_free(values);
 */
#include <stddef.h>
#include <stdlib.h>

#include "utils.h"

#ifndef BUF_INIT_CAPACITY
    #define BUF_INIT_CAPACITY 32
#endif

struct buf {
    size_t capacity; /**< total capacity of buffer */
    size_t size;     /**< size / number of elements */
    char buffer[];   /**< buffer */
};

#define buf_ptr(v) \
    ((struct buf *)((char *)(v) - offsetof(struct buf, buffer)))

/** destroy and free the buffer */
#define buf_free(v) \
    do { \
        if (v) { \
            free(buf_ptr((v))); \
            (v) = 0; \
        } \
    } while (0)

/** return the number of elements in the buffer (size_t) */
#define buf_size(v) \
    ((v) ? buf_ptr((v))->size : 0)

/** return the total capacity of the buffer (size_t) */
#define buf_capacity(v) \
    ((v) ? buf_ptr((v))->capacity : 0)

/** append an element E to the end */
#define buf_push(v, e) \
    do { \
        if (buf_capacity((v)) == buf_size((v))) { \
            (v) = buf_grow1(v, sizeof(*(v)), \
                            !buf_capacity((v)) ? \
                              BUF_INIT_CAPACITY : \
                              buf_capacity((v))); \
        } \
        (v)[buf_ptr((v))->size++] = (e); \
    } while (0)

/** remove and return an element E from the end */
#define buf_pop(v) \
    ((v)[--buf_ptr(v)->size])

/** increase buffer capactity by (ptrdiff_t) N elements */
#define buf_grow(v, n) \
    ((v) = buf_grow1((v), sizeof(*(v)), n))

/** set buffer capactity to exactly (ptrdiff_t) N elements */
#define buf_trunc(v, n) \
    ((v) = buf_grow1((v), sizeof(*(v)), n - buf_capacity(v)))

/** set buffer size to 0 (for push/pop) */
#define buf_clear(v) \
    ((v) ? (buf_ptr((v))->size = 0) : 0)


static void *
buf_grow1(void *v, size_t esize, ptrdiff_t n)
{
    struct buf *p;
    size_t max = (size_t)-1 - sizeof(struct buf);
    if (v) {
        p = buf_ptr(v);
        if (n > 0 && p->capacity + n > max / esize)
            goto fail; /* overflow */
        p = realloc(p, sizeof(struct buf) + esize * (p->capacity + n));
        if (!p)
            goto fail;
        p->capacity += n;
        if (p->size > p->capacity)
            p->size = p->capacity;
    } else {
        if ((size_t)n > max / esize)
            goto fail; /* overflow */
        p = malloc(sizeof(struct buf) + esize * n);
        if (!p)
            goto fail;
        p->capacity = n;
        p->size = 0;
    }
    return p->buffer;
fail:
    TERMINATE(true, "Could not allocate memory.");
    return 0;
}

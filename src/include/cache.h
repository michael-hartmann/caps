#ifndef __CACHE_H
#define __CACHE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

typedef struct {
    uint64_t key;
    double value;
} cache_entry_t;

typedef struct {
    unsigned int num_entries;
    cache_entry_t *table;
} cache_t;

cache_t *cache_new(unsigned int entries);
void cache_free(cache_t *cache);
void cache_insert(cache_t *cache, uint64_t key, double value);
double cache_lookup(cache_t *cache, uint64_t key);

#ifdef __cplusplus
}
#endif

#endif

#ifndef __CACHE_H
#define __CACHE_H

#define CACHE_NUM_ENTRIES 262144
#define CACHE_NUM_LOOKUP 50331653

#include <stdint.h>

typedef struct {
    int head, tail, num_entries, num_lookup;
    uint64_t *keys;
    double *values;
    uint64_t *table;
} cache_t;

cache_t *cache_new(int entries, double filling);
void cache_free(cache_t *cache);
void cache_insert(cache_t *cache, uint64_t key, double value);
double cache_lookup(cache_t *cache, uint64_t key);

#endif

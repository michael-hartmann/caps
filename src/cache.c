#include <stdint.h>
#include <math.h>

#include "utils.h"
#include "cache.h"


cache_t *cache_new(int entries, double filling)
{
    const int primes[] = {
        1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433,
        1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319,
        201326611, 402653189, 805306457, 1610612741
    };

    cache_t *cache = xmalloc(sizeof(cache_t));

    int num_lookup;
    for(unsigned int i = 0; i < sizeof(primes)/sizeof(primes[0]); i++)
    {
        num_lookup = primes[i];
        if(entries < num_lookup*filling)
            break;
    }

    cache->num_entries = entries;
    cache->num_lookup = num_lookup;
    cache->head = entries-1;
    cache->tail = 0;
    cache->table  = xcalloc(num_lookup, sizeof(uint64_t));
    cache->keys   = xcalloc(entries, sizeof(uint64_t));
    cache->values = xcalloc(entries, sizeof(double));

    for(int i = 0; i < entries; i++)
        cache->keys[i] = 0;

    return cache;
}

void cache_free(cache_t *cache)
{
    xfree(cache->values);
    xfree(cache->keys);
    xfree(cache->table);
    xfree(cache);
}

void cache_insert(cache_t *cache, uint64_t key, double value)
{
    const int head = cache->head, tail = cache->tail;
    const int num_entries = cache->num_entries, num_lookup = cache->num_lookup;
    uint64_t *keys = cache->keys;
    double *values = cache->values;
    uint64_t *table = cache->table;

    cache->head = (head+1) % num_entries;
    cache->tail = (tail+1) % num_entries;

    const int64_t oldkey = keys[tail];
    keys[tail]   = key;
    values[tail] = value;

    table[oldkey % num_lookup] = 0;
    table[   key % num_lookup] = tail;
}

double cache_lookup(cache_t *cache, uint64_t key)
{
    const int index = cache->table[key % cache->num_lookup];

    if(cache->keys[index] == key)
        return cache->values[index];

    return NAN;
}

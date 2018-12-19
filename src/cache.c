/**
 * @file   cache.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   December, 2018
 * @brief  implementation of a simple cache using a hash table
 */

#include <stdint.h>
#include <math.h>

#include "utils.h"
#include "cache.h"

/**
 * @brief Create a new cache
 *
 * Create a new cache instance. This cache is quite specific. You specifiy the
 * maximum number of entries and a filling level. The cache is implemented as a
 * hash map that maps keys (uint64_t) to doubles.
 *
 * If the cache cannot contain more elements, the oldest entry will be thrown
 * away, similar to a FIFO. There is no logic to detect collisions. If there is
 * a collision, the old value will be overwritten. You can specifiy a filling
 * level (0 < filling < 1).
 *
 * @param entries maximum number of entries the cache can store
 * @param filling filling level
 * @retval cache cache_t instance
 */
cache_t *cache_new(int entries, double filling)
{
    /* http://planetmath.org/goodhashtableprimes */
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

/**
 * @brief Free cache instance
 *
 * @param cache cache instance
 */
void cache_free(cache_t *cache)
{
    xfree(cache->values);
    xfree(cache->keys);
    xfree(cache->table);
    xfree(cache);
}

/**
 * @brief Insert element into cache
 *
 * Insert the element value with key key to the cache.
 *
 * @param cache cache instance
 * @param key key
 * @param value value
 */
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

/**
 * @brief Find element corresponding to key key
 *
 * @param cache cache instance
 * @param key key
 * @retval element if found
 * @retval NAN otherwise
 */
double cache_lookup(cache_t *cache, uint64_t key)
{
    const int index = cache->table[key % cache->num_lookup];

    if(cache->keys[index] == key)
        return cache->values[index];

    return NAN;
}

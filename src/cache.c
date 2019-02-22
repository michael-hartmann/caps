/**
 * @file   cache.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   February, 2019
 * @brief  implementation of a simple cache using a hash table
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "cache.h"

/**
 * @brief Create a new cache
 *
 * Create a new cache instance.
 *
 * The cache is implemented as a hash map. Collisions are not treated, the
 * value will be overwritten.
 *
 * @param entries maximum number of entries the cache can store
 * @retval cache cache_t instance
 */
cache_t *cache_new(unsigned int entries)
{
    /* http://planetmath.org/goodhashtableprimes */
    const int primes[] = {
        1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433,
        1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319,
        201326611, 402653189, 805306457, 1610612741
    };

    cache_t *cache = xmalloc(sizeof(cache_t));

    /* determine number of entries */
    unsigned int num_entries;
    for(unsigned int i = 0; i < sizeof(primes)/sizeof(primes[0]); i++)
    {
        num_entries = primes[i];
        if(entries < num_entries)
            break;
    }

    cache->num_entries = num_entries;
    cache->table = xmalloc(num_entries*sizeof(cache_entry_t));

    for(unsigned int i = 0; i < num_entries; i++)
    {
        cache->table[i].key   = 0;
        cache->table[i].value = NAN;
    }

    return cache;
}

/**
 * @brief Free cache instance
 *
 * @param cache cache instance
 */
void cache_free(cache_t *cache)
{
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
    const unsigned int index = key % cache->num_entries;

    cache->table[index].key = key;
    cache->table[index].value = value;
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
    const unsigned int index = key % cache->num_entries;

    cache_entry_t entry = cache->table[index];
    if(entry.key == key)
        return entry.value;

    return NAN;
}

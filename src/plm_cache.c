/**
 * @file   plm_cache.c
 * @author Erik Buchenau <e.buchenau@live.de>
 * @date   March, 2016
 * @brief  Integration for legendre polynomials
 */
#include <stddef.h>
#include <assert.h>
#include <math.h>

#include "floattypes.h"
#include "integration_drude.h"
#include "libcasimir.h"
#include "utils.h"
#include "sfunc.h"
#include "plm_cache.h"

/*
 * The cache will save the values up to "max_iteration"
 */
#define max_iteration 15

/*
 * The cache will not save the values below "min_iteration"
 */
#define min_iteration 3

/*
 * Don't use this cache if the mean value of l1 and l2 is smaller than min_l
 */
#define min_l         2


/*
 * Count the number of cache-hits and cache-misses.
 * This is only useful for debugging and performance testing.
 */
//#define CACHE_STATS


/*
 * Debug this file
 */
//#define DEBUG_PLM_CACHE


enum cache_flags
{
    /*
     * The value of this entry is valid
     */
    CACHE_FLAG_VALID = 1U,
};


/*
 * This structure stores the legendre polynoms.
 */
struct cache_entry
{
    /*
     * m % 2. We need this, when we calculate the combinations.
     * They share the prefactor for their sign.
     */
    sign_t            common_sign;


    /*
     * We will save the logarithms and signs of legendre polynomials.
     * This will contain polynomials up the order MAX(l1, l2) + 1.
     */
    float80*          plm_array;
    sign_t*           plm_array_signs;
};


/*
 * This structure holds the values (and signs) of the legendre polynomials for a specific
 * iteration of the integration. Each iteraration i holds 2^i new points. They are
 * stored in the dynamically allocated space "entry".
 */
struct cache_iteration
{
    /*
     * Dynamically allocated memory for the cache_entries.
     * Iteration "i" should have 2^i cache_entries.
     */
    struct cache_entry* entry;
    
    /*
     * Flags from "enum cache_flags"
     */
    unsigned char       flags;
};


/*
 * This structure stores all the values (and signs) for one single integration.
 * The values can be accessed by indexing
 * e.g.:
 *     struct cache_values v;
 *     struct cache_entry e = &(v.iterations[iteration].entry[index])
 */
struct cache_values
{
    /*
     * We don't use the lower iterations, but we have a full array of cache_iterations
     * in memory, because this makes the indexing easier. And it doesn't take much
     * memory since we don't allocate space for the "struct cache_entry"s.
     */
    struct cache_iteration iterations[max_iteration + 1];

    /*
     * We will allocate the space for "plm_array" and "plm_array_signs" in
     * "struct cache_entry" at the beginning of the integration. The pointer of the
     * allocated memory space will be saved here, in order to free it up later.
     * We will allocate one big bulk of memory at the beginning of the integration.
     * This will save many (thousands) of small allocations, which can be problematic
     * for some memory allocators.
     * Don't use theese pointers directly (use "plm_array" and "plm_array_signs" in
     * "struct cache_entry" instead). Use this only to free memory at the end.
     */
    float80* plm_array_space;
    sign_t*  plm_sign_space;
};


/*
 * This is one complete cache. Each thread needs his own "struct plm_cache".
 * Don't create it on your own. Use plm_create_cache instead.
 * Use plm_destroy_cache if you don't need it anymore.
 */
struct plm_cache
{
    casimir_t* casimir;
    /*
     * The value for m. A cache can only save the values for one
     * m at once.
     */
    int m;

    /*
     * l1 + l2. We can reuse the previous results of plm if l1_plus_l2 is kept constant
     */
    int l1_plus_l2;
    
    /*
     * Matsubara term. Each thread should build up one object
     * of struct plm_cache
     */
    double nT;

    /*
     * Values of the legendre polynomials.
     */
    struct cache_values values;

#ifdef CACHE_STATS
    /*
     * Number of successfull cache accesses
     */
    unsigned long cache_hits;

    /*
     * Number of unsuccessfull cache accesses.
     * This is slow, because the values of plm need to be calculated using
     * the recurrence relations.
     */
    unsigned long cache_misses;
#endif
};


/* ****************************** Prototypes ****************************** */

static inline struct cache_entry* get_cache_value(struct integ_context* context,
                                                  float80 x,
                                                  struct plm_cache* cache,
                                                  unsigned int index,
                                                  unsigned int iteration);

static inline void calculate_cache_entry(struct integ_context* context,
                                         float80 x,
                                         struct cache_entry* entry);

static inline void plm_combination(struct integ_context* context,
                                   const float80 x,
                                   struct cache_entry* entry,
                                   plm_combination_t* comb);

/* ****************************** End Prototypes ****************************** */

/*
 * Calculate base^exp for integer values.
 * exp must be greater-equal zero
 */
static inline int pow_i(int base, int exp)
{
    int res = 1;
    while(exp > 0)
    {
        res *= base;
        exp--;
    }
    return res;
}


/*
 * Allocate the space for one object of the type "struct cache_values".
 * You need to call free_cache_values if you don't need it anymore.
 */
static void alloc_cache_values(struct cache_values* cv,
                               integration_drude_t* int_drude)
{
    size_t i;
    struct cache_iteration* ci;

    /*
     * We need all the polynomials from order m to order lmax + 1.
     * Here is the memory it takes.
     */
    size_t array_size = (int_drude->lmax - int_drude->m + 1) * sizeof(float80);
    size_t sign_size  = (int_drude->lmax - int_drude->m + 1) * sizeof(sign_t);

    /*
     * 2nd iteration will take 1 entry. 3 rd iteration will take 2. 4th 4, 5th 8 etc.
     * The total number of entries is:
     *    2^{min_iteration} + 2^{min_iteration + 1} + ... + 2^{max_iteration}
     * This is equal to:
     *    2^{max_iteration + 1} - 2^{min_iteration}
     */
    size_t total_entries = pow_i(2, max_iteration + 1) - pow_i(2, min_iteration);
    
    /*
     * Allocate memory for all the plm_arrays and plm_signs at once
     */
    float80* arrays = xmalloc(total_entries * array_size);
    sign_t*   signs = xmalloc(total_entries * sign_size);

    /*
     * We will save the pointers here, to free them up later.
     */
    cv->plm_array_space = arrays;
    cv->plm_sign_space  = signs;
    
    for(i = min_iteration; i <= max_iteration; ++i)
    {
        size_t entries = pow_i(2, i - 2);
        ci = &cv->iterations[i];
        /*
         * The uninitialized values at the moment are useless. So we won't set the
         * flag CACHE_FLAG_VALID.
         */
        ci->flags = 0U;
        /*
         * In each iteration i the romberg integration uses 2^(i-2) new points.
         * We will allocate the memory for it now.
         */
        ci->entry = xmalloc(entries * sizeof(struct cache_entry));

        
        for(size_t j = 0; j < entries; ++j)
        {
            ci->entry[j].plm_array = arrays;
            ci->entry[j].plm_array_signs = signs;
            arrays += (int_drude->lmax - int_drude->m + 1);
            signs += (int_drude->lmax - int_drude->m + 1);
        }

    }
}


/*
 * Free the space for one object of the type "struct cache_values".
 */
static void free_cache_values(struct cache_values* cv,
                              integration_drude_t* int_drude)
{
    struct cache_iteration* ci;

    /*
     * First, free up the big bulk. All the plm_arrays and signs in all cache_entrys
     * are going to be "free"d now.
     */
    xfree(cv->plm_array_space);
    xfree(cv->plm_sign_space);

    for(size_t i = max_iteration; i >= min_iteration; --i)
    {
        ci = &cv->iterations[i];
        xfree(ci->entry);
    }
}


/*
 * Invalidate the cache values. We need this if we are in a new cycle or if we
 * start to calculate a new matsubara term.
 */
static inline void inval_cache_values(struct cache_values* cv)
{
    size_t i;
    struct cache_iteration* ci;

    for(i = min_iteration; i <= max_iteration; ++i)
    {
        ci = &cv->iterations[i];
        ci->flags &= ~CACHE_FLAG_VALID;
    }
}


/*
 * Create a new "struct plm_cache". Each thread has to call this function in order
 * to use the plm_cache later.
 * In order to cleanup the allocated space, it is neccessary to call "plm_destroy_cach".
 */
void plm_create_cache(integration_drude_t* int_drude)
{
    int_drude->plm_cache = xmalloc(sizeof(struct plm_cache));
    struct plm_cache* cache = int_drude->plm_cache;    
    
    cache->m  = int_drude->m;
    cache->nT = int_drude->nT;
    
#ifdef CACHE_STATS
    cache->cache_hits   = 0;
    cache->cache_misses = 0;
#endif

    /*
     * Allocate the memory for the last and the first cache access
     */
    alloc_cache_values(&cache->values, int_drude);
}


/*
 * Free the space of the plm_cache. If you call "plm_create_cache" you have to call
 * this function at the end to avoid "dangling pointers".
 */
void plm_destroy_cache(integration_drude_t* int_drude)
{
    struct plm_cache* cache = int_drude->plm_cache;

    free_cache_values(&cache->values, int_drude);
    free(cache);
}


/*
 * Call this function before one complete integration. After that you can access the
 * polynomials with plm_cache_PlmPlm. At the end call plm_cache_free.
 * context is the integration context (see integration_drude.h).
 */
void plm_cache_init(struct integ_context* context)
{
    struct plm_cache* cache = context->int_drude->plm_cache;

    /*
     * If l2 == m is true, we have the first integration in one cycle.
     * We need to set up the cache for the new l1_plus_l2 etc.
     * If l1 + l2 or n or m changes, the last values are useless. We need to invalidate
     * the last values.
     */
    if(context->l2 == context->m
       || context->m != cache->m
       || context->l2 + context->l1 != cache->l1_plus_l2
       || cache->nT != context->nT)
    {
        cache->l1_plus_l2 = context->l1 + context->l2;
        cache->m          = context->m;
        cache->nT         = context->nT;
                
        /*
         * The last values are useless. We need to invalidate them
         */
        inval_cache_values(&cache->values);
    }
}


/*
 * Call this function after one single integration. 
 */
void plm_cache_free(integration_drude_t* int_drude)
{
#ifdef CACHE_STATS
    const struct plm_cache* cache = int_drude->plm_cache;
    unsigned long ratio;

    if(cache->cache_hits + cache->cache_misses == 0)
        ratio = 0UL;
    else
        ratio = 100UL * cache->cache_hits / (cache->cache_hits + cache->cache_misses);

    printf("Number:  %lu\nSuccess: %lu\nFailed:  %lu\nRatio:   %lu%%\n",
           cache->cache_hits + cache->cache_misses,
           cache->cache_hits,
           cache->cache_misses,
           ratio);
#endif
}


#ifdef DEBUG_PLM_CACHE

/*
 * Check the values of the cache access. It will calculate the "polynom combinations"
 * using plm_PlmPlm from sfunc.c and compare the result.
 * Only useful for debugging, because it will calculate the polynomials with
 * plm_PlmPlm. This is something we wan't to avoid with this cache.
 * x is the argument of the legendre polynomials.
 * comb is the result of the "polynom combinations" we got from the cache.
 * real_values is the "correct result" calculated with plm_PlmPlm
 */
static void check_values(int l1, int l2, int m, float80 x,
                         plm_combination_t* comb, plm_combination_t* real_values)
{
#define ALMOST_EQUAL(x,y) ((x) == (y) || fabs80(1 - (x)/(y)) < 1e-9)

    assert(ALMOST_EQUAL(real_values->lnPl1mPl2m, comb->lnPl1mPl2m));
    assert(ALMOST_EQUAL(real_values->lnPl1mdPl2m, comb->lnPl1mdPl2m));
    assert(ALMOST_EQUAL(real_values->lndPl1mPl2m, comb->lndPl1mPl2m));
    assert(ALMOST_EQUAL(real_values->lndPl1mdPl2m, comb->lndPl1mdPl2m));

    assert(real_values->sign_Pl1mPl2m == comb->sign_Pl1mPl2m);
    assert(real_values->sign_Pl1mdPl2m == comb->sign_Pl1mdPl2m);
    assert(real_values->sign_dPl1mPl2m == comb->sign_dPl1mPl2m);
    assert(real_values->sign_dPl1mdPl2m == comb->sign_dPl1mdPl2m);
}

#endif


/*
 * This function calculates combinations of legendre polynomials.
 * It tries to fetch the value from the cache. Otherwise it calculates it using
 * recurrence formulas.
 * Before using it, call plm_cache_init. Then access the values with this function
 * Finally call plm_cache_free.
 * Important: If you call this function for one specific iteration and one index,
 * you have to call it for the other indexes too. Otherwise the cache will contain
 * wrong values. E.g. you call index=1, iteration=4. You have to call index=2, iteration=4
 * and index=3, iteration=4 too. "romberg_In1" takes care of this.
 *
 * context   [in]:  The integration context of this integration
 * x         [in]:  Argument of the legendre polynomials
 * comb     [out]:  The result (combination of legendre polynomials)
 * index     [in]:  Index of the point for this iteration (see romberg_In1 in integration_drude.c)
 * iteration [in]:  Index of the iteration (see romberg_In1 in integration_drude.c)
 */
void plm_cache_PlmPlm(struct integ_context* context, float80 x, plm_combination_t* comb,
                      unsigned int index, unsigned int iteration)
{
    int l1, l2, m;
    struct plm_cache* cache;
    struct cache_entry* entry = NULL;

    l1 = context->l1;
    l2 = context->l2;
    m  = context->m;
    
    cache = context->int_drude->plm_cache;

    /*
     * If we calculate a bigger iteration than "max_iteration" something has gone wrong.
     * We should choose "max_iteration" as big as in "integration_drude.c".
     */
    assert(iteration <= max_iteration);
    
    if(l1 + l2 < 2 * min_l)
    {
        /*
         * Don't use the cache for the small multipoles. It won't speed up the calculation
         * either.
         */
        plm_PlmPlm(l1, l2, m, x, comb);
        return;
    }

    if(iteration >= min_iteration && iteration <= max_iteration)
    {
        /*
         * This is the "fast-path". Most calls of this function should be able to get the
         * values from the cache. This is much faster, than calling "plm_PlmPlm".
         */
        
        entry = get_cache_value(context, x, cache, index, iteration);
        plm_combination(context, x, entry, comb);

#ifdef DEBUG_PLM_CACHE
        /*
         * Theese are the "correct values" calculated with plm_PlmPlm.
         * We will check if they agree with the values from the cache.
         */
        plm_combination_t real_values;
        plm_PlmPlm(l1, l2, m, x, &real_values);
        
        check_values(l1, l2, m, x, comb, &real_values);
#endif
    }
    else
    {
        plm_PlmPlm(l1, l2, m, x, comb);
    }
}


/*
 * Returns the cache_entry from the cache.
 */
static inline struct cache_entry* get_cache_value(struct integ_context* context,
                                                  float80 x,
                                                  struct plm_cache* cache,
                                                  unsigned int index,
                                                  unsigned int iteration)
{
    struct cache_iteration* current_i;
    struct cache_entry* current;
    
    current_i  = &cache->values.iterations[iteration];

    if(!(current_i->flags & CACHE_FLAG_VALID))
    {
        /*
         * Can't use the last values. Need to calculate them.
         */
        calculate_cache_entry(context, x, &current_i->entry[index-1]);

#ifdef CACHE_STATS
        /*
         * Update the statistics of this cache
         */
        cache->cache_misses++;
#endif
    }
    else
    {
#ifdef CACHE_STATS
        /*
         * Update the statistics of this cache
         */
        cache->cache_hits++;
#endif
    }
    
    current  = &current_i->entry[index-1];
    current->common_sign = MPOW(context->m % 2);

    /*
     * Since the integration function calculates all the points for one iteration,
     * we can set all points in this iteration to valid (regardless of the index).
     */
    if((int)index == pow_i(2, iteration - 2))
        current_i->flags |= CACHE_FLAG_VALID;
    
    return current;
}


/*
 * Calculate Legendre Polynomials and its derivatives.
 * They are stored in a "struct cache_entry".
 * See "plm_PlmPlm" in sfunc.c.
 */
static inline void calculate_cache_entry(struct integ_context* context,
                                         float80 x,
                                         struct cache_entry* entry)
{
    const int l1   = context->l1;
    const int l2   = context->l2;
    const int lmax = MAX(l1, l2) + 1;
    const int m    = context->m;
    float80* lnPlm = entry->plm_array;
    sign_t*  signs = entry->plm_array_signs;

    plm_lnPlm_array(lmax, m, x, lnPlm, signs);
}


/*
 * Calculate combinations of legendre polynomials using a cache_entry.
 *
 * context   [in]:  Integration context
 * x         [in]:  Argument of the legendre polynomials
 * entry     [in]:  Cache entry with the single legendre plynomials
 * comb     [out]:  Combination of polynomials
 */
static inline void plm_combination(struct integ_context* context,
                                   const float80 x,
                                   struct cache_entry* entry,
                                   plm_combination_t* comb)
{
    float80 lndPl1;
    float80 lndPl2;

    sign_t  sign_lndPl1;
    sign_t  sign_lndPl2;

    const int l1              = context->l1;
    const int l2              = context->l2;
    const int m               = context->m;
    const float80 logx2m1     = log80(pow_2(x) - 1);
    const float80 logx        = log80(x);

    const float80* plm_array  = entry->plm_array;
    const sign_t*  plm_signs  = entry->plm_array_signs;


    lndPl1 = logadd_s( log80(l1-m+1) + plm_array[l1 - m + 1],
                       plm_signs[l1 - m + 1],
                       log80(l1+1) + logx + plm_array[l1 - m],
                       -plm_signs[l1 - m],
                       &sign_lndPl1) - logx2m1;

    lndPl2 = logadd_s( log80(l2-m+1) + plm_array[l2 - m + 1],
                       plm_signs[l2 - m + 1],
                       log80(l2+1) + logx + plm_array[l2 - m],
                       -plm_signs[l2 - m],
                       &sign_lndPl2) - logx2m1;
    
    /* Pl1m*Pl2m */
    comb->lnPl1mPl2m      = plm_array[l1 - m] + plm_array[l2 - m];
    comb->sign_Pl1mPl2m   = entry->common_sign * plm_signs[l1 - m] * plm_signs[l2 - m];

    /* Pl1m*dPl2m */
    comb->lnPl1mdPl2m     = plm_array[l1 - m] + lndPl2;
    comb->sign_Pl1mdPl2m  = entry->common_sign * plm_signs[l1 - m] * sign_lndPl2;

    /* dPl1m*Pl2m */
    comb->lndPl1mPl2m     = lndPl1 +plm_array[l2 - m];
    comb->sign_dPl1mPl2m  = entry->common_sign * sign_lndPl1 * plm_signs[l2 - m];

    /* dPl1m*dPl2m */
    comb->lndPl1mdPl2m    = lndPl1 + lndPl2;
    comb->sign_dPl1mdPl2m = entry->common_sign * sign_lndPl1 * sign_lndPl2;
}

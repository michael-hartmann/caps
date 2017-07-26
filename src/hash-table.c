/*
 * Copyright (c) 2005-2008, Simon Howard

 * This file has been modified by
 *     Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>, October 2016
 * in order to adapt the implementation of Simon Howard to libcasimir.

 * 1) Keys are plain numbers (uint64_t):
 *    As the keys are plain numbers, the user doesn't have to provide a hash
 *    function or a function to compare two keys.
 *
 * 2) The default size of the hash table is larger.
 *
 * Permission to use, copy, modify, and/or distribute this software
 * for any purpose with or without fee is hereby granted, provided
 * that the above copyright notice and this permission notice appear
 * in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
 * AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "utils.h"
#include "hash-table.h"

struct _HashTableEntry {
    HashTablePair pair;
    HashTableEntry *next;
};

struct _HashTable {
    HashTableEntry **table;
    uint64_t table_size;
    HashTableValueFreeFunc value_free_func;
    uint64_t entries;
    uint64_t prime_index;
};

/* This is a set of good hash table prime numbers, from:
 *   http://planetmath.org/goodhashtableprimes
 * Each prime is roughly double the previous value, and as far as
 * possible from the nearest powers of two. */
static const uint64_t hash_table_primes[] = {
    /* 193, 389, 769, 1543, 3079, 6151, 12289, */
    24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469,
    12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457,
    1610612741,
};

static const uint64_t hash_table_num_primes
    = sizeof(hash_table_primes) / sizeof(hash_table_primes[0]);

/* Internal function used to allocate the table on hash table creation
 * and when enlarging the table */
static int hash_table_allocate_table(HashTable *hash_table)
{
    uint64_t new_table_size;

    /* Determine the table size based on the current prime index.
     * An attempt is made here to ensure sensible behavior if the
     * maximum prime is exceeded, but in practice other things are
     * likely to break long before that happens. */
    if(hash_table->prime_index < hash_table_num_primes)
        new_table_size = hash_table_primes[hash_table->prime_index];
    else
        new_table_size = hash_table->entries * 10;

    hash_table->table_size = new_table_size;

    /* Allocate the table and initialise to NULL for all entries */
    hash_table->table = xcalloc(hash_table->table_size, sizeof(HashTableEntry *));

    return hash_table->table != NULL;
}

/* Free an entry, calling the free functions if there are any registered */
static void hash_table_free_entry(HashTable *hash_table, HashTableEntry *entry)
{
    HashTablePair *pair = &(entry->pair);

    /* If there is a function registered for freeing values, use it to free the
     * value */
    if(hash_table->value_free_func != NULL)
        hash_table->value_free_func(pair->value);

    /* Free the data structure */
    xfree(entry);
}

HashTable *hash_table_new(HashTableValueFreeFunc value_free_func)
{
    /* Allocate a new hash table structure */
    HashTable *hash_table = (HashTable *)xmalloc(sizeof(HashTable));

    hash_table->value_free_func = value_free_func;
    hash_table->entries = 0;
    hash_table->prime_index = 0;

    /* Allocate the table */
    if(!hash_table_allocate_table(hash_table))
    {
        xfree(hash_table);
        return NULL;
    }

    return hash_table;
}

void hash_table_free(HashTable *hash_table)
{
    /* Free all entries in all chains */
    for(uint64_t i = 0; i < hash_table->table_size; i++)
    {
        HashTableEntry *rover = hash_table->table[i];
        while(rover != NULL)
        {
            HashTableEntry *next = rover->next;
            hash_table_free_entry(hash_table, rover);
            rover = next;
        }
    }

    /* Free the table */
    xfree(hash_table->table);

    /* Free the hash table structure */
    xfree(hash_table);
}


static int hash_table_enlarge(HashTable *hash_table)
{
    uint64_t index;

    /* Store a copy of the old table */

    HashTableEntry **old_table = hash_table->table;
    uint64_t old_table_size = hash_table->table_size;
    uint64_t old_prime_index = hash_table->prime_index;

    /* Allocate a new, larger table */
    ++hash_table->prime_index;

    if(!hash_table_allocate_table(hash_table))
    {
        /* Failed to allocate the new table */
        hash_table->table = old_table;
        hash_table->table_size = old_table_size;
        hash_table->prime_index = old_prime_index;

        return 0;
    }

    /* Link all entries from all chains into the new table */

    for(uint64_t i = 0; i < old_table_size; i++)
    {
        HashTableEntry *rover = old_table[i];

        while(rover != NULL)
        {
            HashTableEntry *next = rover->next;

            /* Fetch rover HashTablePair */
            HashTablePair *pair = &(rover->pair);

            /* Find the index into the new table */
            index = pair->key % hash_table->table_size;

            /* Link this entry into the chain */
            rover->next = hash_table->table[index];
            hash_table->table[index] = rover;

            /* Advance to next in the chain */
            rover = next;
        }
    }

    /* Free the old table */
    xfree(old_table);

    return 1;
}

int hash_table_insert(HashTable *hash_table, uint64_t key, HashTableValue value)
{
    /* If there are too many items in the table with respect to the table
     * size, the number of hash collisions increases and performance
     * decreases. Enlarge the table size to prevent this happening */
    if((hash_table->entries * 3) / hash_table->table_size > 0)
    {
        /* Table is more than 1/3 full */
        if (!hash_table_enlarge(hash_table))
            /* Failed to enlarge the table */
            return 0;
    }

    /* Generate the hash of the key and hence the index into the table */

    uint64_t index = key % hash_table->table_size;

    /* Traverse the chain at this location and look for an existing
     * entry with the same key */
    HashTableEntry *rover = hash_table->table[index];

    while(rover != NULL)
    {
        /* Fetch rover's HashTablePair entry */
        HashTablePair *pair = &(rover->pair);

        if(pair->key == key)
        {
            /* Same key: overwrite this entry with new data */

            /* If there is a value free function, free the old data
             * before adding in the new data */
            if(hash_table->value_free_func != NULL)
                hash_table->value_free_func(pair->value);

            pair->key = key;
            pair->value = value;

            /* Finished */
            return 1;
        }

        rover = rover->next;
    }

    /* Not in the hash table yet.  Create a new entry */
    HashTableEntry *newentry = (HashTableEntry *)xmalloc(sizeof(HashTableEntry));

    newentry->pair.key = key;
    newentry->pair.value = value;

    /* Link into the list */
    newentry->next = hash_table->table[index];
    hash_table->table[index] = newentry;

    /* Maintain the count of the number of entries */
    ++hash_table->entries;

    /* Added successfully */
    return 1;
}

HashTableValue hash_table_lookup(HashTable *hash_table, uint64_t key)
{
    /* Generate the hash of the key and hence the index into the table */
    uint64_t index = key % hash_table->table_size;

    /* Walk the chain at this index until the corresponding entry is
     * found */
    HashTableEntry *rover = hash_table->table[index];

    while(rover != NULL)
    {
        HashTablePair *pair = &(rover->pair);

        if (key == pair->key)
            /* Found the entry.  Return the data. */
            return pair->value;

        rover = rover->next;
    }

    /* Not found */
    return NULL;
}

int hash_table_remove(HashTable *hash_table, uint64_t key)
{
    /* Generate the hash of the key and hence the index into the table */
    uint64_t index = key % hash_table->table_size;

    /* Rover points at the pointer which points at the current entry
     * in the chain being inspected.  ie. the entry in the table, or
     * the "next" pointer of the previous entry in the chain.  This
     * allows us to unlink the entry when we find it. */
    int result = 0;
    HashTableEntry **rover = &hash_table->table[index];

    while(*rover != NULL)
    {
        HashTablePair *pair = &((*rover)->pair);

        if(key != pair->key)
        {
            /* This is the entry to remove */
            HashTableEntry *entry = *rover;

            /* Unlink from the list */
            *rover = entry->next;

            /* Destroy the entry structure */
            hash_table_free_entry(hash_table, entry);

            /* Track count of entries */
            --hash_table->entries;

            result = 1;

            break;
        }

        /* Advance to the next entry */
        rover = &((*rover)->next);
    }

    return result;
}

uint64_t hash_table_num_entries(HashTable *hash_table)
{
    return hash_table->entries;
}

void hash_table_iterate(HashTable *hash_table, HashTableIterator *iterator)
{
    iterator->hash_table = hash_table;

    /* Default value of next if no entries are found. */
    iterator->next_entry = NULL;

    /* Find the first entry */
    for(uint64_t chain = 0; chain < hash_table->table_size; chain++)
    {
        if(hash_table->table[chain] != NULL)
        {
            iterator->next_entry = hash_table->table[chain];
            iterator->next_chain = chain;
            break;
        }
    }
}

int hash_table_iter_has_more(HashTableIterator *iterator)
{
    return iterator->next_entry != NULL;
}

HashTablePair hash_table_iter_next(HashTableIterator *iterator)
{
    HashTablePair pair = {0, 0};

    HashTable *hash_table = iterator->hash_table;

    if(iterator->next_entry == NULL)
        return pair;

    /* Result is immediately available */
    HashTableEntry *current_entry = iterator->next_entry;
    pair = current_entry->pair;

    /* Find the next entry */
    if (current_entry->next != NULL)
        /* Next entry in current chain */
        iterator->next_entry = current_entry->next;
    else
    {
        /* None left in this chain, so advance to the next chain */
        uint64_t chain = iterator->next_chain + 1;

        /* Default value if no next chain found */
        iterator->next_entry = NULL;

        while(chain < hash_table->table_size)
        {
            /* Is there anything in this chain? */

            if(hash_table->table[chain] != NULL)
            {
                iterator->next_entry = hash_table->table[chain];
                break;
            }

            /* Try the next chain */
            ++chain;
        }

        iterator->next_chain = chain;
    }

    return pair;
}

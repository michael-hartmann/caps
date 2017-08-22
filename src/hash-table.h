/*

Copyright (c) 2005-2008, Simon Howard

Permission to use, copy, modify, and/or distribute this software
for any purpose with or without fee is hereby granted, provided
that the above copyright notice and this permission notice appear
in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

 */

/**
 * @file hash-table.h
 *
 * @brief Hash table.
 *
 * A hash table stores a set of values which can be addressed by a
 * key.  Given the key, the corresponding value can be looked up
 * quickly.
 *
 * To create a hash table, use @ref hash_table_new.  To destroy a
 * hash table, use @ref hash_table_free.
 *
 * To insert a value into a hash table, use @ref hash_table_insert.
 *
 * To remove a value from a hash table, use @ref hash_table_remove.
 *
 * To look up a value by its key, use @ref hash_table_lookup.
 *
 * To iterate over all values in a hash table, use
 * @ref hash_table_iterate to initialise a @ref HashTableIterator
 * structure.  Each value can then be read in turn using
 * @ref hash_table_iter_next and @ref hash_table_iter_has_more.
 */

#ifndef ALGORITHM_HASH_TABLE_H
#define ALGORITHM_HASH_TABLE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * A hash table structure.
 */
typedef struct _HashTable HashTable;

/**
 * Structure used to iterate over a hash table.
 */
typedef struct _HashTableIterator HashTableIterator;

/**
 * Internal structure representing an entry in a hash table.
 */
typedef struct _HashTableEntry HashTableEntry;

/**
 * A value stored in a @ref HashTable.
 */
typedef void *HashTableValue;

/**
 * Internal structure representing an entry in hash table
 * used as @ref HashTableIterator next result.
 */
typedef struct _HashTablePair{
    uint64_t key;
    HashTableValue value;
} HashTablePair;

/**
 * Type of function used to free values when entries are removed from a
 * hash table.
 */
typedef void (*HashTableValueFreeFunc)(HashTableValue value);

/**
 * Create a new hash table.
 *
 * @param value_free_func      Function used to free values.
 * @return                     A new hash table structure, or NULL if it
 *                             was not possible to allocate the new hash
 *                             table.
 */
HashTable *hash_table_new(HashTableValueFreeFunc value_free_func);

/**
 * Destroy a hash table.
 *
 * @param hash_table           The hash table to destroy.
 */
void hash_table_free(HashTable *hash_table);

/**
 * Insert a value into a hash table, overwriting any existing entry
 * using the same key.
 *
 * @param hash_table           The hash table.
 * @param key                  The key for the new value.
 * @param value                The value to insert.
 * @return                     Non-zero if the value was added successfully,
 *                             or zero if it was not possible to allocate
 *                             memory for the new entry.
 */
int hash_table_insert(HashTable *hash_table, uint64_t key, HashTableValue value);

/**
 * Look up a value in a hash table by key.
 *
 * @param hash_table          The hash table.
 * @param key                 The key of the value to look up.
 * @return                    The value, or NULL if there is no value with that
 *                            key in the hash table.
 */
HashTableValue hash_table_lookup(HashTable *hash_table, uint64_t key);

/**
 * Remove a value from a hash table.
 *
 * @param hash_table          The hash table.
 * @param key                 The key of the value to remove.
 * @return                    Non-zero if a key was removed, or zero if the
 *                            specified key was not found in the hash table.
 */
int hash_table_remove(HashTable *hash_table, uint64_t key);

/**
 * Retrieve the number of entries in a hash table.
 *
 * @param hash_table          The hash table.
 * @return                    The number of entries in the hash table.
 */
uint64_t hash_table_num_entries(HashTable *hash_table);

#ifdef __cplusplus
}
#endif

#endif /* #ifndef ALGORITHM_HASH_TABLE_H */

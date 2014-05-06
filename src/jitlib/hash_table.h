#ifndef STRING_CACHE_H
#define STRING_CACHE_H

typedef struct hash_table_t *hash_table;
hash_table hash_table_create();

unsigned hash_table_size(hash_table table);
unsigned hash_table_capacity(hash_table table);
void hash_table_clear(hash_table table);

int* hash_table_get_n(hash_table table, char *key, int n);
char* hash_table_get_key_n(hash_table table, char *key, int n);
int* hash_table_insert_n(hash_table table, char *key, int n);

int* hash_table_get(hash_table table, char *key);
char* hash_table_get_key(hash_table table, char *key);
int* hash_table_insert(hash_table table, char *key);

void hash_table_destroy(hash_table table);

#endif

#include <stdlib.h>
#include <string.h>

static unsigned jenkins_one_at_a_time_hash(char *key, unsigned len) {
    unsigned hash, i;
    for(hash = i = 0;i<len && key[i] != '\0';++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

typedef struct bucket_t {
    unsigned hash;
    int value;
    char *key;
} bucket;

typedef struct hash_table_t {
    unsigned table_capacity, table_size;
    bucket *table;
} *hash_table;

unsigned hash_table_size(hash_table table) {
    return table->table_size;
}

unsigned hash_table_capacity(hash_table table) {
    return table->table_capacity;
}

static char* hash_table_add_string(hash_table table, char *str, int len) {
    (void)table;
    char *result = malloc(len+1);
    if(result == NULL) {
        return result;
    } else {
        memcpy(result, str, len);
        result[len] = '\0';
        return result;
    }
}

static void hash_table_allocate_content(hash_table table, unsigned capacity) {
    table->table_size = 0;
    table->table_capacity = capacity;
    table->table = calloc(table->table_capacity, sizeof(bucket));
}

static void hash_table_deallocate_content(hash_table table) {
    for(unsigned i = 0;i<table->table_capacity;++i) {
        free(table->table[i].key);
    }
    free(table->table);
}

hash_table hash_table_create() {
    hash_table table = malloc(sizeof(struct hash_table_t));
    hash_table_allocate_content(table, 15);
    return table;
}

void hash_table_clear(hash_table table) {
    for(unsigned i = 0;i<table->table_capacity;++i) {
        free(table->table[i].key);
    }
    memset(table->table, 0, sizeof(bucket)*table->table_capacity);
    table->table_size = 0;
}

static int hash_table_find(hash_table table, char *key, int len, int insert, int nocopy);

static void hash_table_grow(hash_table table) {
    struct hash_table_t newtable = *table;
    hash_table_allocate_content(&newtable, 2*(table->table_capacity+1)-1);
    for(unsigned i = 0;i<table->table_capacity;++i) {
        if(table->table[i].key != NULL) {
            char *key = table->table[i].key;
            int index = hash_table_find(&newtable, key, strlen(key), 1, 1);
            newtable.table[-index-1].value = table->table[i].value;
        }
    }
    free(table->table);
    *table = newtable;
}

static int hash_table_find(hash_table table, char *key, int len, int insert, int nocopy) {
    unsigned hash = jenkins_one_at_a_time_hash(key, len);
    int index = hash%table->table_capacity;
    for(;;) {
        if(table->table[index].key == NULL) {
            if(insert==1) {
                if(3*table->table_size > 2*table->table_capacity) {
                    hash_table_grow(table);
                    return hash_table_find(table, key, len, insert, nocopy);
                }
                table->table_size += 1;
                table->table[index].hash = hash;
                table->table[index].key = nocopy?key:hash_table_add_string(table, key, len);
                table->table[index].value = 0;
            }
            return -index-1;
        } else if(  table->table[index].hash == hash && 
                    strncmp(key, table->table[index].key, len) == 0 &&
                    table->table[index].key[len] == '\0') {
            return index;
        } else {
            index = (index+1)%table->table_capacity;
        }
    }
}

int* hash_table_get_n(hash_table table, char *key, int len) {
    int index = hash_table_find(table, key, len, 0, 0);
    if(index >= 0) {
        return &(table->table[index].value);
    } else {
        return NULL;
    }
}

char* hash_table_get_key_n(hash_table table, char *key, int len) {
    int index = hash_table_find(table, key, len, 0, 0);
    if(index >= 0) {
        return table->table[index].key;
    } else {
        return NULL;
    }
}

int* hash_table_insert_n(hash_table table, char *key, int len) {
    int index = hash_table_find(table, key, len, 1, 0);
    if(index >= 0) {
        return &(table->table[index].value);
    } else {
        return &(table->table[-index-1].value);
    }
}

int* hash_table_get(hash_table table, char *key) {
    return hash_table_get_n(table, key, strlen(key));
}

char* hash_table_get_key(hash_table table, char *key) {
    return hash_table_get_key_n(table, key, strlen(key));
}

int* hash_table_insert(hash_table table, char *key) {
    return hash_table_insert_n(table, key, strlen(key));
}

void hash_table_destroy(hash_table table) {
    hash_table_deallocate_content(table);
    free(table);
}

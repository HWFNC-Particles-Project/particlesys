#include <stdio.h>
#include <stdlib.h>
#include "malloc_align.h"

void *malloc_align(size_t size, size_t align, void **out_mem) {
    void *mem = malloc(size + align);
    if (mem == NULL) {
        if (out_mem != NULL) *out_mem = NULL;
        return NULL;
    }
    // align:
    uint64_t mem_i = (uint64_t)mem;
    uint64_t a_mask = (uint64_t)(align - 1);
    mem_i = (mem_i + a_mask) & ~a_mask;
    if (out_mem != NULL) *out_mem = mem;
    return (void *)mem_i;
}


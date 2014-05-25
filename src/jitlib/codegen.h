#ifndef CODEGEN_H
#define CODEGEN_H

#include <stdint.h>
#include "ssa_def.h"

typedef struct {
    uint32_t chunksize;
    uint32_t size;
} binary_header;

void free_code(void *fun);
void dump_code(void *fun);
void (*gencode_avx_ss(ssa_block *block, unsigned increment))(float*, float*, int);
void (*gencode_avx_ps(ssa_block *block, unsigned increment))(float*, float*, int);
void (*gencode_avx8_ps(ssa_block *block, unsigned increment))(float*, float*, int);

#endif

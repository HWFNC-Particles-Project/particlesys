#ifndef SSA_DEF_H
#define SSA_DEF_H

#include <stdint.h>
#include <stdlib.h>

static inline uint64_t float2uint(float x) {
    return *((uint32_t*)&x);
}

static inline float uint2float(uint64_t x) {
    x &= 0xFFFFFFFFull;
    return *((float*)&x);
}

enum ops_t {
    SSA_NONE, SSA_ADD, SSA_MUL, SSA_MIN, SSA_MAX, SSA_AND, SSA_XOR, SSA_OR,
    SSA_DIV, SSA_SUB, SSA_CMPEQ, SSA_CMPLT, SSA_CMPLE, SSA_CMPNEQ, SSA_CMPNLT, SSA_CMPNLE, SSA_ANDN,
    SSA_LOAD, SSA_CONST, SSA_STORE, SSA_SQRT, SSA_RSQRT, SSA_RCP,
    SSA_BLEND
};

#define DEAD_BIT (0x1000ull)
#define FLAG_BITS (0xF000ull)
#define OP_MASK (0x0FFFull)

#define GET_OP(x) ((x)&OP_MASK)
#define GET_CONST(x) (uint2float((x)>>16u))
#define GET_ARG1(x) (((x)>>16u)&0xFFFFu)
#define GET_ARG2(x) (((x)>>32u)&0xFFFFu)
#define GET_ARG3(x) (((x)>>48u)&0xFFFFu)
#define GET_ARGn(x, n) (((x)>>(16u*n))&0xFFFFu)

#define IS_UNARY(x) (GET_OP(x)>=SSA_STORE && GET_OP(x)<=SSA_RCP)
#define IS_BINARY(x) (GET_OP(x)>=SSA_ADD && GET_OP(x)<=SSA_ANDN)
#define IS_TERNARY(x) (GET_OP(x)>=SSA_BLEND && GET_OP(x)<=SSA_BLEND)
#define IS_COMMUTATIVE(x) (GET_OP(x)>=SSA_ADD && GET_OP(x)<=SSA_OR)
#define IS_LDSTCNST(x) (GET_OP(x)>=SSA_LOAD && GET_OP(x)<=SSA_STORE)
#define IS_ARITHMETIC(x) (!IS_LDSTCNST(x))
#define IS_COMPARISON(x) (GET_OP(x)>=SSA_CMPEQ && GET_OP(x)<=SSA_CMPNLE)

static inline uint64_t ssa_node(uint64_t op, uint64_t a, uint64_t b, uint64_t c) {
    return op | (a<<16u) | (b<<32u) | (c<<48u);
}
static inline uint64_t ssa_const(float x) {
    return SSA_CONST | (float2uint(x)<<16u);
}

typedef struct {
    int latency, throughput, port;
} ssa_scheduling_info;

typedef struct {
    size_t size, capacity;
    uint64_t *buffer;
    uint32_t *remap;
} ssa_block;


void ssa_block_load(ssa_block *block, uint64_t *code, size_t size);
uint64_t ssa_get_index(ssa_block *block, uint64_t i);
void ssa_remap_index(ssa_block *block, uint64_t i, uint64_t j);
void ssa_block_reserve(ssa_block *block, size_t capacity);
uint64_t ssa_apply_offset(uint64_t op, uint64_t offset);
uint64_t ssa_apply_remap_location(uint64_t op, const uint64_t *remap);
void ssa_block_append(ssa_block *block, ssa_block *block2, const uint64_t *remap);
void ssa_block_destroy(ssa_block *block);
uint64_t ssa_canonicalize_op(ssa_block *block, uint64_t op);
void ssa_apply_remap(ssa_block *block);
void ssa_apply_remap_raw(ssa_block *block);
void ssa_resolve(ssa_block *block);

#endif

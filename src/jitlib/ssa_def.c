#include <ssa_def.h>
#include <string.h>

void ssa_block_load(ssa_block *block, uint64_t *code, size_t size) {
    block->size = size;
    block->capacity = size;
    block->buffer = malloc(size*sizeof(uint64_t));
    block->remap = malloc(size*sizeof(uint32_t));
    memcpy(block->buffer, code, size*sizeof(uint64_t));
    for(uint32_t i = 0;i<size;++i) {
        block->remap[i] = i;
    }
}

uint64_t ssa_get_index(ssa_block *block, uint64_t i) {
    return block->remap[i];
}

void ssa_remap_index(ssa_block *block, uint64_t i, uint64_t j) {
    block->remap[i] = block->remap[j];
}

void ssa_remap_index_raw(ssa_block *block, uint64_t i, uint64_t j) {
    block->remap[i] = j;
}

void ssa_block_reserve(ssa_block *block, size_t capacity) {
    if(capacity < block->capacity) {
        return;
    }
    block->buffer = realloc(block->buffer, capacity*sizeof(uint64_t));
    block->remap = realloc(block->remap, capacity*sizeof(uint32_t));
    block->capacity = capacity;
}

uint64_t ssa_apply_offset(uint64_t op, uint64_t offset) {
    if(IS_UNARY(op)) {
        op += offset<<16u;
    } else if(IS_BINARY(op)) {
        op += offset<<16u | offset<<32u;
    } else if(IS_TERNARY(op)) {
        op += offset<<16u | offset<<32u | offset<<48u;
    }
    return op;
}

uint64_t ssa_apply_remap_location(uint64_t op, const uint64_t *remap) {
    if(GET_OP(op) == SSA_LOAD) {
        return ssa_node(GET_OP(op), remap[GET_ARG1(op)], 0, 0);
    } else if(GET_OP(op) == SSA_STORE) {
        return ssa_node(GET_OP(op), GET_ARG1(op), remap[GET_ARG2(op)], 0);
    } else {
        return op;
    }
}

void ssa_block_append(ssa_block *block, ssa_block *block2, const uint64_t *remap) {
    size_t size = block2->size;
    uint64_t *code = block2->buffer;
    ssa_block_reserve(block, block->size + size);
    for(uint32_t i = 0;i<size;++i) {
        block->buffer[block->size+i] = ssa_apply_offset(ssa_apply_remap_location(code[i], remap), block->size);
        block->remap[block->size+i] = block->size+i;
    }
    block->size += size;
}


void ssa_block_destroy(ssa_block *block) {
    free(block->buffer);
    free(block->remap);
}

int ssa_compare_op(ssa_block *block, uint64_t a, uint64_t b) {
    if(GET_OP(a) != GET_OP(b)) {
        return 0;
    }
    if(GET_OP(a) == SSA_CONST || GET_OP(a) == SSA_LOAD) {
        return a == b;
    }
    uint64_t a1 = ssa_get_index(block, GET_ARG1(a));
    uint64_t a2 = ssa_get_index(block, GET_ARG2(a));
    uint64_t a3 = ssa_get_index(block, GET_ARG3(a));

    uint64_t b1 = ssa_get_index(block, GET_ARG1(b));
    uint64_t b2 = ssa_get_index(block, GET_ARG2(b));
    uint64_t b3 = ssa_get_index(block, GET_ARG3(b));

    if(GET_OP(a) == SSA_STORE) {
        return (a1 == b1 && a2 == b2);
    } else if(IS_UNARY(a)) {
        return a1 == b1;
    } else if(IS_BINARY(a)) {
        if(IS_COMMUTATIVE(a)){
            return (a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1);
        } else {
            return (a1 == b1 && a2 == b2);
        }
    } else if(IS_TERNARY(a)) {
        return (a1 == b1 && a2 == b2 && a3 == b3);
    } else {
        return 0; // huh?
    }
}

void ssa_apply_remap(ssa_block *block) {
    for(size_t i = 0;i<block->size;++i) {
        if(IS_UNARY(block->buffer[i])) {
            uint64_t value = ssa_get_index(block, GET_ARG1(block->buffer[i]));
            block->buffer[i] &= ~(0xFFFFull<<16u);
            block->buffer[i] |= value<<16u;
        }

        if(IS_BINARY(block->buffer[i])) {
            uint64_t value = ssa_get_index(block, GET_ARG1(block->buffer[i]));
            block->buffer[i] &= ~(0xFFFFull<<16u);
            block->buffer[i] |= value<<16u;

            value = ssa_get_index(block, GET_ARG2(block->buffer[i]));
            block->buffer[i] &= ~(0xFFFFull<<32u);
            block->buffer[i] |= value<<32u;
        }

        if(IS_TERNARY(block->buffer[i])) {
            uint64_t value = ssa_get_index(block, GET_ARG1(block->buffer[i]));
            block->buffer[i] &= ~(0xFFFFull<<16u);
            block->buffer[i] |= value<<16u;

            value = ssa_get_index(block, GET_ARG2(block->buffer[i]));
            block->buffer[i] &= ~(0xFFFFull<<32u);
            block->buffer[i] |= value<<32u;

            value = ssa_get_index(block, GET_ARG3(block->buffer[i]));
            block->buffer[i] &= ~(0xFFFFull<<48u);
            block->buffer[i] |= value<<48u;
        }
    }
    for(uint32_t i = 0;i<block->size;++i) {
        block->remap[i] = i;
    }
}

void ssa_resolve(ssa_block *block) {
    ssa_apply_remap(block);
    size_t insert = 0;
    for(size_t i = 0;i<block->size;++i) {
        if(!(block->buffer[i] & DEAD_BIT)) {
            block->buffer[insert] = block->buffer[i];
            ssa_remap_index_raw(block, i, insert);
            ++insert;
        }
    }
    block->size = insert;
    ssa_apply_remap(block);
}

#include <stdio.h>
#include <limits.h>

#include <ssa_def.h>
#include <ssa_print.h>
#include <stdlib.h>
#include <math.h>

static void ssa_recurse_undead(ssa_block *block, size_t pos) {
    pos = ssa_get_index(block, pos);
    if(!(block->buffer[pos]&DEAD_BIT)) {
        return;
    }

    block->buffer[pos] &= ~DEAD_BIT;

    if(IS_UNARY(block->buffer[pos])) {
        ssa_recurse_undead(block, GET_ARG1(block->buffer[pos]));
    } else if(IS_BINARY(block->buffer[pos])) {
        ssa_recurse_undead(block, GET_ARG1(block->buffer[pos]));
        ssa_recurse_undead(block, GET_ARG2(block->buffer[pos]));
    } else if(IS_TERNARY(block->buffer[pos])) {
        ssa_recurse_undead(block, GET_ARG1(block->buffer[pos]));
        ssa_recurse_undead(block, GET_ARG2(block->buffer[pos]));
        ssa_recurse_undead(block, GET_ARG3(block->buffer[pos]));
    }
}

void ssa_mark_dead(ssa_block *block) {
    size_t maxstore = 0;
    for(size_t i = 0;i<block->size;++i) {
        if(GET_OP(block->buffer[i]) == SSA_STORE && GET_ARG2(block->buffer[i]) > maxstore) {
            maxstore = GET_ARG2(block->buffer[i]);
        }
        block->buffer[i] |= DEAD_BIT;
    }

    int stored[maxstore+1];
    for(size_t j = 0;j<=maxstore;++j) {
        stored[j] = 0;
    }

    for(size_t j = 1;j<=block->size;++j) {
        size_t i = block->size-j;
        if(GET_OP(block->buffer[i]) == SSA_STORE && !stored[GET_ARG2(block->buffer[i])]) {
            ssa_recurse_undead(block, i);
            stored[GET_ARG2(block->buffer[i])] = 1;
        }
    }
}

size_t find_index(ssa_block *block, int *hash_map, size_t C, uint64_t op) {
    const uint64_t p = 9223372036854775783ull;
    uint64_t x = ssa_canonicalize_op(block, op);
    size_t index = (p*x)%C;
    while(hash_map[index]>=0 && x != ssa_canonicalize_op(block, block->buffer[ssa_get_index(block, hash_map[index])])) {
        index = (index+1)%C;
    }
    return index;
}

void ssa_remap_duplicates(ssa_block *block, size_t window) {
    size_t C = 2*block->size;
    int hash_map[C];
    for(size_t i = 0;i<C;++i) {
        hash_map[i] = -1;
    }

    for(size_t i = 0;i<block->size;++i) {
        size_t i1 = ssa_get_index(block, i);
        size_t index = find_index(block, hash_map, C, block->buffer[i1]);
        if(hash_map[index]>=0) {
            ssa_remap_index(block, i1, ssa_get_index(block, hash_map[index]));
        } else {
            hash_map[index] = i1;
        }
    }
}

void ssa_fuse_load_store(ssa_block *block) {
    size_t maxinout = 0;
    for(size_t i = 0;i<block->size;++i) {
        uint64_t op = GET_OP(block->buffer[i]);
        uint64_t arg1 = GET_ARG1(block->buffer[i]);
        uint64_t arg2 = GET_ARG2(block->buffer[i]);

        if(op == SSA_LOAD && arg1 > maxinout) {
            maxinout = arg1;
        }
        if(op == SSA_STORE && arg2 > maxinout) {
            maxinout = arg2;
        }
    }

    int scratch[maxinout+1];
    for(size_t i = 0;i<maxinout+1;++i) {
        scratch[i] = -1;
    }
    for(size_t i = 0;i<block->size;++i) {
        uint64_t op = GET_OP(block->buffer[i]);
        uint64_t arg1 = GET_ARG1(block->buffer[i]);
        uint64_t arg2 = GET_ARG2(block->buffer[i]);

        if(op == SSA_LOAD) {
            if(scratch[arg1] != -1) {
                ssa_remap_index(block, i, scratch[arg1]);
            }
        } else if(op == SSA_STORE) {
            scratch[arg2] = arg1;
        }
    }
}

void ssa_fold_constants(ssa_block *block) {
    for(size_t i = 0;i<block->size;++i) {
        uint64_t op = block->buffer[i];
        if(IS_BINARY(op)) {
            uint64_t left = block->buffer[ssa_get_index(block, GET_ARG1(op))];
            uint64_t right = block->buffer[ssa_get_index(block, GET_ARG2(op))];
            if(GET_OP(left) == SSA_CONST && GET_OP(right) == SSA_CONST) {
                switch(GET_OP(op)) {
                    case SSA_ADD:
                        block->buffer[i] = ssa_const(GET_CONST(left) + GET_CONST(right));
                        break;
                    case SSA_SUB:
                        block->buffer[i] = ssa_const(GET_CONST(left) - GET_CONST(right));
                        break;
                    case SSA_MUL:
                        block->buffer[i] = ssa_const(GET_CONST(left) * GET_CONST(right));
                        break;
                    case SSA_DIV:
                        block->buffer[i] = ssa_const(GET_CONST(left) / GET_CONST(right));
                        break;
                    case SSA_MIN:
                        block->buffer[i] = ssa_const(fminf(GET_CONST(left),GET_CONST(right)));
                        break;
                    case SSA_MAX:
                        block->buffer[i] = ssa_const(fmaxf(GET_CONST(left),GET_CONST(right)));
                        break;
                }
            }
        } else if(IS_UNARY(op)){
            uint64_t left = block->buffer[ssa_get_index(block, GET_ARG1(op))];
            if(GET_OP(left) == SSA_CONST) {
                switch(GET_OP(op)) {
                    case SSA_SQRT:
                        block->buffer[i] = ssa_const(sqrtf(GET_CONST(left)));
                        break;
                    case SSA_RSQRT:
                        break;
                    case SSA_RCP:
                        break;
                    default:
                        break;
                }
            }
        }
    }
}

typedef struct  {
    int *data;
    int size, capacity;
} set;

void set_reserve(set *s, int cap) {
    if(s->capacity<cap) {
        s->capacity = cap;
        s->data = realloc(s->data, sizeof(int)*s->capacity);
    }
}

void set_init(set *s, int capacity) {
    s->size = 0;
    s->capacity = 0;
    s->data = NULL;
    set_reserve(s, capacity);
}

void set_free(set *s) {
    free(s->data);
}

int set_find(set *s, int value) {
   for(int i = 0;i<s->size;++i) {
        if(s->data[i] == value) {
            return 1;
        }
   }
   return 0;
}

void set_insert(set *s, int value) {
    if(set_find(s, value)) {
        return;
    } else {
        if(s->size == s->capacity) {
            set_reserve(s, s->capacity*2);
        }
        s->data[s->size++] = value;
    }
}

void set_remove(set *s, int value) {
    int i ,j;
    for(i = 0, j = 0;i<s->size;++i) {
        if(s->data[i] != value) {
            s->data[j] = s->data[i];
            ++j;
        }
    }
    s->size = j;
}

static void ssa_add_latencies(ssa_block *block, int *latencies, int i, int latency, const ssa_scheduling_info *info) {
    uint64_t op = block->buffer[i];
    latency += info[GET_OP(op)].latency;
    if(latencies[i] < latency) {
        latencies[i] = latency;
        if(IS_UNARY(op)) {
            ssa_add_latencies(block, latencies, GET_ARG1(op), latency, info);
        }
        if(IS_BINARY(op)) {
            ssa_add_latencies(block, latencies, GET_ARG1(op), latency, info);
            ssa_add_latencies(block, latencies, GET_ARG2(op), latency, info);
        }
        if(IS_TERNARY(op)) {
            ssa_add_latencies(block, latencies, GET_ARG1(op), latency, info);
            ssa_add_latencies(block, latencies, GET_ARG2(op), latency, info);
            ssa_add_latencies(block, latencies, GET_ARG3(op), latency, info);
        }
    }
}

void ssa_schedule(ssa_block *block, const ssa_scheduling_info *info) {
    ssa_apply_remap(block);

    uint64_t *newbuf = malloc(block->capacity*sizeof(uint64_t));
    size_t scheduled[block->size];
    int latencies[block->size];
    int uses[block->size];
    int registers = 0;

    int total_uses = 0;
    for(size_t i = 0;i<block->size;++i) {
        uint64_t op = block->buffer[i];
        scheduled[i] = block->size;
        latencies[i] = -1;
        uses[i] = 0;
        if(GET_OP(op)==SSA_STORE) {
            ssa_add_latencies(block, latencies, i, 0, info);
            ++uses[GET_ARG1(op)];
            total_uses += 1;
        } else if(IS_ARITHMETIC(op)) {
            if(IS_UNARY(op)) {
                ++uses[GET_ARG1(op)];
                total_uses += 1;
            } else if(IS_BINARY(op)) {
                ++uses[GET_ARG1(op)];
                ++uses[GET_ARG2(op)];
                total_uses += 2;
            } else if(IS_TERNARY(op)) {
                ++uses[GET_ARG1(op)];
                ++uses[GET_ARG2(op)];
                ++uses[GET_ARG3(op)];
                total_uses += 3;
            }
        }
    }

    int start[block->size+1];
    start[block->size] = total_uses;

    int dependencies[block->size];
    int graph[total_uses];
    int sum = 0;
    int independents = 0;
    for(size_t i = 0;i<block->size;++i) {
        sum += uses[i];
        start[i] = sum;
        dependencies[i] = 0;
        uint64_t op = block->buffer[i];
        if(GET_OP(op)==SSA_STORE) {
            graph[--start[GET_ARG1(op)]] = i;
            dependencies[i] = 1;
        } else if(IS_ARITHMETIC(op)) {
            if(IS_UNARY(op)) {
                graph[--start[GET_ARG1(op)]] = i;
                dependencies[i] = 1;
            } else if(IS_BINARY(op)) {
                graph[--start[GET_ARG1(op)]] = i;
                graph[--start[GET_ARG2(op)]] = i;
                dependencies[i] = 2;
            } else if(IS_TERNARY(op)) {
                graph[--start[GET_ARG1(op)]] = i;
                graph[--start[GET_ARG2(op)]] = i;
                graph[--start[GET_ARG3(op)]] = i;
                dependencies[i] = 3;
            }
        } else {
            ++independents;
        }
    }

    set open;
    set_init(&open, independents);
    for(size_t i = 0;i<block->size;++i) {
        if(dependencies[i] == 0){
            set_insert(&open, i);
        }
    }

    for(size_t i = 0;i<block->size;++i) {
        size_t select = 0;
        int best = INT_MIN;
        int dreg = 100;
        for(size_t k = 0;k<open.size;++k) {
            size_t j = open.data[k];
            int score = 0;
            int registerchange = 0;
            uint64_t op = block->buffer[j];
            if(IS_ARITHMETIC(op)) {
                if(IS_UNARY(op)) {
                    if(scheduled[GET_ARG1(op)] < i) {
                        score += 1+latencies[j];
                        registerchange = 1 - (uses[GET_ARG1(op)]==1?1:0);
                    }
                } else if(IS_BINARY(op)) {
                    if(scheduled[GET_ARG1(op)] < i && scheduled[GET_ARG2(op)] < i) {
                        score += 1+latencies[j];
                        registerchange = 1 - (uses[GET_ARG1(op)]==1?1:0) - (uses[GET_ARG2(op)]==1?1:0);
                    }
                } else if(IS_TERNARY(op)) {
                    if(scheduled[GET_ARG1(op)] < i && scheduled[GET_ARG2(op)] < i && scheduled[GET_ARG3(op)] < i) {
                        score += 1+latencies[j];
                        registerchange = 1 - (uses[GET_ARG1(op)]==1?1:0) - (uses[GET_ARG2(op)]==1?1:0) - (uses[GET_ARG3(op)]==1?1:0);
                    }
                }
            } else if(GET_OP(op)==SSA_STORE) {
                if(scheduled[GET_ARG1(op)] < i) {
                    score += 1+latencies[j];
                    registerchange = -(uses[GET_ARG1(op)]==1?1:0);
                }
            } else if(GET_OP(op)==SSA_LOAD || GET_OP(op)==SSA_CONST) {
                score += 1+latencies[j];
                registerchange = 1;
            }
            if(registers >= 15) {
                score += 1000*(1-registerchange);
            }
            if(best<score) {
                select = j;
                dreg = registerchange;
                best = score;
            }
        }

        registers += dreg;
        uint64_t op = block->buffer[select];

        if(IS_UNARY(op)) {
            --uses[GET_ARG1(op)];
        } else if(IS_BINARY(op)) {
            --uses[GET_ARG1(op)];
            --uses[GET_ARG2(op)];
        } else if(IS_TERNARY(op)) {
            --uses[GET_ARG1(op)];
            --uses[GET_ARG2(op)];
            --uses[GET_ARG3(op)];
        }

        set_remove(&open, select);
        for(int j = start[select];j<start[select+1];++j) {
            int k = graph[j];
            if(--dependencies[k] == 0) {
                set_insert(&open, k);
            }
        }

        scheduled[select] = i;
        ssa_remap_index(block, select, i);
        newbuf[i] = op;
    }

    set_free(&open);

    free(block->buffer);
    block->buffer = newbuf;
    ssa_apply_remap_raw(block);
}

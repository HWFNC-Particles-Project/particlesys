#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include <sys/mman.h>

#include <plasm.h>
#include "codegen.h"

static mnemonic_t ssa2asm_ss[] = {
    [SSA_NONE]   = UNKNOWN,
    [SSA_ADD]    = VADDPS,
    [SSA_MUL]    = VMULPS,
    [SSA_MIN]    = VMINPS,
    [SSA_MAX]    = VMAXPS,
    [SSA_AND]    = VANDPS,
    [SSA_ANDN]   = VANDNPS,
    [SSA_XOR]    = VXORPS,
    [SSA_OR]     = VORPS,
    [SSA_DIV]    = VDIVPS,
    [SSA_SUB]    = VSUBPS,
    [SSA_SQRT]   = VSQRTPS,
    [SSA_RSQRT]  = VRSQRTPS,
    [SSA_RCP]    = VRCPPS,
    [SSA_BLEND]  = VBLENDVPS,
};

static uint8_t cmpimm[] = {
    [SSA_CMPEQ]  = 0x00,
    [SSA_CMPLT]  = 0x01,
    [SSA_CMPLE]  = 0x02,
    [SSA_CMPNEQ] = 0x04,
    [SSA_CMPNLT] = 0x05,
    [SSA_CMPNLE] = 0x06,
};

typedef struct {
    size_t pos;
    int *regs, *stackpos, *uses, *free_regs, *locked, *owner;
    int top, spillpos;
    unsigned increment;
    plasm *as;
    ssa_block *block;
    int *start, *graph;
} codegen_data;

static opspec_t allocate_temp_register(codegen_data *data);

static void gather8(codegen_data *data, opspec_t target, opspec_t idx, int offset, int stride) {
    opspec_t staging = allocate_temp_register(data);
    int regidx = target.type >> 16ull;
    int stageidx = staging.type >> 16ull;

    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+0*stride), IMM8(0x00));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+1*stride), IMM8(0x10));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+2*stride), IMM8(0x20));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+3*stride), IMM8(0x30));
    plasm_put_op(data->as, VINSERTF128, YMM(regidx), YMM(regidx), XMM(stageidx), IMM8(0));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+4*stride), IMM8(0x00));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+5*stride), IMM8(0x10));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+6*stride), IMM8(0x20));
    plasm_put_op(data->as, VINSERTPS, XMM(stageidx), XMM(stageidx), MEM(idx, offset+7*stride), IMM8(0x30));
    plasm_put_op(data->as, VINSERTF128, YMM(regidx), YMM(regidx), XMM(stageidx), IMM8(1));
}

static void scatter8(codegen_data *data, opspec_t source, opspec_t idx, int offset, int stride) {
    int regidx = source.type >> 16ull;

    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+0*stride), XMM(regidx), IMM8(0));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+1*stride), XMM(regidx), IMM8(1));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+2*stride), XMM(regidx), IMM8(2));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+3*stride), XMM(regidx), IMM8(3));
    plasm_put_op(data->as, VEXTRACTF128, XMM(regidx), YMM(regidx), IMM8(1));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+4*stride), XMM(regidx), IMM8(0));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+5*stride), XMM(regidx), IMM8(1));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+6*stride), XMM(regidx), IMM8(2));
    plasm_put_op(data->as, VEXTRACTPS, MEM(idx, offset+7*stride), XMM(regidx), IMM8(3));
}

static int next_use(codegen_data *data, int i, int pos) {
    for(size_t j = data->start[i];j<data->start[i+1];++j) {
        if(data->graph[j]>pos) {
            return data->graph[j];
        }
    }
    return -1;
}

static opspec_t evict_operand(codegen_data *data, int i) {
    opspec_t src = YMM(data->regs[i]-1);

    data->free_regs[data->top++] = data->regs[i];
    data->owner[data->regs[i]-1] = -1;

    if(data->stackpos[i] == 0) {
        data->stackpos[i] = data->spillpos--;
        data->regs[i] = data->stackpos[i];
        opspec_t dest = MEM(RSP, 32*(data->regs[i]));
        plasm_put_op(data->as, VMOVUPS, dest, src);
    } else {
        data->regs[i] = data->stackpos[i];
    }
    return src;
}

static opspec_t spill_operand(codegen_data *data) {
    int i = 0;
    int next = 0;
    for(int j = 0;j<16;++j) {
        if(data->owner[j]>=0 && !data->locked[j]) {
            int use = next_use(data, data->owner[j], data->pos);
            if(data->stackpos[data->owner[j]] != 0) {
                use += 1000;
            }
            if(use>next) {
                next = use;
                i = data->owner[j];
            }
        }
    }

    return evict_operand(data, i);
}

static opspec_t allocate_temp_register(codegen_data *data) {
    if(data->top <= 0) {
        return spill_operand(data);
    }
    return YMM(data->free_regs[data->top-1]-1);
}

static opspec_t allocate_operand(codegen_data *data, int i) {
    if(data->top <= 0) {
        spill_operand(data);
    }
    data->regs[i] = data->free_regs[--data->top];
    data->owner[data->regs[i]-1] = i;
    return YMM(data->regs[i]-1);
}

static opspec_t lock_operand(codegen_data *data, int i, int mem) {
    if(data->regs[i]>0) {
        data->locked[data->regs[i]-1] = 1;
        return YMM(data->regs[i]-1);
    } else {
        uint64_t src = data->block->buffer[i];
        uint64_t location = data->regs[i];

        opspec_t target;
        if(GET_OP(src) == SSA_LOAD) {
            target = allocate_operand(data, i);
            if(GET_ARG1(src)<data->increment) {
                gather8(data, target, RDI, 4*GET_ARG1(src), 4*data->increment);
            } else {
                plasm_put_op(data->as, VBROADCASTSS, target, MEM(RSI, 4*(GET_ARG1(src)-data->increment)));
            }
        } else {
            if(mem) {
                return MEM(RSP, 32*location);
            } else {
                target = allocate_operand(data, i);
                plasm_put_op(data->as, VMOVUPS, target, MEM(RSP, 32*location));
            }
        }
        data->locked[data->regs[i]-1] = 1;
        return target;
    }
}

static void release_operand(codegen_data *data, int i) {
    if(data->regs[i]>0) {
        data->locked[data->regs[i]-1] = 0;
    }
    if(--data->uses[i] == 0 && data->regs[i]>0) {
        data->free_regs[data->top++] = data->regs[i];
        data->owner[data->regs[i]-1] = -1;
    }
}

void (*gencode_avx8_ps(ssa_block *block, unsigned increment))(float*, float*, int)  {
    int uses[block->size];
    int total_uses = 0;
    for(size_t i = 0;i<block->size;++i) {
        uint64_t op = block->buffer[i];
        uses[i] = 0;
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

    int start[block->size+1];
    start[block->size] = total_uses;

    int dependencies[block->size];
    int graph[total_uses];
    int sum = 0;
    for(size_t i = 0;i<block->size;++i) {
        start[i] = sum;
        sum += uses[i];
        dependencies[i] = 0;
        uint64_t op = block->buffer[i];
        if(IS_UNARY(op)) {
            graph[start[GET_ARG1(op)]++] = i;
            dependencies[i] = 1;
        } else if(IS_BINARY(op)) {
            graph[start[GET_ARG1(op)]++] = i;
            graph[start[GET_ARG2(op)]++] = i;
            dependencies[i] = 2;
        } else if(IS_TERNARY(op)) {
            graph[start[GET_ARG1(op)]++] = i;
            graph[start[GET_ARG2(op)]++] = i;
            graph[start[GET_ARG3(op)]++] = i;
            dependencies[i] = 3;
        }
    }
    sum = 0;
    for(size_t i = 0;i<block->size;++i) {
        start[i] = sum;
        sum += uses[i];
    }

    int regs[block->size];
    int stackpos[block->size];
    int free_regs[16] = {16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    int locked[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int owner[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int top = 16;

    uint32_t chunksize = ((block->size*24u)+4095u)&(~4095u);
    binary_header *chunk = mmap(NULL, chunksize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    uint8_t *buffer = (uint8_t*)(chunk+1);
    chunk->chunksize = chunksize;

    plasm as;
    plasm_init(&as, buffer, chunksize-sizeof(binary_header));
    void (*fun)(float*, float*, int) = (void(*)(float*, float*, int))plasm_get_current_ptr(&as);

    codegen_data data = {0, regs, stackpos, uses, free_regs, locked, owner, top, -1, increment, &as, block, start, graph};

    float *constbase = (float*)(((uint8_t*)chunk) + chunksize);
    int constindex = -1;

    plasm_put_op(&as, PUSH, RBP);
    plasm_put_op(&as, MOV, RBP, RSP);
    plasm_put_op(&as, SHL, RDX, IMM8(2));
    plasm_put_op(&as, ADD, RDX, RDI);
    int32_t *endptr;
    plasm_put_op(&as, LEA, RAX, VALUE_PTR(DISP_RIP(0), &endptr));
    *endptr = (uint8_t*)constbase-plasm_get_current_ptr(&as);

    uint8_t *loop = plasm_get_current_ptr(&as);

    for(size_t i = 0;i<block->size;++i) {
        regs[i] = 0;
        stackpos[i] = 0;
        uint64_t op = block->buffer[i];
        if(GET_OP(op) == SSA_LOAD) {

        } else if(GET_OP(op) == SSA_STORE) {
            opspec_t source = lock_operand(&data, GET_ARG1(op), 0);

            scatter8(&data, source, RDI, 4*GET_ARG2(op), 4*increment);

            release_operand(&data, GET_ARG1(op));
        } else if(GET_OP(op) == SSA_CONST) {
            opspec_t target = allocate_operand(&data, i);

            *(constbase+constindex) = GET_CONST(op);

            plasm_put_op(&as, VBROADCASTSS, target, MEM(RAX, 4*constindex));
            --constindex;

        } else if(IS_COMPARISON(op)) {
            int arg1 = GET_ARG1(op);
            int arg2 = GET_ARG2(op);
            opspec_t src1 = lock_operand(&data, arg1, 0);
            opspec_t src2 = lock_operand(&data, arg2, 1);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, VCMPPS, target, src1, src2, IMM8(cmpimm[GET_OP(op)]));

            release_operand(&data, arg2);
            release_operand(&data, arg1);
        } else if(IS_UNARY(op)) {
            int arg1 = GET_ARG1(op);
            opspec_t src1 = lock_operand(&data, arg1, 1);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, ssa2asm_ss[GET_OP(op)], target, src1);

            release_operand(&data, arg1);
        } else if(IS_BINARY(op)) {
            int arg1 = GET_ARG1(op);
            int arg2 = GET_ARG2(op);
            if(IS_COMMUTATIVE(op) && data.regs[arg1]<=0) {
                int tmp = arg1; arg1 = arg2; arg2 = tmp;
            }
            opspec_t src1 = lock_operand(&data, arg1, 0);
            opspec_t src2 = lock_operand(&data, arg2, 1);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, ssa2asm_ss[GET_OP(op)], target, src1, src2);

            release_operand(&data, arg2);
            release_operand(&data, arg1);
        } else if(IS_TERNARY(op)) {
            int arg1 = GET_ARG1(op);
            int arg2 = GET_ARG2(op);
            int arg3 = GET_ARG3(op);

            opspec_t src1 = lock_operand(&data, arg1, 0);
            opspec_t src2 = lock_operand(&data, arg2, 1);
            opspec_t src3 = lock_operand(&data, arg3, 0);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, ssa2asm_ss[GET_OP(op)], target, src3, src2, src1);

            release_operand(&data, arg3);
            release_operand(&data, arg2);
            release_operand(&data, arg1);
        }
    }

    plasm_put_op(&as, ADD, RDI, IMM32(8*4*increment));
    plasm_put_op(&as, CMP, RDI, RDX);
    int32_t *reloff_ptr;
    plasm_put_op(&as, JNE, VALUE_PTR(REL32(0), &reloff_ptr));
    *reloff_ptr = loop - plasm_get_current_ptr(&as);

    plasm_put_op(&as, POP, RBP);
    plasm_put_op(&as, RET);

    chunk->size = as.position;

    mprotect(chunk, chunksize, PROT_EXEC|PROT_READ);

    return fun;
}

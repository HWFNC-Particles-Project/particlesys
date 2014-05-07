#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include <sys/mman.h>

#include <plasm.h>
#include <codegen.h>
#include <ssa_print.h>

static mnemonic_t ssa2asm_ss[] = {
    [SSA_NONE]   = UNKNOWN,
    [SSA_ADD]    = VADDSS,
    [SSA_MUL]    = VMULSS,
    [SSA_MIN]    = VMINSS,
    [SSA_MAX]    = VMAXSS,
    [SSA_AND]    = VANDPS,
    [SSA_ANDN]   = VANDNPS,
    [SSA_XOR]    = VXORPS,
    [SSA_OR]     = VORPS,
    [SSA_DIV]    = VDIVSS,
    [SSA_SUB]    = VSUBSS,
    [SSA_SQRT]   = VSQRTSS,
    [SSA_RSQRT]  = VRSQRTSS,
    [SSA_RCP]    = VRCPSS,
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
    int *regs, *uses, *free_regs, *locked, *owner;
    int top, spillpos;
    unsigned increment;
    plasm *as;
    ssa_block *block;
    int *start, *graph;
} codegen_data;

int next_use(codegen_data *data, int i, int pos) {
    for(size_t j = data->start[i];j<data->start[i+1];++j) {
        if(data->graph[j]>pos) {
            return data->graph[j];
        }
    }
    return -1;
}

static void spill_operand(codegen_data *data) {
    int i = 0;
    int next = 0;
    for(int j = 0;j<16;++j) {
        if(data->owner[j]>=0 && !data->locked[j]) {
            int use = next_use(data, data->owner[j], data->pos);
            if(use>next) {
                next = use;
                i = data->owner[j];
            }
        }
    }

    opspec_t src = XMM(data->regs[i]-1);

    data->free_regs[data->top++] = data->regs[i];
    data->owner[data->regs[i]-1] = -1;

    data->regs[i] = data->spillpos--;

    opspec_t dest = MEM(RSP, 4*(data->regs[i]));
    plasm_put_op(data->as, VMOVSS, dest, src);
}

static opspec_t allocate_operand(codegen_data *data, int i) {
    if(data->top <= 0) {
        spill_operand(data);
    }
    data->regs[i] = data->free_regs[--data->top];
    data->owner[data->regs[i]-1] = i;
    return XMM(data->regs[i]-1);
}

static opspec_t calc_mem(codegen_data *data, int i) {
    uint64_t src = data->block->buffer[i];
    if(GET_OP(src) == SSA_LOAD) {
        if(GET_ARG1(src)<data->increment) {
            return MEM(RDI, 4*(GET_ARG1(src)));
        } else {
            return MEM(RSI, 4*(GET_ARG1(src)-data->increment));
        }
    } else if(GET_OP(src) == SSA_STORE) {
        if(GET_ARG2(src)<data->increment) {
            return MEM(RDI, 4*(GET_ARG2(src)));
        } else {
            return MEM(RSI, 4*(GET_ARG2(src)-data->increment));
        }
    } else {
        return MEM(RSP, 4*(data->regs[i]));
    }
}

static opspec_t lock_operand(codegen_data *data, int i, int mem) {
    if(data->regs[i]>0) {
        data->locked[data->regs[i]-1] = 1;
        return XMM(data->regs[i]-1);
    } else if(mem) {
        return calc_mem(data, i);
    } else {
        opspec_t source = calc_mem(data, i);
        opspec_t target = allocate_operand(data, i);
        plasm_put_op(data->as, VMOVSS, target, source);
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

void (*gencode_avx_ss(ssa_block *block, unsigned increment))(float*, float*, int)  {
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

    //~ for(size_t i = 0;i<block->size;++i) {
        //~ printf("%d %d %d ", (int)i, dependencies[i], uses[i]);
        //~ ssa_print_op(block, i);
        //~ printf("\t\t\t");
//~
        //~ for(size_t j = start[i];j<start[i+1];++j) {
            //~ printf(" %d", graph[j]);
        //~ }
        //~ printf("\n");
    //~ }

    int regs[block->size];
    int free_regs[16] = {16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    int locked[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int owner[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int top = 16;

    uint32_t chunksize = ((block->size*16u)+4095u)&(~4095u);
    binary_header *chunk = mmap(NULL, chunksize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    uint8_t *buffer = (uint8_t*)(chunk+1);
    chunk->chunksize = chunksize;

    plasm as;
    plasm_init(&as, buffer, chunksize-sizeof(binary_header));
    void (*fun)(float*, float*, int) = (void(*)(float*, float*, int))plasm_get_current_ptr(&as);

    codegen_data data = {0, regs, uses, free_regs, locked, owner, top, -1, increment, &as, block, start, graph};

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
        data.pos = i;
        regs[i] = 0;
        uint64_t op = block->buffer[i];
        if(GET_OP(op) == SSA_LOAD) {

        } else if(GET_OP(op) == SSA_STORE) {
            opspec_t source = lock_operand(&data, GET_ARG1(op), 0);

            plasm_put_op(&as, VMOVSS, calc_mem(&data, i), source);

            release_operand(&data, GET_ARG1(op));
        } else if(GET_OP(op) == SSA_CONST) {
            opspec_t target = allocate_operand(&data, i);

            *(constbase+constindex) = GET_CONST(op);

            plasm_put_op(&as, VMOVSS, target, MEM(RAX, 4*constindex));
            --constindex;

        } else if(IS_COMPARISON(op)) {
            int arg1 = GET_ARG1(op);
            int arg2 = GET_ARG2(op);
            opspec_t src1 = lock_operand(&data, arg1, 0);
            opspec_t src2 = lock_operand(&data, arg2, 1);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, VCMPSS, target, src1, src2, IMM8(cmpimm[GET_OP(op)]));

            release_operand(&data, arg2);
            release_operand(&data, arg1);
        } else if(IS_UNARY(op)) {
            int arg1 = GET_ARG1(op);
            opspec_t src1 = lock_operand(&data, arg1, 1);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, ssa2asm_ss[GET_OP(op)], target, target, src1);

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
            opspec_t src2 = lock_operand(&data, arg2, 0);
            opspec_t src3 = lock_operand(&data, arg3, 0);
            opspec_t target = allocate_operand(&data, i);

            plasm_put_op(&as, ssa2asm_ss[GET_OP(op)], target, src3, src2, src1);

            release_operand(&data, arg3);
            release_operand(&data, arg2);
            release_operand(&data, arg1);
        }
    }

    plasm_put_op(&as, ADD, RDI, IMM8(4*increment));
    plasm_put_op(&as, CMP, RDI, RDX);
    int32_t *reloff_ptr;
    plasm_put_op(&as, JL, VALUE_PTR(REL32(0), &reloff_ptr));
    *reloff_ptr = loop - plasm_get_current_ptr(&as);

    plasm_put_op(&as, POP, RBP);
    plasm_put_op(&as, RET);

    chunk->size = as.position;

    mprotect(chunk, chunksize, PROT_EXEC|PROT_READ);

    return fun;
}

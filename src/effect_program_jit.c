#include "effect_program_jit.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "jitlib/ssa_def.h"
#include "jitlib/ssa_scheduling.h"
#include "jitlib/ssa_parser.h"
#include "jitlib/ssa_print.h"
#include "jitlib/ssa_optimizer_passes.h"
#include "jitlib/codegen.h"

#include "effect_program.h"

uint64_t nanotime();

static struct {
    int parameters;
    int inputs;
    char *source;
} sources [] = {
    [EFFECT_TYPE_LINEAR_ACCEL]      = {4, 1,
        "dt = [8]\n"
        "[4] = [4] + [9]*dt\n"
        "[5] = [5] + [10]*dt\n"
        "[6] = [6] + [11]*dt\n"
        },

    [EFFECT_TYPE_LINEAR_FORCE]      = {4, 1,
        "dt = [8]\n"
        "[4] = [4] + [9]*(dt/[3])\n"
        "[5] = [5] + [10]*(dt/[3])\n"
        "[6] = [6] + [11]*(dt/[3])\n"},

    [EFFECT_TYPE_CENTRAL_FORCE]     = {5, 1,
        "dt = [8]\n"
        "x = [9]\n"
        "y = [10]\n"
        "z = [11]\n"
        "mu = [12]\n"
        "diffx = [0]-x\n"
        "diffy = [1]-y\n"
        "diffz = [2]-z\n"

        "r = 1.0/sqrt(diffx*diffx + diffy*diffy + diffz*diffz)\n"

        //~ "r0 = diffx*diffx + diffy*diffy + diffz*diffz\n"
        //~ "r1 = rsqrt(r0)\n"
        //~ "r = (0.5)*r1*(3-r0*r1*r1)\n"
        //~ "r3 = r*r*r\n"

        //~ "r0 = sqrt(diffx*diffx + diffy*diffy + diffz*diffz)\n"
        //~ "r1 = rcp(r0)\n"
        //~ "r = 2*r1-r0*r1*r1\n"

        "r3 = r*r*r\n"
        "[4] = [4] + dt*mu*diffx*r3\n"
        "[5] = [5] + dt*mu*diffy*r3\n"
        "[6] = [6] + dt*mu*diffz*r3\n"
        },

    [EFFECT_TYPE_PLANE_BOUNCE]      = {6, 1,
        "dt = [8]\n"
        "normalx = [9]\n"
        "normaly = [10]\n"
        "normalz = [11]\n"
        "d = [12]\n"
        "a = [13]\n"
        "dist = normalx*[0] + normaly*[1] + normalz*[2]\n"
        "vn = normalx*[4] + normaly*[5] + normalz*[6]\n"
        "mask = dist<d & vn<0\n"

        "dv = mask & vn\n"
        "factor = mask ? a:1\n"
        "[4] = ([4]-2*normalx*dv)*factor\n"
        "[5] = ([5]-2*normaly*dv)*factor\n"
        "[6] = ([6]-2*normalz*dv)*factor\n"
        "offset = mask & (d-dist)\n"
        "[0] = [0] + normalx*offset\n"
        "[1] = [1] + normaly*offset\n"
        "[2] = [2] + normalz*offset\n"
        },

    [EFFECT_TYPE_SPHERE_BOUNCE]     = {6, 1,
        "dt = [8]\n"
        "x0 = [9]\n"
        "y0 = [10]\n"
        "z0 = [11]\n"
        "r = [12]\n"
        "a = [13]\n"
        "diffx = [0]-x0\n"
        "diffy = [1]-y0\n"
        "diffz = [2]-z0\n"
        "dot = diffx*diffx + diffy*diffy + diffz*diffz\n"
        "dist = sqrt(dot)\n"
        "normalx = diffx*(1.0/dist)\n"
        "normaly = diffy*(1.0/dist)\n"
        "normalz = diffz*(1.0/dist)\n"
        "vx = [4]\n"
        "vy = [5]\n"
        "vz = [6]\n"
        "vn = normalx*vx + normaly*vy + normalz*vz\n"
        "mask = dist<r & vn<0\n"
        "dv = mask & vn\n"
        "factor = mask ? a:1\n"
        "[4] = ([4]-2*normalx*dv)*factor\n"
        "[5] = ([5]-2*normaly*dv)*factor\n"
        "[6] = ([6]-2*normalz*dv)*factor\n"
        "offset = mask & (r-dist)\n"
        "[0] = [0] + normalx*offset\n"
        "[1] = [1] + normaly*offset\n"
        "[2] = [2] + normalz*offset\n"
        },

    [EFFECT_TYPE_NEWTON_STEP]       = {1, 1,
        "dt = [8]\n"
        "[0] = [0] + dt*[4]\n"
        "[1] = [1] + dt*[5]\n"
        "[2] = [2] + dt*[6]\n"
        },

    [EFFECT_TYPE_GRAVITY_FORCE]     = {0, 2, ""},
    [EFFECT_TYPE_SPHERE_COLLISION]  = {3, 2,
        "dt = [16]\n"
        "r = [17]\n"
        "a = [18]\n"
        "x0 = [0]\n"
        "y0 = [1]\n"
        "z0 = [2]\n"
        "x1 = [8]\n"
        "y1 = [9]\n"
        "z1 = [10]\n"


        "[4] = (x1-x0)*dt\n"
        "[5] = (y1-y0)*dt\n"
        "[6] = (z1-z0)*dt\n"


        /*
        "diffx = x1-x0\n"
        "diffy = y1-y0\n"
        "diffz = z1-z0\n"
        "dot = diffx*diffx + diffy*diffy + diffz*diffz\n"
        "dist = sqrt(dot)\n"
        "normalx = diffx*(1.0/dist)\n"
        "normaly = diffy*(1.0/dist)\n"
        "normalz = diffz*(1.0/dist)\n"
        "[4] = [4] + normalx\n"
        "[5] = [5] + normaly\n"
        "[6] = [6] + normalz\n"
        "[12] = [12] - normalx\n"
        "[13] = [13] - normaly\n"
        "[14] = [14] - normalz\n"
        */
        },
};

typedef struct effect_jit_t {
    int parameters, inputs;
    ssa_block block;
} effect_jit;

typedef struct effect_program_jit_t {
    int level;
    effect_jit effects[EFFECT_TYPE_COUNT];
    performance_count pc;
    float *const_input;
    void (*fun)(float*, float*, int);
} effect_program_jit;

static void effect_program_jit_compile(effect_program *self, const effect_desc *desc);
static void effect_program_jit_execute(const effect_program *self, particle_array *arr, float dt);
static void effect_program_jit_destroy(effect_program *self);
static void effect_program_jit_perf (const effect_program *self, const particle_array *arr, float dt, performance_count *out);

int effect_program_create_jit(effect_program *p, int level) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_jit_compile;
    p->execute = effect_program_jit_execute;
    p->destroy = effect_program_jit_destroy;
    p->perf_c = effect_program_jit_perf;
    effect_program_jit *program = (p->usr = malloc(sizeof(effect_program_jit)));

    uint64_t t0 = nanotime();
    for(int i = 0;i<EFFECT_TYPE_COUNT;++i) {
        program->effects[i].parameters = sources[i].parameters;
        program->effects[i].inputs = sources[i].inputs;
        program->effects[i].block = ssa_parse(sources[i].source);
    }
    uint64_t t1 = nanotime();

    //~ fprintf(stderr, "parse:     %f\n", (t1-t0)*1.0e-6);
    program->level = level;
    program->const_input = NULL;
    program->fun = NULL;

    return 0;
}

void effect_program_jit_destroy(effect_program *self) {
    effect_program_jit *program = self->usr;
    for(int i = 0;i<EFFECT_TYPE_COUNT;++i) {
        ssa_block_destroy(&program->effects[i].block);
    }

    free_code(program->fun);
    free(program->const_input);
    free(program);
}

void count_ops(ssa_block *block, performance_count *pc) {
    memset(pc, 0, sizeof(performance_count));
    for(size_t i = 0;i<block->size;++i) {
        switch(GET_OP(block->buffer[i])) {
            case SSA_NONE: break;
            case SSA_LOAD: pc->loads++; break;
            case SSA_CONST: break;
            case SSA_STORE: pc->stores++; break;
            case SSA_ADD: pc->add++; break;
            case SSA_MUL: pc->mul++; break;
            case SSA_MIN: break;
            case SSA_MAX: break;
            case SSA_AND: break;
            case SSA_ANDN: break;
            case SSA_XOR: break;
            case SSA_OR: break;
            case SSA_DIV: pc->div++; break;
            case SSA_SUB: pc->add++; break;
            case SSA_SQRT: pc->sqrt++; break;
            case SSA_RSQRT: break;
            case SSA_RCP: break;
            case SSA_BLEND: break;
            case SSA_CMPEQ: pc->cmp++; break;
            case SSA_CMPLT: pc->cmp++; break;
            case SSA_CMPLE: pc->cmp++; break;
            case SSA_CMPNEQ: pc->cmp++; break;
            case SSA_CMPNLT: pc->cmp++; break;
            case SSA_CMPNLE: pc->cmp++; break;
        }
    }
}

void effect_program_jit_compile(effect_program *self, const effect_desc *desc) {
    effect_program_jit *program = self->usr;

    ssa_block block;
    ssa_block_load(&block, NULL, 0);

    program->const_input = realloc(program->const_input, sizeof(float)*128);

    int paramcount = 1;
    for(size_t i = 0;i<desc->size;++i) {
        uint64_t input_remap[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
        const effect_desc_ele *el = &desc->elements[i];

        int effect_params = program->effects[el->type].parameters;
        int effect_inputs = program->effects[el->type].inputs;

        for(int j = 1;j<effect_params;++j) {
            program->const_input[paramcount] = el->float_usr[j-1];
            input_remap[8*effect_inputs + j] = 8*effect_inputs + paramcount;
            paramcount++;
        }
        ssa_block_append(&block, &program->effects[el->type].block, input_remap);
    }


    uint64_t t0 = nanotime();
    if(program->level >= JIT_O1) {
        ssa_fuse_load_store(&block);
    }
    uint64_t t1 = nanotime();

    if(program->level >= JIT_O2) {
        ssa_fold_constants(&block);
    }
    uint64_t t2 = nanotime();

    if(program->level >= JIT_O2) {
        ssa_remap_duplicates(&block, 1024);
    }
    uint64_t t3 = nanotime();

    if(program->level >= JIT_O2) {
        ssa_mark_dead(&block);
    }
    uint64_t t4 = nanotime();

    ssa_resolve(&block);
    uint64_t t5 = nanotime();

    if(program->level >= JIT_O3) {
        switch(program->level&0x0F) {
            case JIT_AVX8: ssa_schedule(&block, ivybridge8); break;
            case JIT_AVX4: ssa_schedule(&block, ivybridge4); break;
            default:       ssa_schedule(&block, ivybridge1); break;
        }

    }
    uint64_t t6 = nanotime();

    count_ops(&block, &program->pc);

    switch(program->level&0x0F) {
        case JIT_AVX8: program->fun = gencode_avx8_ps(&block, 8); break;
        case JIT_AVX4: program->fun = gencode_avx_ps(&block, 8); break;
        default:
        case JIT_AVX1: program->fun = gencode_avx_ss(&block, 8); break;
    }
    uint64_t t7 = nanotime();
/*
    fprintf(stderr, "fuse:      %f\n", (t1-t0)*1.0e-6);
    fprintf(stderr, "fold:      %f\n", (t2-t1)*1.0e-6);
    fprintf(stderr, "dedup:     %f\n", (t3-t2)*1.0e-6);
    fprintf(stderr, "mark:      %f\n", (t4-t3)*1.0e-6);
    fprintf(stderr, "resolve:   %f\n", (t5-t4)*1.0e-6);
    fprintf(stderr, "schedule:  %f\n", (t6-t5)*1.0e-6);
    fprintf(stderr, "generate:  %f\n", (t7-t6)*1.0e-6);
*/

    //~ dump_code(program->fun);
    ssa_block_destroy(&block);
}

void effect_program_jit_execute(const effect_program *self, particle_array *arr, float dt) {
    effect_program_jit *program = self->usr;
    program->const_input[0] = dt;
    program->fun((float*)(arr->particles), program->const_input, 8*arr->size);
}

static void effect_program_jit_perf (const effect_program *self, const particle_array *arr, float dt, performance_count *out) {
    effect_program_jit *program = self->usr;
    memset(out, 0, sizeof(performance_count));
    size_t particles = particle_array_size(arr);

    out->add = program->pc.add * particles;
    out->cmp = program->pc.cmp * particles;
    out->mul = program->pc.mul * particles;
    out->div = program->pc.div * particles;
    out->sqrt = program->pc.sqrt * particles;
    out->loads = program->pc.loads * particles;
    out->stores = program->pc.stores * particles;
}

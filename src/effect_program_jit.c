#include "effect_program_jit.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <jitlib/ssa_def.h>
#include <jitlib/ssa_scheduling.h>
#include <jitlib/ssa_parser.h>
#include <jitlib/ssa_print.h>
#include <jitlib/ssa_optimizer_passes.h>
#include <jitlib/codegen.h>

uint64_t nanotime();

static struct {
    int parameters;
    char *source;
} sources [] = {
    [EFFECT_TYPE_LINEAR_ACCEL]      = {4,
        "dt = [8]\n"
        "[4] = [4] + [9]*dt\n"
        "[5] = [5] + [10]*dt\n"
        "[6] = [6] + [11]*dt\n"
        },

    [EFFECT_TYPE_LINEAR_FORCE]      = {4,
        "dt = [8]\n"
        "[4] = [4] + [9]*(dt/[3])\n"
        "[5] = [5] + [10]*(dt/[3])\n"
        "[6] = [6] + [11]*(dt/[3])\n"},

    [EFFECT_TYPE_CENTRAL_FORCE]     = {5,
        "dt = [8]\n"
        "x = [9]\n"
        "y = [10]\n"
        "z = [11]\n"
        "mu = [12]\n"
        "diffx = [0]-x\n"
        "diffy = [1]-y\n"
        "diffz = [2]-z\n"
        "r = 1.0/sqrt(diffx*diffx + diffy*diffy + diffz*diffz)\n"
        "r3 = r*r*r\n"
        "[4] = [4] + dt*mu*diffx*r3\n"
        "[5] = [5] + dt*mu*diffy*r3\n"
        "[6] = [6] + dt*mu*diffz*r3\n"
        },

    [EFFECT_TYPE_PLANE_BOUNCE]      = {6,
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

    [EFFECT_TYPE_SPHERE_BOUNCE]     = {6,
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

    [EFFECT_TYPE_NEWTON_STEP]       = {1,
        "dt = [8]\n"
        "[0] = [0] + dt*[4]\n"
        "[1] = [1] + dt*[5]\n"
        "[2] = [2] + dt*[6]\n"
        },

    [EFFECT_TYPE_GRAVITY_FORCE]     = {0, ""},
    [EFFECT_TYPE_SPHERE_COLLISION]  = {0, ""},
};

typedef struct effect_jit_t {
    int parameters;
    ssa_block block;
} effect_jit;

typedef struct effect_program_jit_t {
    effect_jit effects[EFFECT_TYPE_COUNT];
    float *const_input;
    void (*fun)(float*, float*, int);
} effect_program_jit;

void effect_program_jit_compile(effect_program *self, const effect_desc *desc);
void effect_program_jit_execute(const effect_program *self, particle_array *arr, float dt);
void effect_program_jit_destroy(effect_program *self);


int effect_program_create_jit(effect_program *p) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_jit_compile;
    p->execute = effect_program_jit_execute;
    p->destroy = effect_program_jit_destroy;
    effect_program_jit *program = (p->usr = malloc(sizeof(effect_program_jit)));

    uint64_t t0 = nanotime();
    for(int i = 0;i<EFFECT_TYPE_COUNT;++i) {
        program->effects[i].parameters = sources[i].parameters;
        program->effects[i].block = ssa_parse(sources[i].source);
    }
    uint64_t t1 = nanotime();

    fprintf(stderr, "parse:     %f\n", (t1-t0)*1.0e-6);
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

void effect_program_jit_compile(effect_program *self, const effect_desc *desc) {
    effect_program_jit *program = self->usr;

    ssa_block block;
    ssa_block_load(&block, NULL, 0);

    program->const_input = realloc(program->const_input, sizeof(float)*128);

    uint64_t input_remap[16] = {0,1,2,3,4,5,6,7,8};
    int paramcount = 1;
    for(size_t i = 0;i<desc->size;++i) {
        const effect_desc_ele *el = &desc->elements[i];

        int effect_params = program->effects[el->type].parameters;

        for(int j = 1;j<effect_params;++j) {
            program->const_input[paramcount] = el->float_usr[j-1];
            input_remap[8 + j] = 8 + paramcount;
            paramcount++;
        }
        ssa_block_append(&block, &program->effects[el->type].block, input_remap);
    }


    uint64_t t0 = nanotime();
    ssa_fuse_load_store(&block);
    uint64_t t1 = nanotime();

    ssa_fold_constants(&block);
    uint64_t t2 = nanotime();

    ssa_remap_duplicates(&block, 1024);
    uint64_t t3 = nanotime();

    ssa_mark_dead(&block);
    uint64_t t4 = nanotime();

    ssa_resolve(&block);
    uint64_t t5 = nanotime();

    ssa_schedule(&block, ivybridge);
    uint64_t t6 = nanotime();

    program->fun = gencode_avx_ss(&block, 8);
    uint64_t t7 = nanotime();

    fprintf(stderr, "fuse:      %f\n", (t1-t0)*1.0e-6);
    fprintf(stderr, "fold:      %f\n", (t2-t1)*1.0e-6);
    fprintf(stderr, "dedup:     %f\n", (t3-t2)*1.0e-6);
    fprintf(stderr, "mark:      %f\n", (t4-t3)*1.0e-6);
    fprintf(stderr, "resolve:   %f\n", (t5-t4)*1.0e-6);
    fprintf(stderr, "schedule:  %f\n", (t6-t5)*1.0e-6);
    fprintf(stderr, "generate:  %f\n", (t7-t6)*1.0e-6);

    //~ dump_code(program->fun);

    ssa_block_destroy(&block);
}

void effect_program_jit_execute(const effect_program *self, particle_array *arr, float dt) {
    effect_program_jit *program = self->usr;
    program->const_input[0] = dt;
    program->fun((float*)(arr->particles), program->const_input, 8*arr->size);
}


#ifndef __effect_program_H
#define __effect_program_H

#include <stdint.h>
#include <stddef.h>
#include "effect_desc.h"
#include "particle_array.h"

struct effect_program_t;
typedef struct effect_program_t effect_program;

typedef struct performance_count_t {
    size_t add;
    size_t cmp;
    size_t mul;
    size_t div;
    size_t sqrt;
    size_t loads;
    size_t stores;
} performance_count;

struct effect_program_t {
    void (*compile)(      effect_program *self,    const effect_desc *desc);
    void (*execute)(const effect_program *self,       particle_array *arr, float dt);
    void (*destroy)(      effect_program *self);
    void (*perf_c) (const effect_program *self, const particle_array *arr, float dt, performance_count *out);
    void *usr;
};

void effect_program_compile(effect_program *p, const effect_desc *desc);

void effect_program_perf_c (const effect_program *p, const particle_array *arr, float dt, performance_count *out);
void effect_program_execute(const effect_program *p, particle_array *arr, float dt);

void effect_program_destroy(effect_program *p);

#endif // __effect_program_H

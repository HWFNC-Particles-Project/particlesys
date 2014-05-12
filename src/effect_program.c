#include "effect_program.h"


void effect_program_compile(effect_program *p, const effect_desc *desc) {
    if (p->compile) p->compile(p, desc);
}

void effect_program_perf_c (const effect_program *p, const particle_array *arr, float dt, performance_count *out) {
    if (p->perf_c) p->perf_c(p, arr, dt, out);
}

void effect_program_execute(const effect_program *p, particle_array *arr, float dt) {
    if (p->execute) p->execute(p, arr, dt);
}

void effect_program_destroy(effect_program *p) {
    if (p->destroy) p->destroy(p);
}




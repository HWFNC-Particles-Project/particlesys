#include <stdio.h>
#include <stdlib.h>
#include "particle_array.h"
#include "effect_program_naive.h"

int main(int argc, char *argv[]) {
    (void) argc; (void) argv;
    
    particle_array arr;
    particle_array_create(&arr);
    for(size_t i = 0;i<128;++i) {
        particle p;
        p.position[0] = i; p.position[1] = 0; p.position[2] = 0;
        p.velocity[0] = 0; p.velocity[1] = 0; p.velocity[2] = 0;
        p.mass = 1;
        p.charge = 1;
        particle_array_add(&arr, p);
    }
    
    effect_desc effects;
    effect_desc_init(&effects);
    // construct effect array:
    effect_desc_add_linear_accel(&effects, 0, -10, 0);
    effect_desc_add_newton_step (&effects);
    // create program:
    effect_program naive_program;
    effect_program_create_naive(&naive_program);
    // compile
    effect_program_compile(&naive_program, &effects);
    // execute
    effect_program_execute(&naive_program, &arr, 0.01);
    
    for(size_t i = 0; i < arr.size; ++i) {
        printf("%.3f, %.3f, %.3f\n", arr.particles[i].position[0], arr.particles[i].position[1], arr.particles[i].position[2]);
    }
    
    // clean up:
    effect_desc_destroy(&effects);
    effect_program_destroy(&naive_program);
    particle_array_destroy(&arr);
    
    return 0;
}

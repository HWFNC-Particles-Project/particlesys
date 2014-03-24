#include <stdio.h>
#include <stdlib.h>
#include <particle_array.h>
#include <particle_effects.h>

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
    
    particle_effect effects[] = {
        linear_force_effect(0,-10,0),
        newton_step_effect(),
        {NULL, NULL}
    };
    particle_array_apply_effects(&arr, effects, 0.01);
    particle_effect_free(effects);
    
    particle_array_destroy(&arr);
    
    return 0;
}

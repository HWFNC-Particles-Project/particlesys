#include <stdio.h>
#include <stdlib.h>
#include "test_effects.h"
#include "particle_array.h"
#include "effect_desc.h"
#include "performance_measurement.h"

#include <time.h>
#include <string.h>

#define DT 0.01f

static float randf(float min, float max) {
    return ((float)rand() / (float)RAND_MAX) * (max - min) + min;
}

static void add_random_particles(particle_array *arr, size_t N) {
    for(size_t i = 0;i<N;++i) {
        particle p;
        p.position[0] = randf(-1, 1); p.position[1] = randf(-1, 1); p.position[2] = randf(-1, 1);
        p.velocity[0] = randf(-1, 1); p.velocity[1] = randf(-1, 1); p.velocity[2] = randf(-1, 1);
        p.mass = randf(0.5, 1.5);
        p.charge = randf(-1, 1);
        particle_array_add(arr, p);
    }
}

static void add_random_particles_cond(particle_array *arr, size_t N, int (*pc)(particle, int, float *), int test_case, float *cond_data) {
    for(size_t i = 0;i<N;++i) {
        particle p;
        do {
            p.position[0] = randf(-1, 1); p.position[1] = randf(-1, 1); p.position[2] = randf(-1, 1);
            p.velocity[0] = randf(-1, 1); p.velocity[1] = randf(-1, 1); p.velocity[2] = randf(-1, 1);
            p.mass = randf(0.5, 1.5);
            p.charge = randf(-1, 1);
        } while(pc(p, test_case, cond_data) != 0);
        particle_array_add(arr, p);
    }
}

int particle_condition_function(particle p, int test_case, float *data) {
    switch(test_case) {
        case 3:
        case 4: {
                // plane bounce.
                float dist    = data[0]*p.position[0] + data[1]*p.position[1] + data[2]*p.position[2];
                float vnormal = data[0]*p.velocity[0] + data[1]*p.velocity[1] + data[2]*p.velocity[2];
                if (test_case == 3) { 
                    // all particles in front:
                    if(dist<data[3]) {
                        return -1;
                    } else {
                        return 0;
                    }
                } else if (test_case == 4) {
                    // all particles in back and moving to the back:
                    if(dist<data[3] && vnormal < 0.0f) {
                        return 0;
                    } else {
                        return -1;
                    }
                }
            }
            break;
    }
}

int test_effects_prepare(int test_case, effect_desc *effects, particle_array *arr) {
    effect_desc_init(effects);
    switch(test_case) {
        case 0:
            // linear accel:
            effect_desc_add_linear_accel(effects,  randf(-1,1), randf(-1,1), randf(-1,1));
            particle_array_create(arr);
            add_random_particles(arr, 1024);
            return 1;
            break;
        case 1:
            // linear force:
            effect_desc_add_linear_force(effects,  randf(-1,1), randf(-1,1), randf(-1,1));
            particle_array_create(arr);
            add_random_particles(arr, 1024);
            return 1;
            break;
        case 2:
            // central force:
            effect_desc_add_central_force(effects,  randf(-1,1), randf(-1,1), randf(-1,1), randf(-1,1));
            particle_array_create(arr);
            add_random_particles(arr, 1024);
            return 1;
            break;
        case 3:
        case 4: {
                float plane_data[4] = {randf(-1,1), randf(-1,1), randf(-1,1), randf(-1,1)};
                // plane bounce, all particles in front:
                effect_desc_add_plane_bounce(effects, plane_data[0], plane_data[1], plane_data[2], plane_data[3], randf(0,1));
                particle_array_create(arr);
                add_random_particles_cond(arr, 1024, particle_condition_function, test_case, plane_data);
                return 0;
            }
            break;
    }
    return -1;
}

void test_effects_all(effect_program *test_program) {
    // first, determine approximate number of iterations to perform
    int test_case = 0;
    while(1) {
        // allocate performance array:
        particle_array initial_arr;
        effect_desc effects;
        int tep_ret = test_effects_prepare(test_case, &effects, &initial_arr);
        if (tep_ret < 0) break;
        if (tep_ret > 0) {
            test_case++;
            continue;
        }
        printf("TESTCASE %d:\n", test_case);
        
        // compile:
        effect_program_compile(test_program, &effects);
        
        performance_count perf;
        
        particle_array arr;
        particle_array_create(&arr);
        // determine performance:
        particle_array_copy(&arr, &initial_arr);
        float dt = DT;
        // get performance:
        effect_program_perf_c(test_program, &arr, dt, &perf);
        
        double freq = 3.4e9;
        size_t cost = perf.add * 3 + perf.mul * 5 + perf.div * 11 + perf.cmp * 3 + perf.sqrt * 16 + perf.loads * 4 + perf.stores * 4;
        size_t repeats = 5 * freq / cost;
        
        // run it:
        clock_t start = clock();
        for (size_t r = 0; r < repeats; ++r) {
            particle_array_copy(&arr, &initial_arr);
            effect_program_execute(test_program, &arr, dt);
        }
        clock_t end = clock();
        printf("time: %d\n", (int)(end-start));
        double cycles = (end-start) * 0.001 * freq;
        
        performance_count total = perf;
        double part_iteration = (double)particle_array_size(&initial_arr);
        printf("cycles:      %8.2f\n", cycles / part_iteration / (double)repeats);
        //printf("cycles/add:  %8.2f\n", cycles / (double)total.add);
        //printf("cycles/cmp:  %8.2f\n", cycles / (double)total.cmp);
        //printf("cycles/mul:  %8.2f\n", cycles / (double)total.mul);
        //printf("cycles/div:  %8.2f\n", cycles / (double)total.div);
        //printf("cycles/sqrt: %8.2f\n", cycles / (double)total.sqrt);
        //printf("cycles/lds:  %8.2f\n", cycles / (double)total.loads);
        //printf("cycles/sts:  %8.2f\n", cycles / (double)total.stores);
        printf("add:         %8.2f\n", (double)total.add / part_iteration);
        printf("cmp:         %8.2f\n", (double)total.cmp / part_iteration);
        printf("mul:         %8.2f\n", (double)total.mul / part_iteration);
        printf("div:         %8.2f\n", (double)total.div / part_iteration);
        printf("sqrt:        %8.2f\n", (double)total.sqrt / part_iteration);
        printf("lds:         %8.2f\n", (double)total.loads / part_iteration);
        printf("sts:         %8.2f\n", (double)total.stores / part_iteration);
        
        particle_array_destroy(&arr);
        effect_desc_destroy(&effects);
        particle_array_destroy(&initial_arr);
        
        test_case++;
    }
}

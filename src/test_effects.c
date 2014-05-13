#include <stdio.h>
#include <stdlib.h>
#include "test_effects.h"
#include "particle_array.h"
#include "effect_desc.h"
#include "effect_program_naive.h"
#include "performance_measurement.h"

#include <time.h>
#include <string.h>
#include <math.h>

#define DT 0.01f

void verify(const effect_program *test_program, const effect_program *ref_program, const particle_array *initial_arr);

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
        case 2: {
                // avoid particles too near to the central force point:
                float d0 = p.position[0] - data[0];
                float d1 = p.position[1] - data[1];
                float d2 = p.position[2] - data[2];
                if (d0 * d0 + d1 * d1 + d2 * d2 < 0.001f * 0.001f) {
                    return -1;
                } else {
                    return 0;
                }
            }
            break;
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
        case 5:
        case 6: {
                // sphere bounce.
                float normal[3];
                normal[0] = p.position[0]-data[0];
                normal[1] = p.position[1]-data[1];
                normal[2] = p.position[2]-data[2];
                float r = sqrtf(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
                normal[0] /= r;
                normal[1] /= r;
                normal[2] /= r;
                float vnormal = normal[0]*p.velocity[0] + normal[1]*p.velocity[1] + normal[2]*p.velocity[2];
                float d = data[3];
                if (test_case == 5) { 
                    // all particles in front:
                    if(r<d) {
                        return -1;
                    } else {
                        return 0;
                    }
                } else if (test_case == 6) {
                    // all particles in back and moving to the back:
                    if(r<d && vnormal<0.0f) {
                        return 0;
                    } else {
                        return -1;
                    }
                }
            }
            break;
    }
    return 0;
}

int test_effects_prepare(int test_case, effect_desc *effects, particle_array *arr, const char **out_description) {
    effect_desc_init(effects);
    const char *dummy_p = NULL;
    if (out_description == NULL) {
        out_description = &dummy_p;
    }
    *out_description = "none";
    switch(test_case) {
        case 0:
            // linear accel:
            particle_array_create(arr);
            //return 1;
            *out_description = "linear acceleration";
            effect_desc_add_linear_accel(effects,  randf(-1,1), randf(-1,1), randf(-1,1));
            add_random_particles(arr, 1024);
            return 0;
            break;
        case 1:
            // linear force:
            particle_array_create(arr);
            //return 1;
            *out_description = "linear force";
            effect_desc_add_linear_force(effects,  randf(-1,1), randf(-1,1), randf(-1,1));
            add_random_particles(arr, 1024);
            return 0;
            break;
        case 2: {
                // central force:
                particle_array_create(arr);
                //return 1;
                *out_description = "central force";
                float cf_data[4] = {randf(-1,1), randf(-1,1), randf(-1,1)};
                effect_desc_add_central_force(effects,  cf_data[0], cf_data[1], cf_data[2], randf(0,1));
                add_random_particles_cond(arr, 1024, particle_condition_function, test_case, cf_data);
                return 0;
            }
            break;
        case 3:
        case 4: {
                particle_array_create(arr);
                return 1;
                if (test_case == 3) {
                    *out_description = "plane bounce, all particles in front";
                } else if (test_case == 4) {
                    *out_description = "plane bounce, all particles in the back, moving away from back";
                }
                float plane_data[4] = {randf(-1,1), randf(-1,1), randf(-1,1), randf(-1,1)};
                // plane bounce, all particles in front:
                effect_desc_add_plane_bounce(effects, plane_data[0], plane_data[1], plane_data[2], plane_data[3], randf(0,1));
                add_random_particles_cond(arr, 1024, particle_condition_function, test_case, plane_data);
                return 0;
            }
            break;
        case 5:
        case 6: {
                particle_array_create(arr);
                return 1;
                if (test_case == 5) {
                    *out_description = "sphere bounce, all particles outside";
                } else if (test_case == 6) {
                    *out_description = "sphere bounce, all particles inside, moving to the center";
                }
                float sphere_data[4] = {randf(-1,1), randf(-1,1), randf(-1,1), randf(0,1)};
                // sphere bounce, all particles in front:
                effect_desc_add_sphere_bounce(effects, sphere_data[0], sphere_data[1], sphere_data[2], sphere_data[3], randf(0,1));
                add_random_particles_cond(arr, 1024, particle_condition_function, test_case, sphere_data);
                return 0;
            }
            break;
        case 7:
            particle_array_create(arr);
            return 1;
            // pairwise gravity:
            *out_description = "pairwise gravity";
            effect_desc_add_gravity_force(effects,  randf(0,1));
            add_random_particles(arr, 200);
            return 0;
            break;
        case 8:
        case 9:
            particle_array_create(arr);
            return 1;
            // pairwise collision:
            *out_description = "pairwise collision";
            effect_desc_add_sphere_collision(effects, randf(0,0.1), randf(0,1));
            add_random_particles(arr, 200);
            return 0;
            break;
        case 10:
            particle_array_create(arr);
            return 1;
            // linear force:
            *out_description = "newton step";
            effect_desc_add_newton_step(effects);
            add_random_particles(arr, 1024);
            return 0;
            break;
    }
    return -1;
}

void test_effects_all(effect_program *test_program) {
    // first, determine approximate number of iterations to perform
    int test_case = 0;
    effect_program ref_program;
    effect_program_create_naive(&ref_program);
    while(1) {
        // allocate performance array:
        particle_array initial_arr;
        effect_desc effects;
        const char *description_text = NULL;
        int tep_ret = test_effects_prepare(test_case, &effects, &initial_arr, &description_text);
        if (tep_ret < 0) break;
        if (tep_ret == 0) {
            printf("TESTCASE %d: %s\n", test_case, description_text);
            
            // compile:
            effect_program_compile(test_program, &effects);
            effect_program_compile(&ref_program, &effects);
            
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
            size_t repeats = 10 * freq / cost;
            
            // verify:
            verify(test_program, &ref_program, &arr);
            
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
        }
        
        effect_desc_destroy(&effects);
        particle_array_destroy(&initial_arr);
        
        test_case++;
    }
}

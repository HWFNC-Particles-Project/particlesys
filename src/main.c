#include <stdio.h>
#include <stdlib.h>
#include "particle_array.h"
#include "effect_program_naive.h"
#include "effect_program_c.h"
#include "effect_program_jit.h"
#include "particle_vis.h"
#include "performance_measurement.h"
#include "test_effects.h"

#include <time.h>
#include <string.h>

#define LOGGING 1
#define LOGFILE "perf_log"

#define DT 0.01f

void test_performance(const effect_program *test_program, const particle_array *initial_arr) {
    // first, determine approximate number of iterations to perform
    size_t iterations = 100;
    size_t repeats = 1;
    // allocate performance array:
    performance_count *perf_array = malloc(iterations * sizeof(performance_count));
    particle_array arr;
    particle_array_create(&arr);
    // determine performance at every timestep:
    particle_array_copy(&arr, initial_arr);
    float dt = DT;
    for (size_t i = 0; i < iterations; ++i) {
        // get performance:
        effect_program_perf_c(test_program, &arr, dt, &perf_array[i]);
        // step to next timestep:
        effect_program_execute(test_program, &arr, dt);
        //if(i%50==0){printf("%d-th iteration\n",i);}
    }
    // run it:
	clock_t start = clock();
	for (size_t r = 0; r < repeats; ++r) {
        particle_array_copy(&arr, initial_arr);
        for (size_t i = 0; i < iterations; ++i) {
            effect_program_execute(test_program, &arr, dt);
            //if(i%50==0){printf("%d-th iteration\n",i);}
        }
	}
    clock_t end = clock();
    printf("time: %d\n", (int)(end-start));
    double cycles = (end-start) * 0.001 * 3.4e9;
    
    performance_count total;
    memset(&total, 0, sizeof(total));
    for (size_t i = 0; i < iterations; ++i) {
        total.add += perf_array[i].add * repeats;
        total.cmp += perf_array[i].cmp * repeats;
        total.mul += perf_array[i].mul * repeats;
        total.div += perf_array[i].div * repeats;
        total.sqrt += perf_array[i].sqrt * repeats;
        total.loads += perf_array[i].loads * repeats;
        total.stores += perf_array[i].stores * repeats;
    }
    double part_iteration = (double)particle_array_size(initial_arr) * (double)iterations * (double)repeats;
    printf("cycles:      %8.2f\n", cycles / part_iteration);
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
    free(perf_array);
}

void verify(const effect_program *test_program, const effect_program *ref_program, const particle_array *initial_arr) {
    int iterations = 50;
    particle_array test_arr;
    particle_array ref_arr;
    particle_array_create(&test_arr);
    particle_array_create(&ref_arr);
    particle_array_copy(&test_arr, initial_arr);
    float dt = DT;
    for(int i = 0;i<iterations;++i) {
        // at each iteration copy test arr to ref arr and compare errors only after one timestep:
        particle_array_copy(&ref_arr, &test_arr);
        effect_program_execute(ref_program,  &ref_arr,  dt);
        effect_program_execute(test_program, &test_arr, dt);
        particle_vis_draw(&test_arr);
        if (particle_array_compare(&ref_arr, &test_arr, 1e-4) != 0) {
            fprintf(stderr, "iteration %d\n", (int)i);
            break;
        }
        //if(i%50==0){printf("%d-th iteration\n",i);}
    }
    particle_array_destroy(&test_arr);
    particle_array_destroy(&ref_arr);
}

int main(int argc, char *argv[]) {
    (void) argc; (void) argv;

	FILE* logfile;
	logfile = fopen(LOGFILE, "w");

    perf_measure perf_default = {"default measurement", 0, 0, LOGGING, logfile};
    perf_measure perf_program_creation = perf_default;
    perf_program_creation.name = "creation";
    perf_measure perf_program_execution = perf_default;
    perf_program_execution.name = "execution";
    perf_measurement_init();

    particle_vis_init();

    particle_array initial_arr;
    particle_array_create(&initial_arr);
    unsigned int N = 32;
    for(size_t i = 0;i<N;++i) {
        for(size_t j = 0;j<N;++j) {
            particle p;
            p.position[0] = (i+0.5f)/(float)N-0.5f; p.position[1] = 1; p.position[2] = (j+0.5f)/(float)N-0.5f;
            p.velocity[0] = 0; p.velocity[1] = 0; p.velocity[2] = 0;
            p.mass = 1;
            p.charge = 1;
            particle_array_add(&initial_arr, p);
        }
    }

    effect_desc effects;
    effect_desc_init(&effects);
    // construct effect array:
    effect_desc_add_linear_accel(&effects,  0, -10, 0);

    effect_desc_add_sphere_collision(&effects,  0.015, 0.7);

    //effect_desc_add_sphere_bounce(&effects, 0,  0,  0, 0.84, 0.8);
    //effect_desc_add_plane_bounce(&effects,  0,  1,  0, 0.8,    0.8);
    //effect_desc_add_plane_bounce(&effects,  0,  -1,  0, -0.85,    0.8);
    effect_desc_add_sphere_bounce(&effects, 0,  0,  0, 0.707, 0.8);
    effect_desc_add_plane_bounce(&effects,  0,  1,  0, -1,    0.8);
    effect_desc_add_plane_bounce(&effects,  0, -1,  0, -1,    0.8);
    effect_desc_add_plane_bounce(&effects,  0, -1,  0, -1,    0.8);
    effect_desc_add_plane_bounce(&effects,  1,  0,  0, -1,    0.8);
    effect_desc_add_plane_bounce(&effects, -1,  0,  0, -1,    0.8);
    effect_desc_add_plane_bounce(&effects,  0,  0,  1, -1,    0.8);
    effect_desc_add_plane_bounce(&effects,  0,  0, -1, -1,    0.8);

//    effect_desc_add_central_force(&effects, 1, 0, 1, -1);
//    effect_desc_add_central_force(&effects,-1, 0, 1, -1);
//    effect_desc_add_central_force(&effects, 1, 0,-1, -1);
//    effect_desc_add_central_force(&effects,-1, 0,-1, -1);
//
//    effect_desc_add_sphere_bounce(&effects,  1,  0,  1, 0.5, 0.9);
//    effect_desc_add_sphere_bounce(&effects,  1,  0, -1, 0.5, 0.9);
//    effect_desc_add_sphere_bounce(&effects, -1,  0,  1, 0.5, 0.9);
//    effect_desc_add_sphere_bounce(&effects, -1,  0, -1, 0.5, 0.9);


    effect_desc_add_newton_step (&effects);

    // create reference program:
    effect_program ref_program;
    effect_program_create_naive(&ref_program);
    effect_program_compile(&ref_program, &effects);

    // create program:

    perf_start_measurement(&perf_program_creation);

    effect_program test_program_1;
    effect_program test_program_0;

    effect_program_create_naive(&test_program_0);
    effect_program_create_c_optimze1(&test_program_1);
    //effect_program_create_jit(&test_program);

    test_effects_all(&test_program_1);
    test_effects_all(&test_program_0);

    // compile
    effect_program_compile(&test_program_0, &effects);
    effect_program_compile(&test_program_1, &effects);

    perf_stop_measurement(&perf_program_creation);

    // execute
	perf_start_measurement(&perf_program_execution);
	
	//test_performance(&test_program_1, &initial_arr);
	//test_performance(&test_program_0, &initial_arr);
    //verify(&test_program_1, &ref_program, &initial_arr);
    
    perf_stop_measurement(&perf_program_execution);

    // print mesurements
    perf_print_measurement(&perf_program_creation);
    perf_print_measurement(&perf_program_execution);
    // clean up:
    effect_desc_destroy(&effects);
    effect_program_destroy(&ref_program);
    effect_program_destroy(&test_program_0);
    effect_program_destroy(&test_program_1);
    particle_array_destroy(&initial_arr);
    particle_vis_deinit();
	fclose(logfile);
    return 0;
}

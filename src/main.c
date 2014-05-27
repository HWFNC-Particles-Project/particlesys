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

#define DT 0.001f

uint64_t nanotime();

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
    uint64_t t_total = 0;
	for (size_t r = 0; r < repeats; ++r) {
        particle_array_copy(&arr, initial_arr);
        uint64_t start = nanotime();
        for (size_t i = 0; i < iterations; ++i) {
            effect_program_execute(test_program, &arr, dt);
            //if(i%50==0){printf("%d-th iteration\n",i);}
        }
        uint64_t end = nanotime();
        t_total += end-start;
	}
    //~ printf("time: %d\n", (int)(t_total));
    double cycles = (t_total) * 3.2;

    performance_count total;
    memset(&total, 0, sizeof(total));
    for (size_t i = 0; i < iterations; ++i) {
        total.add += perf_array[i].add * repeats;
        total.cmp += perf_array[i].cmp * repeats;
        total.mul += perf_array[i].mul * repeats;
        total.div += perf_array[i].div * repeats;
        total.rcp += perf_array[i].rcp * repeats;
        total.sqrt += perf_array[i].sqrt * repeats;
        total.rsqrt += perf_array[i].rsqrt * repeats;
        total.loads += perf_array[i].loads * repeats;
        total.stores += perf_array[i].stores * repeats;
    }
    double part_iteration = (double)particle_array_size(initial_arr) * (double)iterations * (double)repeats;
    //~ printf("%8.2f ", cycles / part_iteration);
    printf("cycles:      %8.2f\n", cycles / part_iteration);
    printf("add:         %8.2f\n", (double)total.add / part_iteration);
    printf("cmp:         %8.2f\n", (double)total.cmp / part_iteration);
    printf("mul:         %8.2f\n", (double)total.mul / part_iteration);
    printf("div:         %8.2f\n", (double)total.div / part_iteration);
    printf("rcp:         %8.2f\n", (double)total.rcp / part_iteration);
    printf("sqrt:        %8.2f\n", (double)total.sqrt / part_iteration);
    printf("rsqrt:       %8.2f\n", (double)total.rsqrt / part_iteration);
    printf("lds:         %8.2f\n", (double)total.loads / part_iteration);
    printf("sts:         %8.2f\n", (double)total.stores / part_iteration);
    printf("flops/cycle: %8.2f\n", (total.add+total.cmp+total.mul+total.div+total.sqrt+total.rsqrt+total.rcp)/(double)cycles);

    particle_array_destroy(&arr);
    free(perf_array);
}

void print_array(particle_array *arr) {
    for(size_t i = 0;i<particle_array_size(arr);++i) {
        particle p = particle_array_get(arr, i);
        printf("%f %f %f %f %f %f %f %f\n", p.array[0], p.array[1], p.array[2], p.array[3], p.array[4], p.array[5], p.array[6], p.array[7]);
    }
}

void verify(const effect_program *test_program, const effect_program *ref_program, const particle_array *initial_arr) {
    int iterations = 5000;
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

        if (particle_array_compare(&ref_arr, &test_arr, 1e-5, 1e-5) != 0) {
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
    unsigned int N = 16;
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
    //~ effect_desc_add_linear_accel(&effects,  0, -10, 0);
    //~ effect_desc_add_sphere_bounce(&effects, 0,  0,  0, 0.707, 0.8);
    //~ effect_desc_add_plane_bounce(&effects,  0,  1,  0, -1,    0.8);
    //~ effect_desc_add_plane_bounce(&effects,  0, -1,  0, -1,    0.8);
    //~ effect_desc_add_plane_bounce(&effects,  1,  0,  0, -1,    0.8);
    //~ effect_desc_add_plane_bounce(&effects, -1,  0,  0, -1,    0.8);
    //~ effect_desc_add_plane_bounce(&effects,  0,  0,  1, -1,    0.8);
    //~ effect_desc_add_plane_bounce(&effects,  0,  0, -1, -1,    0.8);
    //~ effect_desc_add_central_force(&effects, 1, 0, 1, -1);
    //~ effect_desc_add_central_force(&effects,-1, 0, 1, -1);
    //~ effect_desc_add_central_force(&effects, 1, 0,-1, -1);
    //~ effect_desc_add_central_force(&effects,-1, 0,-1, -1);
    //~ effect_desc_add_sphere_bounce(&effects,  1,  0,  1, 0.1, 0.9);
    //~ effect_desc_add_sphere_bounce(&effects,  1,  0, -1, 0.1, 0.9);
    //~ effect_desc_add_sphere_bounce(&effects, -1,  0,  1, 0.1, 0.9);
    //~ effect_desc_add_sphere_bounce(&effects, -1,  0, -1, 0.1, 0.9);
    effect_desc_add_central_force(&effects, 1, 0, 1, -1);
    effect_desc_add_sphere_bounce(&effects,  1,  0,  1, 0.1, 0.9);
    effect_desc_add_central_force(&effects,-1, 0, 1, -1);
    effect_desc_add_sphere_bounce(&effects, -1,  0,  1, 0.1, 0.9);
    effect_desc_add_central_force(&effects, 1, 0,-1, -1);
    effect_desc_add_sphere_bounce(&effects,  1,  0, -1, 0.1, 0.9);
    effect_desc_add_central_force(&effects,-1, 0,-1, -1);
    effect_desc_add_sphere_bounce(&effects, -1,  0, -1, 0.1, 0.9);



    effect_desc_add_newton_step (&effects);

    // create reference program:
    effect_program ref_program;
    effect_program_create_naive(&ref_program);
    effect_program_compile(&ref_program, &effects);

    // create program:

    perf_start_measurement(&perf_program_creation);

    effect_program test_program_0;
    effect_program test_program_1;
    effect_program test_program_2;
    effect_program test_program_3;
    effect_program jit_program_0;
    effect_program jit_program_1;
    effect_program jit_program_2;
    effect_program jit_program_3;

    effect_program_create_naive(&test_program_0);
    effect_program_create_c_optimze1(&test_program_1);
    effect_program_create_c_optimze2(&test_program_2);
    effect_program_create_jit(&jit_program_0, JIT_AVX4 | JIT_O0);
    effect_program_create_jit(&jit_program_1, JIT_AVX4 | JIT_O1);
    effect_program_create_jit(&jit_program_2, JIT_AVX4 | JIT_O2);
    effect_program_create_jit(&jit_program_3, JIT_AVX4 | JIT_O3);

    //~ test_effects_all(&test_program_3);
    //~ test_effects_all(&test_program_1);

    // compile
    effect_program_compile(&test_program_0, &effects);
    effect_program_compile(&test_program_1, &effects);
    effect_program_compile(&test_program_2, &effects);
    effect_program_compile(&test_program_3, &effects);
    effect_program_compile(&jit_program_0, &effects);
    effect_program_compile(&jit_program_1, &effects);
    effect_program_compile(&jit_program_2, &effects);
    effect_program_compile(&jit_program_3, &effects);

    perf_stop_measurement(&perf_program_creation);

    // execute
	perf_start_measurement(&perf_program_execution);

	test_performance(&test_program_0, &initial_arr);
	test_performance(&test_program_1, &initial_arr);
	test_performance(&test_program_2, &initial_arr);
	test_performance(&jit_program_0, &initial_arr);
	test_performance(&jit_program_1, &initial_arr);
	test_performance(&jit_program_2, &initial_arr);
	test_performance(&jit_program_3, &initial_arr);
    printf("\n");
    //verify(&test_program_1, &ref_program, &initial_arr);
	verify(&jit_program_2, &ref_program, &initial_arr);

    perf_stop_measurement(&perf_program_execution);

    // print mesurements
    //~ perf_print_measurement(&perf_program_creation);
    //~ perf_print_measurement(&perf_program_execution);
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

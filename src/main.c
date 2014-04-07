#include <stdio.h>
#include <stdlib.h>
#include "particle_array.h"
#include "effect_program_naive.h"
#include "particle_vis.h"
#include "performance_measurement.h"

#define LOGGING 1
#define LOGFILE "perf_log"

int main(int argc, char *argv[]) {
    (void) argc; (void) argv;

	FILE* logfile;
	logfile = fopen(LOGFILE, "w");

    perf_measure perf_default = {"default measurement", 0.0, LOGGING, logfile};
    perf_measure perf_program_creation = perf_default;
    perf_program_creation.name = "creation";
    perf_measure perf_program_execution = perf_default;
    perf_program_execution.name = "execution";
    perf_measurement_init();

    particle_vis_init();

    particle_array arr;
    particle_array_create(&arr);
    unsigned int N = 16;
    for(size_t i = 0;i<N;++i) {
        for(size_t j = 0;j<N;++j) {
            particle p;
            p.position[0] = i/(float)N-0.5; p.position[1] = 1; p.position[2] = j/(float)N-0.5;
            p.velocity[0] = 0; p.velocity[1] = 0; p.velocity[2] = 0;
            p.mass = 1;
            p.charge = 1;
            particle_array_add(&arr, p);
        }
    }
    effect_desc effects;
    effect_desc_init(&effects);
    // construct effect array:
    effect_desc_add_linear_accel(&effects,  0, -10, 0);
    effect_desc_add_sphere_bounce(&effects, 0,  0,  0, 0.707, 0.9);
    effect_desc_add_plane_bounce(&effects,  0,  1,  0, -1,    0.9);
    effect_desc_add_plane_bounce(&effects,  0, -1,  0, -1,    0.9);
    effect_desc_add_plane_bounce(&effects,  1,  0,  0, -1,    0.9);
    effect_desc_add_plane_bounce(&effects, -1,  0,  0, -1,    0.9);
    effect_desc_add_plane_bounce(&effects,  0,  0,  1, -1,    0.9);
    effect_desc_add_plane_bounce(&effects,  0,  0, -1, -1,    0.9);
    effect_desc_add_newton_step (&effects);
    // create program:

    perf_start_measurement(&perf_program_creation);

    effect_program naive_program;
    effect_program_create_naive(&naive_program);
    // compile
    effect_program_compile(&naive_program, &effects);

    perf_stop_measurement(&perf_program_creation);

    // execute
	perf_start_measurement(&perf_program_execution);

    for(int i = 0;i<10000;++i) {
        effect_program_execute(&naive_program, &arr, 0.005);
        particle_vis_draw(&arr);
    }

    perf_stop_measurement(&perf_program_execution);

    // clean up:
    effect_desc_destroy(&effects);
    effect_program_destroy(&naive_program);
    particle_array_destroy(&arr);
    particle_vis_deinit();
	fclose(logfile);
    return 0;
}

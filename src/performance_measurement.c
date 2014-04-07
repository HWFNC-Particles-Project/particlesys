#define _POSIX_C_SOURCE 199309L
#include "performance_measurement.h"

#ifdef NO_PERF

int perf_measurement_init(){
	return 0;
}

void perf_start_measurement(perf_measure* pm) {
    (void) pm;
}
void perf_stop_measurement(perf_measure* pm) {
    (void) pm;
}

#else
#include <time.h>
#include <sys/time.h>

uint64_t nanotime() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return 1000000000ull*t.tv_sec + t.tv_nsec;
}

int perf_measurement_init(){
	return 0;
}

void perf_start_measurement(perf_measure* pm){
	// get start measurement
	pm->start = nanotime();
}

void perf_stop_measurement(perf_measure* pm){
	// calcualte time/cycles etc, and print/log it with to the according measurement name
	uint64_t elapsed = nanotime() - pm->start;

	if(pm->bool_log){
		fprintf(pm->logfile, "%s: %f", pm->name, elapsed*1.e-9);
	}
	printf("%s: %f\n", pm->name, elapsed*1.e-9);
}

#endif

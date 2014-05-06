#ifndef PERF_MEASUREMENT
#define PERF_MEASUREMENT
#include <stdio.h>
#include <stdint.h>

#define CPU_FREQ 1.8e9
#define NUM_PART (32*32)
#define COST_LOW 26.5*NUM_PART*NUM_PART+47.0*NUM_PART
#define COST_HIGH 39.0*NUM_PART*NUM_PART+92.0*NUM_PART
#define NUM_ITER 10000

int perf_measurement_init();

typedef struct perf_measure perf_measure;
struct perf_measure{
	char* name;
	uint64_t start,finish;
	int bool_log;
	FILE* logfile;
};

void perf_start_measurement(perf_measure* pm);
void perf_stop_measurement(perf_measure* pm);
void perf_print_measurement(perf_measure* pm);

#endif

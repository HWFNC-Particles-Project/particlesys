#ifndef PERF_MEASUREMENT
#define PERF_MEASUREMENT
#include <stdio.h>

int perf_measurement_init();

typedef struct perf_measure perf_measure;
struct perf_measure{
	char* name;
	double start;
	int bool_log;
	FILE* logfile;
};

void perf_start_measurement(perf_measure* pm);
void perf_stop_measurement(perf_measure* pm);



#endif
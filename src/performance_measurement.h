#ifndef PERF_MEASUREMENT
#define PERF_MEASUREMENT
#include <stdio.h>

int perf_measurement_init();

typedef struct perf_measure perf_measure;
struct perf_measure{
	int bool_log;
	double start;
	char* name;
	FILE* logfile;
};

void start_measurement(perf_measure pm);
void stop_measurement(perf_measure pm);



#endif
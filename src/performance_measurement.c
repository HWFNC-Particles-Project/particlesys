#include "performance_measurement.h"
#include "ftimer.h"

#ifdef NO_PERF

int perf_measurement_init(){
	return 0;
}

void start_measurement(perf_measure pm){}
void stop_measurement(perf_measure pm){}

#else

int perf_measurement_init(){
	init_etime();
//	log = false;
	return 0;
}


void start_measurement(perf_measure pm){
	// get start measurement
	pm.start = get_etime();
}

void stop_measurement(perf_measure pm){
	// calcualte time/cycles etc, and print/log it with to the according measurement name
	double elapsed = get_etime() - pm.start;

	if(pm.bool_log){
		fprintf(pm.logfile, "%s:, %d", pm.name, elapsed);
	}
	printf("%s:, %d", pm.name, elapsed);
}

#endif
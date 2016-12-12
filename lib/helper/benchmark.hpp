#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <time.h>
class Benchmark {
	clock_t start;
public:
	Benchmark() { reset(); }
	void reset() { start = clock(); }
	double seconds_from_start() { return ((double) (clock() - start)) / CLOCKS_PER_SEC; }
};
#endif

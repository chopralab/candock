#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <time.h>
class Benchmark {
	static clock_t start;
public:
	static void reset() { start = clock(); }
	static double seconds_from_start() { return ((double) (clock() - start)) / CLOCKS_PER_SEC; }
};
#endif

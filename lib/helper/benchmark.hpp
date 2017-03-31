#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <ctime>
#include <string>

class __declspec(dllexport) Benchmark {
        clock_t start;
public:
        Benchmark() { reset(); }
        void reset() { start = std::clock(); }
        double seconds_from_start() const { return ((double) (std::clock() - start)) / CLOCKS_PER_SEC; }
        void display_time (const std::string& what) const;
};

#endif

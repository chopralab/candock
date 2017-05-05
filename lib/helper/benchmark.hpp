#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <ctime>
#include <string>

#include <chrono>
#include <iostream>

class Benchmark {
        std::chrono::time_point<std::chrono::system_clock> start;
public:
        Benchmark() { reset(); }
        void reset() { start = std::chrono::system_clock::now(); }
        double seconds_from_start() const {
                std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
                return elapsed_seconds.count();
                
        }

        void display_time(const std::string& what) const {
                std::time_t my_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                std::cout << what << " on " <<  std::ctime(&my_time) << std::endl;
                std::cout << "Clock is at " << seconds_from_start() << std::endl;
        }
};

#endif

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <ctime>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>

class Benchmark {
        clock_t start;
public:
        Benchmark() { reset(); }
        void reset() { start = std::clock(); }
        double seconds_from_start() const { return ((double) (std::clock() - start)) / CLOCKS_PER_SEC; }
		void display_time(const std::string& what) const {
			std::cout << what << " on " << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << std::endl;
			std::cout << "Clock is at " << seconds_from_start() << std::endl;
		}
};

#endif

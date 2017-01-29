#include "benchmark.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>

void Benchmark::display_time(const std::string& what) const
{
        std::cout << what << " on " << boost::posix_time::to_simple_string (boost::posix_time::second_clock::local_time()) << std::endl;
        std::cout << "Clock is at " << seconds_from_start() << std::endl;
}

/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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

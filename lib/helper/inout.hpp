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

#ifndef INOUT_H
#define INOUT_H
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace Inout {
        enum f_not_found {panic=0, no_panic=1};
        int  file_size(const std::string &name);
        void read_file(const std::string &name, std::vector<std::string>&,
                       f_not_found=panic, const int num_occur=-1, 
                       const std::string &pattern=""); // throws Error
        void read_file(const std::string &name, std::string&,
                       f_not_found=panic, const int num_occur=-1, 
                       const std::string &pattern=""); // throws Error
        void read_file(const std::string &name, std::vector<std::string>&,
                       std::streampos &pos_in_file, f_not_found=panic, const int num_occur=-1, 
                       const std::string &pattern=""); // throws Error
        void read_file(const std::string &name, std::string&,
                       std::streampos &pos_in_file, f_not_found=panic, const int num_occur=-1, 
                       const std::string &pattern=""); // throws Error
        void file_open_put_contents(const std::string &name,
                       const std::vector<std::string> &v, std::ios_base::openmode=std::ios_base::out);
        void file_open_put_stream(const std::string &name,
                       const std::stringstream &ss, std::ios_base::openmode=std::ios_base::out);
        std::vector<std::string> files_matching_pattern(const std::string&, const std::string&);

        template<class T>
        void output_file(T &anything, const std::string &filename, std::ios_base::openmode mode=std::ios_base::out) {
                std::stringstream ss;
                ss << anything;
                file_open_put_stream(filename, ss, mode);
        }
        template<class T>
        void output_file(const T &anything, const std::string &filename, std::ios_base::openmode mode=std::ios_base::out) {
                std::stringstream ss;
                ss << anything;
                file_open_put_stream(filename, ss, mode);
	}
};
#endif

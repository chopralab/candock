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
        void read_stream(std::istream &in, std::vector<std::string> &s,
                       std::streampos &pos_in_file, const int num_occur=-1,
                       const std::string &pattern="");
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

#ifndef NOSQLREADER_H
#define NOSQLREADER_H
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <streambuf>
#include <string>

namespace candock {

class NosqlReader {
    std::map<std::string, unsigned int> __hash;
    unsigned int __num;
    unsigned int __hash_num(const std::string& name) {
        auto i = __hash.find(name);
        if (i == __hash.end()) {
            __hash[name] = ++__num;
            return __num;
        } else {
            return i->second;
        }
    }

   public:
    NosqlReader() : __num(0) {}
    void parse_NOSQL(const std::string);
    void parse_dir_of_NOSQL(const std::string);
};
}

#endif

#ifndef NOSQLREADER_H
#define NOSQLREADER_H
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <math.h>
#include <fstream>
#include <streambuf>
#include <string>
#include <map>
#include <iostream>
using namespace std;

namespace candock {

class NosqlReader {
	map<string, unsigned int> __hash;
	unsigned int __num;
	unsigned int __hash_num(const string &name) { auto i = __hash.find(name); if (i == __hash.end()) { __hash[name] = ++__num; return __num; } else { return i->second; } }
public:
	NosqlReader() : __num(0) {}
	void parse_NOSQL(const string);
	void parse_dir_of_NOSQL(const string);
};

}

#endif

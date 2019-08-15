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

class NosqlReader {
	map<string, unsigned int> __hash;
	unsigned int __num;
	unsigned int __hash_num(const string &name) { auto i = __hash.find(name); if (i == __hash.end()) { __hash[name] = ++__num; return __num; } else { return i->second; } }
public:
	NosqlReader() : __num(0) {}
	void parse_NOSQL(const string);
	void parse_dir_of_NOSQL(const string);
};
#endif

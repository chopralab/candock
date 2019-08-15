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

#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <sstream>
#include <algorithm>
#include "nosqlreader.hpp"
#include "helper/error.hpp"
#include "helper/inout.hpp"

void NosqlReader::parse_NOSQL(const string NOSQL_file) {
	vector<string> vec;
	boost::smatch m;
	Inout::read_file(NOSQL_file, vec);
	string bname = boost::filesystem::basename(NOSQL_file);
	const unsigned int num1 = __hash_num(bname);
	for (string &line : vec) {
		if (boost::regex_search(line, m, boost::regex("(\\S+)\\s+\\S+\\s+[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,([^,]+).*"))) {
			if (m[1].matched && m[2].matched) {
				const unsigned int num2 = __hash_num(m[1].str());
				cout << num1 << "\t" << num2 << "\t" << bname << "\t" << m[1].str() << "\t" << m[2].str() << endl;
			}
		}
	}
}

void NosqlReader::parse_dir_of_NOSQL(const string NOSQL_dir) {
	//~ cout << NOSQL_dir << endl;
	vector<string> files = Inout::files_matching_pattern(NOSQL_dir, ".nosql$");
	for (string &f : files) {
		//~ cout << f << endl;
		parse_NOSQL(f);
	}
}

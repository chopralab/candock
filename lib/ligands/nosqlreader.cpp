#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <sstream>
#include <algorithm>
#include "candock/ligands/nosqlreader.hpp"
#include "candock/helper/error.hpp"
#include "candock/helper/inout.hpp"

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

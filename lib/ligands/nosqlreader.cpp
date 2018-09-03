#include "candock/ligands/nosqlreader.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <streambuf>
#include <string>
#include "statchem/helper/error.hpp"
#include "statchem/fileio/inout.hpp"

namespace candock {
void NosqlReader::parse_NOSQL(const std::string NOSQL_file) {
    std::vector<std::string> vec;
    std::smatch m;
    statchem::fileio::read_file(NOSQL_file, vec);
    std::string bname = boost::filesystem::basename(NOSQL_file);
    const unsigned int num1 = __hash_num(bname);
    for (std::string& line : vec) {
        if (std::regex_search(
                line, m, std::regex("(\\S+)\\s+\\S+\\s+[^,]+,[^,]+,[^,]+,[^,]+,"
                                    "[^,]+,[^,]+,[^,]+,[^,]+,([^,]+).*"))) {
            if (m[1].matched && m[2].matched) {
                const unsigned int num2 = __hash_num(m[1].str());
                std::cout << num1 << "\t" << num2 << "\t" << bname << "\t"
                          << m[1].str() << "\t" << m[2].str() << std::endl;
            }
        }
    }
}

void NosqlReader::parse_dir_of_NOSQL(const std::string NOSQL_dir) {
    //~ cout << NOSQL_dir << endl;
    std::vector<std::string> files =
        statchem::fileio::files_matching_pattern(NOSQL_dir, ".nosql$");
    for (std::string& f : files) {
        //~ cout << f << endl;
        parse_NOSQL(f);
    }
}
}

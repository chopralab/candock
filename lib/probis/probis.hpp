#ifndef PROBIS_H
#define PROBIS_H
#include "helper/debug.hpp"
#include <string>
using namespace std;

namespace probis {
	void compare_against_bslib(const string &receptor_file, const string &surface_file,
		const string &receptor_chain_id, const string &bslib_file, const int ncpu, 
		const string &nosql_file, const string &json_file);
}
#endif

#ifndef PROBIS_H
#define PROBIS_H
#include <string>

namespace probis {
void compare_against_bslib(const std::string& receptor_file,
                           const std::string& surface_file,
                           const std::string& receptor_chain_id,
                           const std::string& bslib_file, const int ncpu,
                           const std::string& nosql_file,
                           const std::string& json_file);
}
#endif

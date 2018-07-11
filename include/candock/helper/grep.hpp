#ifndef GREP_HPP
#define GREP_HPP

#include <istream>
#include <regex>

namespace candock {

namespace Grep {
// TODO: Make more like other inout functions???????
size_t count_matches(std::istream& stream, const std::regex& regex);
size_t find_first_case_greater_than(std::istream& stream,
                                    const std::regex& regex, double limit);
std::vector<std::string> search_stream(std::istream& stream,
                                       const std::regex& regex);
}
}

#endif

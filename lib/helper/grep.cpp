#include "grep.hpp"

size_t Grep::count_matches(std::istream& stream, const boost::regex& regex) {
	size_t count = 0;

	boost::smatch what;
	std::string line;

	while(std::getline(stream, line)) {
		bool result = boost::regex_search(line, what, regex);
		if(result) {
			++count;
		}
	}

	return count;
}

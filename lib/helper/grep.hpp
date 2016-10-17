#ifndef GREP_HPP
#define GREP_HPP

#include <istream>
#include <boost/regex.hpp>

namespace Grep {
	// TODO: Make more like other inout functions???????
	size_t count_matches( std::istream& stream, const boost::regex& regex );
}

#endif

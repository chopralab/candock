#include "candock/helper/path.hpp"
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <boost/filesystem.hpp>

namespace candock{
std::string Path::join(const std::string &str1, const std::string &str2) {
	boost::filesystem::path p1 (str1);
	boost::filesystem::path p2 (str2);
	p1 /= p2;
	return  p1.string();
}
}

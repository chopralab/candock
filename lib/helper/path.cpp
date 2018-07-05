#include "candock/helper/path.hpp"
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <boost/filesystem.hpp>


string Path::join(const string &str1, const string &str2) {
	boost::filesystem::path p1 (str1);
	boost::filesystem::path p2 (str2);
	p1 /= p2;
	return  p1.string();
}

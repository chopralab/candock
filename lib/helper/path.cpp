#include "candock/helper/path.hpp"
#include <math.h>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

namespace candock {
std::string Path::join(const std::string& str1, const std::string& str2) {
    boost::filesystem::path p1(str1);
    boost::filesystem::path p2(str2);
    p1 /= p2;
    return p1.string();
}
}

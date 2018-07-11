#ifndef ERROR_H
#define ERROR_H
#include <exception>
#include <sstream>
#include <string>

namespace candock {

class Error : public std::exception {
    const std::string __msg;

   public:
    Error(const std::string& msg) : __msg(msg) {}
    ~Error() throw() {}
    const char* what() const noexcept { return __msg.c_str(); }
};
}

#endif

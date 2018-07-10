#ifndef ERROR_H
#define ERROR_H
#include <string>
#include <sstream>
#include <exception>

namespace candock {

class Error : public std::exception {
	const std::string __msg;
public:
	Error(const std::string &msg) : __msg(msg) {}
	~Error() throw() {}
	const char* what() const noexcept {return __msg.c_str();}
};

}

#endif

#ifndef ERROR_H
#define ERROR_H
#include <string>
#include <sstream>
#include <exception>
using namespace std;

class Error : public exception {
	const string __msg;
public:
	Error(const string &msg) : __msg(msg) {}
	~Error() throw() {}
	const char* what() const noexcept {return __msg.c_str();}
};

#endif

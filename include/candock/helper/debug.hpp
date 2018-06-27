#ifndef CANDOCK_DEBUG_HPP
#define CANDOCK_DEBUG_HPP

#include <iostream>
#include <iomanip>

#ifdef NDEBUG

#define dbgprint(x)
#define dbgmsg(message)

#else

#define dbgprint(x) std::cerr << #x << ": " << x << " in "  << __FILE__ << ":" << __LINE__ << endl
#define dbgmsg(message) std::cerr << message << " in "  << __FILE__ << ":" << __LINE__ <<endl
#endif

#endif

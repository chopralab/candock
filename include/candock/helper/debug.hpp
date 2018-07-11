#ifndef CANDOCK_DEBUG_HPP
#define CANDOCK_DEBUG_HPP

#include <iomanip>
#include <iostream>

#ifdef NDEBUG

#define dbgprint(x)
#define dbgmsg(message)

#else

#define dbgprint(x)                                                       \
    std::cerr << #x << ": " << x << " in " << __FILE__ << ":" << __LINE__ \
              << std::endl
#define dbgmsg(message) \
    std::cerr << message << " in " << __FILE__ << ":" << __LINE__ << std::endl
#endif

#endif

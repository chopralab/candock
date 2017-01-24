#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#ifdef NDEBUG

#define dbgprint(x)
#define dbgmsg(message)

#else

#define dbgprint(x) cerr << #x << ": " << x << " in "  << __FILE__ << ":" << __LINE__ << endl
#define dbgmsg(message) cerr << message << " in "  << __FILE__ << ":" << __LINE__ << endl
#endif

#endif
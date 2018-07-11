#ifndef CANDOCK_EXPORT_H_
#define CANDOCK_EXPORT_H_

#ifdef _WINDOWS
// We don't want to hear about how sprintf is "unsafe".
#pragma warning(disable : 4996)
// Keep MS VC++ quiet about lack of dll export of private members.
#pragma warning(disable : 4251)
#if defined(CANDOCK_SHARED_LIBRARY)
#define CANDOCK_EXPORT __declspec(dllexport)
#elif defined(CANDOCK_STATIC_LIBRARY) || defined(CANDOCK_USE_STATIC_LIBRARIES)
#define CANDOCK_EXPORT
#else
#define CANDOCK_EXPORT \
    __declspec(dllimport)
#endif
#else
#define CANDOCK_EXPORT
#endif

#endif

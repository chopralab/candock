/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef CANDOCK_EXPORT_H_
#define CANDOCK_EXPORT_H_

/*
 * Shared libraries are messy in Visual Studio. We have to distinguish three
 * cases:
 *   (1) this header is being used to build the OpenMM shared library
 *       (dllexport)
 *   (2) this header is being used by a *client* of the OpenMM shared
 *       library (dllimport)
 *   (3) we are building the OpenMM static library, or the client is
 *       being compiled with the expectation of linking with the
 *       OpenMM static library (nothing special needed)
 * In the CMake script for building this library, we define one of the symbols
 *     KB_BUILDING_{SHARED|STATIC}_LIBRARY
 * Client code normally has no special symbol defined, in which case we'll
 * assume it wants to use the shared library. However, if the client defines
 * the symbol OPENMM_USE_STATIC_LIBRARIES we'll suppress the dllimport so
 * that the client code can be linked with static libraries. Note that
 * the client symbol is not library dependent, while the library symbols
 * affect only the OpenMM library, meaning that other libraries can
 * be clients of this one. However, we are assuming all-static or all-shared.
 */

#ifdef _WINDOWS
    // We don't want to hear about how sprintf is "unsafe".
    #pragma warning(disable:4996)
    // Keep MS VC++ quiet about lack of dll export of private members.
    #pragma warning(disable:4251)
    #if defined(CANDOCK_SHARED_LIBRARY)
        #define CANDOCK_EXPORT __declspec(dllexport)
    #elif defined(CANDOCK_STATIC_LIBRARY) || defined(CANDOCK_USE_STATIC_LIBRARIES)
        #define CANDOCK_EXPORT
    #else
        #define CANDOCK_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define CANDOCK_EXPORT // Linux, Mac
#endif

#endif // CANDOCK_EXPORT_H_


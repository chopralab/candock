/* This is candockexport.hpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

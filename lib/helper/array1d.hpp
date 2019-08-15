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

#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <stdio.h>
#include <string.h>

template<typename T>
struct Array1d {
        T *data;
        size_t sz;

        void init (size_t SZ) {
                sz = SZ;
                data = new T[sz];
                memset (data, 0, sz * sizeof (T));
        }

        Array1d (size_t SZ) : data (nullptr) {
                init (SZ);
        }

        Array1d() : data (nullptr), sz (0) {}

        Array1d (const Array1d &other) { // copy
                dbgmsg ("Array1d copy constructor");
                sz = other.sz;
                data = new T[sz];
                memcpy (data, other.data, sz * sizeof (T));
        }

        ~Array1d() {
                if (data) {
                        delete [] data;
                        data = nullptr;
                }
        }

        void reset() {
                memset (data, 0, sz * sizeof (T));
        }
};

template<typename T>
ostream &operator<< (ostream &stream, const Array1d<T> &s) {
        for (size_t i = 0; i < s.sz; ++i) {
                stream << s.data[i] << endl;
        }

        return stream;
}

#endif

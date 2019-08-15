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

#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <helper/debug.hpp>
using namespace std;

template<typename T>
struct Array2d {
        T **data;
        int szi, szj;

        void init (int SZI, int SZJ) {
                szi = SZI;
                szj = SZJ;
                data = new T*[szi];

                for (int i = 0; i < szi; ++i) {
                        data[i] = new T[szj];
                        memset (data[i], 0, szj * sizeof (T));
                }
        }

        Array2d() : data (nullptr), szi (0), szj (0) {
                dbgmsg ("Array2d empty constructor");

        }

        Array2d (int SZ) : data (nullptr) {
                dbgmsg ("Array2d one size constructor");
                init (SZ, SZ);
        }

        Array2d (int SZI, int SZJ) : data (nullptr) {
                dbgmsg ("Array2d two size constructor");
                init (SZI, SZJ);
        }

        Array2d (const Array2d &other) { // copy
                dbgmsg ("Array2d copy constructor");
                szi = other.szi;
                szj = other.szj;
                data = new T*[szi];

                for (int i = 0; i < szi; ++i) {
                        data[i] = new T[szj];
                        memcpy (data[i], other.data[i], szj * sizeof (T));
                }
        }

        ~Array2d() {
                if (data) {
                        for (int i = 0; i < szi; ++i) {
                                if (data[i]) delete [] data[i];

                                data[i] = nullptr;
                        }

                        delete [] data;
                        data = nullptr;
                }
        }
};

template<typename T>
ostream &operator<< (ostream &stream, const Array2d<T> &s) {
        for (int i = 0; i < s.szi; ++i) {
                for (int j = 0; j < s.szj; ++j) {
                        stream << s.data[i][j];
                }

                stream << endl;
        }

        return stream;
}


/**
 * Specialization for Array2d<bool>
 *
 */

// Check windows
#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif


#ifdef ENVIRONMENT64
#define WORD 32
#else
#define WORD 32
#endif

template<>
struct Array2d<bool> {
private:
        unsigned int **data;
        size_t szi, szj;
public:

        size_t get_szi() const {
                return szi;
        }
        size_t get_szj() const {
                return szj;
        }

        bool get (size_t i, size_t j) const {
                return (data[i/WORD][j]&1<< (i%WORD)) != 0;
        }

        void set (size_t i, size_t j) {
                data[i/WORD][j] = data[i/WORD][j]|1<< (i%WORD);
        }

        void unset (size_t i, size_t j) {
                data[i/WORD][j] = data[i/WORD][j]& (~ (1<< (i%WORD)));
        }

        void init (size_t SZI, size_t SZJ) {
                szi = SZI;
                szj = SZJ;
                size_t byszi = 1 + szi / WORD;
                data = new unsigned int *[byszi];

                for (size_t i = 0; i < byszi; ++i) {
                        data[i] = new unsigned int[szj];
                        memset (data[i], 0, szj * sizeof (unsigned int));
                }
        }

        Array2d() : data (nullptr), szi (0), szj (0) {
                dbgmsg ("Array2d<bool> empty constructor");

        }

        Array2d (size_t SZ) : data (nullptr) {
                dbgmsg ("Array2d<bool> one size constructor");
                init (SZ, SZ);
        }

        Array2d (size_t SZI, size_t SZJ) : data (nullptr) {
                dbgmsg ("Array2d<bool> two size constructor");
                init (SZI, SZJ);
        }

        Array2d (const Array2d &other) { // copy
                dbgmsg ("Array2d copy constructor");
                szi = other.szi;
                szj = other.szj;
                size_t byszi = 1 + szi / WORD;
                data = new unsigned int *[byszi];

                for (size_t i = 0; i < byszi; ++i) {
                        data[i] = new unsigned int[szj];
                        memcpy (data[i], other.data[i], szj * sizeof (unsigned int));
                }
        }

        ~Array2d() {
                if (data) {
                        size_t byszi = 1 + szi / WORD;

                        for (size_t i = 0; i < byszi; ++i) {
                                if (data[i]) delete [] data[i];

                                data[i] = nullptr;
                        }

                        delete [] data;
                        data = nullptr;
                }
        }
        friend ostream &operator<< (ostream &stream, const Array2d<bool> &s) {
                for (size_t i = 0; i < s.szi; ++i) {
                        for (size_t j = 0; j < s.szj; ++j) {
                                stream << s.get (i, j);
                        }

                        stream << endl;
                }

                return stream;
        }

};

#endif

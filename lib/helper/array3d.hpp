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

#ifndef ARRAY3D_H
#define ARRAY3D_H

template<typename T>
struct Array3d {
        T ** *data;
        size_t szi, szj, szk;

        void init (size_t SZI, size_t SZJ, size_t SZK) {
                szi = SZI;
                szj = SZJ;
                szk = SZK;
                data = new T **[szi];

                for (size_t i = 0; i < szi; ++i) {
                        data[i] = new T*[szj];

                        for (size_t j = 0; j < szj; ++j) {
                                data[i][j] = new T[szk];
                                memset (data[i][j], 0, szk * sizeof (T));
                        }
                }
        }

        Array3d() : data (nullptr), szi (0), szj (0), szk (0) {}

        Array3d (size_t SZ) : data (nullptr) {
                init (SZ, SZ, SZ);
        }

        Array3d (size_t SZI, size_t SZJ, size_t SZK) : data (nullptr) {
                init (SZI, SZJ, SZK);
        }

        Array3d (const Array3d &other) { // copy
                dbgmsg ("Array3d copy constructor");
                szi = other.szi;
                szj = other.szj;
                szk = other.szk;
                data = new T*[szi];

                for (size_t i = 0; i < szi; ++i) {
                        data[i] = new T*[szj];

                        for (size_t j = 0; j < szj; ++j) {
                                data[i][j] = new T[szk];
                                memcpy (data[i][j], other.data[i], szk * sizeof (T));
                        }
                }
        }

        ~Array3d() {
                if (data) {
                        for (size_t i = 0; i < szi; ++i) {
                                for (size_t j = 0; j < szj; ++j) {
                                        if (data[i][j]) delete [] data[i][j];

                                        data[i][j] = nullptr;
                                }

                                if (data[i]) delete [] data[i];

                                data[i] = nullptr;
                        }

                        delete [] data;
                        data = nullptr;
                }
        }
};

#endif

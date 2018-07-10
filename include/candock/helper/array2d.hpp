#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <stdio.h>
#include <string.h>
#include <assert.h>

namespace candock {

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
        }

        Array2d (int SZ) : data (nullptr) {
                init (SZ, SZ);
        }

        Array2d (int SZI, int SZJ) : data (nullptr) {
                init (SZI, SZJ);
        }

        Array2d (const Array2d &other) { // copy
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
std::ostream &operator<< (std::ostream &stream, const Array2d<T> &s) {
        for (int i = 0; i < s.szi; ++i) {
                for (int j = 0; j < s.szj; ++j) {
                        stream << s.data[i][j];
                }

                stream << std::endl;
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
        }

        Array2d (size_t SZ) : data (nullptr) {
                init (SZ, SZ);
        }

        Array2d (size_t SZI, size_t SZJ) : data (nullptr) {
                init (SZI, SZJ);
        }

        Array2d (const Array2d &other) { // copy
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
        friend std::ostream &operator<< (std::ostream &stream, const Array2d<bool> &s) {
                for (size_t i = 0; i < s.szi; ++i) {
                        for (size_t j = 0; j < s.szj; ++j) {
                                stream << s.get (i, j);
                        }

                        stream << std::endl;
                }

                return stream;
        }

};

}

#endif

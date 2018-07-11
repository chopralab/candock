#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <stdio.h>
#include <string.h>

namespace candock {

template <typename T>
struct Array1d {
    T* data;
    size_t sz;

    void init(size_t SZ) {
        sz = SZ;
        data = new T[sz];
        memset(data, 0, sz * sizeof(T));
    }

    Array1d(size_t SZ) : data(nullptr) { init(SZ); }

    Array1d() : data(nullptr), sz(0) {}

    Array1d(const Array1d& other) {  // copy
        dbgmsg("Array1d copy constructor");
        sz = other.sz;
        data = new T[sz];
        memcpy(data, other.data, sz * sizeof(T));
    }

    ~Array1d() {
        if (data) {
            delete[] data;
            data = nullptr;
        }
    }

    void reset() { memset(data, 0, sz * sizeof(T)); }
};

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Array1d<T>& s) {
    for (size_t i = 0; i < s.sz; ++i) {
        stream << s.data[i] << std::endl;
    }

    return stream;
}
}

#endif

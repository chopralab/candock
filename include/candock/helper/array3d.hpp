#ifndef ARRAY3D_H
#define ARRAY3D_H

namespace candock {

template <typename T>
struct Array3d {
    T*** data;
    size_t szi, szj, szk;

    void init(size_t SZI, size_t SZJ, size_t SZK) {
        szi = SZI;
        szj = SZJ;
        szk = SZK;
        data = new T**[szi];

        for (size_t i = 0; i < szi; ++i) {
            data[i] = new T*[szj];

            for (size_t j = 0; j < szj; ++j) {
                data[i][j] = new T[szk];
                memset(data[i][j], 0, szk * sizeof(T));
            }
        }
    }

    Array3d() : data(nullptr), szi(0), szj(0), szk(0) {}

    Array3d(size_t SZ) : data(nullptr) { init(SZ, SZ, SZ); }

    Array3d(size_t SZI, size_t SZJ, size_t SZK) : data(nullptr) {
        init(SZI, SZJ, SZK);
    }

    Array3d(const Array3d& other) {  // copy
        szi = other.szi;
        szj = other.szj;
        szk = other.szk;
        data = new T*[szi];

        for (size_t i = 0; i < szi; ++i) {
            data[i] = new T*[szj];

            for (size_t j = 0; j < szj; ++j) {
                data[i][j] = new T[szk];
                memcpy(data[i][j], other.data[i], szk * sizeof(T));
            }
        }
    }

    ~Array3d() {
        if (data) {
            for (size_t i = 0; i < szi; ++i) {
                for (size_t j = 0; j < szj; ++j) {
                    if (data[i][j]) delete[] data[i][j];

                    data[i][j] = nullptr;
                }

                if (data[i]) delete[] data[i];

                data[i] = nullptr;
            }

            delete[] data;
            data = nullptr;
        }
    }
};
}

#endif

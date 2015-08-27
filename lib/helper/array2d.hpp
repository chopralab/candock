#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <assert.h>
#include <helper/debug.hpp>


template<typename T>
struct Array2d {
	T **data;
	int szi, szj;

	void init(int SZI, int SZJ) {
		szi = SZI; szj = SZJ;
		data = new T*[szi];
		for (int i = 0; i < szi; ++i) {
			data[i] = new T[szj];
			memset(data[i], 0, szj * sizeof(T));
		}
	}

	Array2d() : data(nullptr), szi(0), szj(0) {
		dbgmsg("Array2d empty constructor");
	
	}

	Array2d(int SZ) : data(nullptr) {
		dbgmsg("Array2d one size constructor");
		init(SZ, SZ);
	}

	Array2d(int SZI, int SZJ) : data(nullptr) {
		dbgmsg("Array2d two size constructor");
		init(SZI, SZJ);
	}

	template<typename U>
	Array2d(const Array2d<U> &other) { // copy
		dbgmsg("Array2d copy constructor");
		szi = other.szi;
		szj = other.szj;
		data = new T*[szi];
		for (int i = 0; i < szi; ++i) {
			data[i] = new T[szj];
			memcpy(data[i], other.data[i], szj * sizeof(T));
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
ostream& operator<< (ostream& stream, const Array2d<T> s) {
	for (int i = 0; i < s.szi; ++i) {
		for (int j = 0; j < s.szj; ++j) {
			stream << s.data[i][j];
		}
		stream << endl;
	}
	return stream;
}

#endif
	

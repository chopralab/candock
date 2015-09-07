#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <stdio.h>
#include <string.h>

template<typename T>
struct Array1d {
	T *data;
	int sz;

	void init(int SZ) {
		sz = SZ;
		data = new T[sz];
		memset(data, 0, sz * sizeof(T));
	}
	
	Array1d(int SZ) : data(nullptr) {
		init(SZ);
	}
	
	Array1d() : data(nullptr), sz(0) {}
	
	Array1d(const Array1d &other) { // copy
		dbgmsg("Array1d copy constructor");
		sz = other.sz;
		data = new T[sz];
		memcpy(data, other.data, sz * sizeof(T));
	}
	
	~Array1d() {
		if (data) {
			delete [] data;
			data = nullptr;
		}
	}
	
	void reset() {
		memset(data, 0, sz * sizeof(T));
	}
};

template<typename T>
ostream& operator<< (ostream& stream, const Array1d<T> &s) {
	for (int i = 0; i < s.sz; ++i) {
		stream << s.data[i] << endl;
	}
	return stream;
}

#endif

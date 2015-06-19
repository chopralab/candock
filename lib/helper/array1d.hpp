#ifndef ARRAY1D_H
#define ARRAY1D_H

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

#endif

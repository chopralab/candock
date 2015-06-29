#ifndef ARRAY2D_H
#define ARRAY2D_H

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
	Array2d(int SZI, int SZJ) : data(nullptr) {
		init(SZI, SZJ);
	}
	Array2d() : data(nullptr), szi(0), szj(0) {}
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

#endif
	

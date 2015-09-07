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

	Array2d(const Array2d &other) { // copy
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
ostream& operator<< (ostream& stream, const Array2d<T> &s) {
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
	int szi, szj;
public:

	const int get_szi() const { return szi; }
	const int get_szj() const { return szj; }
	
	const bool get(int i, int j) const {
	  return data[i/WORD][j]&1<<(i%WORD);
	}
	
	void set(int i, int j) {
	  data[i/WORD][j] = data[i/WORD][j]|1<<(i%WORD);
	} 
	
	void unset(int i, int j) {
	  data[i/WORD][j] = data[i/WORD][j]&(~(1<<(i%WORD)));
	} 
	
	void init(int SZI, int SZJ) {
		szi = SZI; szj = SZJ;
		int byszi = 1 + szi / WORD;
		data = new unsigned int*[byszi];
		for (int i = 0; i < byszi; ++i) {
			data[i] = new unsigned int[szj];
			memset(data[i], 0, szj * sizeof(unsigned int));
		}
	}

	Array2d() : data(nullptr), szi(0), szj(0) {
		dbgmsg("Array2d<bool> empty constructor");
	
	}

	Array2d(int SZ) : data(nullptr) {
		dbgmsg("Array2d<bool> one size constructor");
		init(SZ, SZ);
	}

	Array2d(int SZI, int SZJ) : data(nullptr) {
		dbgmsg("Array2d<bool> two size constructor");
		init(SZI, SZJ);
	}

	Array2d(const Array2d &other) { // copy
		dbgmsg("Array2d copy constructor");
		szi = other.szi;
		szj = other.szj;
		int byszi = 1 + szi / WORD;
		data = new unsigned int*[byszi];
		for (int i = 0; i < byszi; ++i) {
			data[i] = new unsigned int[szj];
			memcpy(data[i], other.data[i], szj * sizeof(unsigned int));
		}
	}
	
	~Array2d() {
		if (data) {
			int byszi = 1 + szi / WORD;
			for (int i = 0; i < byszi; ++i) {
				if (data[i]) delete [] data[i];
				data[i] = nullptr;
			}
			delete [] data;
			data = nullptr;
		}
	}
	friend ostream& operator<< (ostream& stream, const Array2d<bool> &s);
};

#endif
	

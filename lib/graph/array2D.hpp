#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <stdio.h>
#include <string.h>
#include "helper/error.hpp"
#include "helper/help.hpp"
#include "helper/debug.hpp"

template<class T>
class Array2D {
	T **arr;
	size_t dim;
	void destroy() {
		if (arr != nullptr) {
			dbgmsg("destroying address " << arr);
			for (int i = 0; i < dim; i++)
				if (arr[i]) delete [] arr [i];
			delete [] arr;	
		}
		//~ arr = nullptr;
	}
public:
	class Proxy {
		T* _array;
	public:
		Proxy(T* _array) : _array(_array) { }
		T& operator[](int index) {
			return _array[index];
		}
	};
	Array2D() : dim(0), arr(nullptr) {}
	//~ Array2D(Array2D&& o) : arr(std::move(o.arr)), dim(o.dim) {
	//~ Array2D(Array2D&& other) {
		//~ dbgmsg("in array2D move constructor");
		//~ dim = other.dim;
		//~ arr = new T*[dim];
		//~ for (int i = 0; i < dim; i++) {
			//~ arr[i] = new T[dim]; 
			//~ memset(arr[i], 0, dim * sizeof(T));
		//~ }
		//~ for (int i = 0; i < dim; i++) {
			//~ for (int j = 0; j < dim; j++) {
				//~ arr[i][j] = other.arr[i][j];
			//~ }
		//~ }
	//~ }
	//~ Array2D(const Array2D& other) { // copy constructor
		//~ dbgmsg("in array2D copy constructor");
		//~ dim = other.dim;
		//~ arr = new T*[dim];
		//~ for (int i = 0; i < dim; i++) {
			//~ arr[i] = new T[dim]; 
			//~ memset(arr[i], 0, dim * sizeof(T));
		//~ }
		//~ for (int i = 0; i < dim; i++) {
			//~ for (int j = 0; j < dim; j++) {
				//~ arr[i][j] = other.arr[i][j];
			//~ }
		//~ }
	//~ }
	//~ Array2D operator=(const Array2D &rhs) {
		//~ dbgmsg("in array2D copy operator");
	//~ }
	~Array2D() { dbgmsg("destructor of Array2D"); destroy();	}
	Proxy operator[](int index) {
		return Proxy(arr[index]);
	}
	size_t size() const { return dim; }
	void resize(size_t dim) {
		destroy();
		this->dim = dim;
		arr = new T*[dim];
		for (int i = 0; i < dim; i++) {
			arr[i] = new T[dim]; 
			memset(arr[i], 0, dim * sizeof(T));
		}
	}
	T** data() const { return arr; }
};
#endif

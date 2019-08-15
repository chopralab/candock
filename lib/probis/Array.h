/* MIT License
*
* Copyright (c) 2017 Janez Konc
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#ifndef Array_H_INCLUDED
#define Array_H_INCLUDED


#include <algorithm>
#include <iostream>


template<class T, size_t N>
class Array {
	T data[N];
	
public:
	typedef T ElementType;
	
	const Array& operator= (const Array& src) {
		std::copy(src.begin(), src.end(), begin());
		return *this;
	}
	
	inline T& operator[] (size_t i) 					{return data[i];}
	inline const T& operator[] (size_t i) const 		{return data[i];}
	inline size_t size() const 							{return N;}
	inline T* begin() 									{return data;}
	inline const T* begin() const 						{return data;}
	inline T* end() 									{return data+N;}
	inline const T* end() const 						{return data+N;}
};


template<class T, size_t N>
Array<T, N> operator+ (const Array<T, N>& a, const Array<T, N>& b) {
	Array<T, N> temp;
	for (size_t i = 0; i < N; ++i)
		temp[i] = a[i] + b[i];
	return temp;
}


template<class T, size_t N>
Array<T, N> operator- (const Array<T, N>& a, const Array<T, N>& b) {
	Array<T, N> temp;
	for (size_t i = 0; i < N; ++i)
		temp[i] = a[i] - b[i];
	return temp;
}


template<class T, size_t N>
Array<T, N> operator* (const Array<T, N>& a, const Array<T, N>& b) {
	Array<T, N> temp;
	for (size_t i = 0; i < N; ++i)
		temp[i] = a[i] * b[i];
	return temp;
}


template<class T, size_t N>
Array<T, N> operator* (const Array<T, N>& a, T b) {
	Array<T, N> temp;
	for (size_t i = 0; i < N; ++i)
		temp[i] = a[i] * b;
	return temp;
}


template<class C, class T, size_t N>
std::basic_ostream<C>& operator<< (std::basic_ostream<C>& out, const Array<T, N>& ar) {
	out << "[" << ar[0];
	for (size_t i = 1; i < N; ++i)
		out << "," << ar[i];
	return (out << "]");
}

#endif // Array_H_INCLUDED

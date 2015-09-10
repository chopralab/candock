#include "array2d.hpp"

ostream& operator<< (ostream& stream, const Array2d<bool> &s) {
	for (int i = 0; i < s.szi; ++i) {
		for (int j = 0; j < s.szj; ++j) {
			stream << s.get(i, j);
		}
		stream << endl;
	}
	return stream;
}


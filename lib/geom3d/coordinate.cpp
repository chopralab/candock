#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include "matrix.hpp"
#include "coordinate.hpp"
using namespace std;

namespace Geom3D {
	void Coordinate::rotate_inline(const Matrix &matrix) {
		gsl_vector *vec1 = gsl_vector_alloc(3);
		gsl_vector *vec2 = gsl_vector_alloc(3);
		gsl_vector_set(vec1, 0, __x);
		gsl_vector_set(vec1, 1, __y);
		gsl_vector_set(vec1, 2, __z);
		gsl_blas_dgemv(CblasNoTrans, 1, matrix.rota(), vec1, 0, vec2);
		gsl_vector_add(vec2, matrix.trans());
		__x = gsl_vector_get(vec2, 0);
		__y = gsl_vector_get(vec2, 1);
		__z = gsl_vector_get(vec2, 2);
		gsl_vector_free(vec1);
		gsl_vector_free(vec2);
	}
	
	unique_ptr<Coordinate> Coordinate::rotate(const Matrix &matrix) {
		unique_ptr<Coordinate> c(new Coordinate(*this));
		c->rotate_inline(matrix);
		return c;
	}
	
	void Coordinate::inverse_rotate_inline(const Matrix &matrix) {
		gsl_vector *vec1 = gsl_vector_alloc(3);
		gsl_vector *vec2 = gsl_vector_alloc(3);
		gsl_vector_set(vec1, 0, __x);
		gsl_vector_set(vec1, 1, __y);
		gsl_vector_set(vec1, 2, __z);
		gsl_vector_sub(vec1, matrix.trans());
		gsl_blas_dgemv(CblasTrans, 1, matrix.rota(), vec1, 0, vec2);
		__x = gsl_vector_get(vec2, 0);
		__y = gsl_vector_get(vec2, 1);
		__z = gsl_vector_get(vec2, 2);
		gsl_vector_free(vec1);
		gsl_vector_free(vec2);
	}
	
	unique_ptr<Coordinate> Coordinate::inverse_rotate(const Matrix &matrix) {
		unique_ptr<Coordinate> c(new Coordinate(*this));
		c->inverse_rotate_inline(matrix);
		return c;
	}
};

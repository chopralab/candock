#ifndef KABSCH_H
#define KABSCH_H
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include "candock/helper/debug.hpp"
#include "candock/helper/error.hpp"
#include "candock/geom3d/matrix.hpp"
#include "candock/geom3d/coordinate.hpp"

class Kabsch {
	gsl_matrix *__X, *__Y;
	gsl_matrix *__U;
	gsl_vector *__t;
	int __counter, __sz;
	void __gsl_vector_cross(const gsl_vector*, const gsl_vector*, gsl_vector*);
	int __kabsch(unsigned int, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_vector*, double*);
	static const double __NORM_EPS;
public:
	Kabsch (const int sz=0) : __X(nullptr), __Y(nullptr), __U(nullptr), __t(nullptr), __counter(0), __sz(sz) { resize(sz); }
	~Kabsch () { clear(); }
	void resize(const int sz) { __sz = sz; if (sz>0) { clear(); __X = gsl_matrix_alloc(sz, 3); __Y = gsl_matrix_alloc(sz, 3); __U = gsl_matrix_alloc(3, 3); __t = gsl_vector_alloc(3); }}
	void clear() { __counter = 0; if (__X) gsl_matrix_free(__X); if (__Y) gsl_matrix_free(__Y); if (__U) gsl_matrix_free(__U); if (__t) gsl_vector_free(__t); __X = nullptr; __Y = nullptr; __U = nullptr; __t = nullptr; }
	void add_vertex(const Geom3D::Coordinate &c, const Geom3D::Coordinate &d) {
		gsl_matrix_set (__X, __counter, 0, c.x());
		gsl_matrix_set (__X, __counter, 1, c.y());
		gsl_matrix_set (__X, __counter, 2, c.z());
		gsl_matrix_set (__Y, __counter, 0, d.x());
		gsl_matrix_set (__Y, __counter, 1, d.y());
		gsl_matrix_set (__Y, __counter, 2, d.z());
		__counter++;
	}
	void superimpose() {
		if (__kabsch(__sz, __X, __Y, __U, __t, nullptr) != 1)
			throw Error("die : kabsch superimposition failed");
	}
	Geom3D::Matrix get_rota() const { return Geom3D::Matrix(__U, __t); }
};
#endif // KABSCH_H


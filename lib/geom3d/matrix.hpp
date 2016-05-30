#ifndef MATRIX_H
#define MATRIX_H
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <json/json.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include "helper/error.hpp"
using namespace std;

namespace Geom3D {
	class Matrix {
		typedef pair<gsl_matrix*, gsl_vector*> matrix_pair;
		matrix_pair __matrix;
	public:
		typedef tuple<double, double, double, double> matrix_tuple;
		Matrix() : __matrix(make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3))) {}
		Matrix(const double d[9]) : __matrix(make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3))) { for (int i=0; i<9; i++) { gsl_matrix_set(__matrix.first, i/3, i%3, d[i]); } gsl_vector_set_all(__matrix.second, 0); }
		Matrix(gsl_matrix *U, gsl_vector *t) : __matrix(make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3))) { gsl_matrix_memcpy(__matrix.first, U); gsl_vector_memcpy(__matrix.second, t);}
		Matrix(const Matrix &other) { __matrix = make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3)); gsl_matrix_memcpy(__matrix.first, other.rota()); gsl_vector_memcpy(__matrix.second, other.trans()); } // copy constructor
        Matrix& operator=(const Matrix &other) { gsl_matrix_memcpy(__matrix.first, other.rota()); gsl_vector_memcpy(__matrix.second, other.trans()); return *this; } // copy constructor
		Matrix(Json::Value rota, Json::Value trans) { // c++11 :-)
			__matrix = make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3));
		    gsl_matrix_set(__matrix.first, 0, 0, rota[0][0].asDouble());
		    gsl_matrix_set(__matrix.first, 0, 1, rota[0][1].asDouble());
		    gsl_matrix_set(__matrix.first, 0, 2, rota[0][2].asDouble());
		    gsl_matrix_set(__matrix.first, 1, 0, rota[1][0].asDouble());
		    gsl_matrix_set(__matrix.first, 1, 1, rota[1][1].asDouble());
		    gsl_matrix_set(__matrix.first, 1, 2, rota[1][2].asDouble());
		    gsl_matrix_set(__matrix.first, 2, 0, rota[2][0].asDouble());
		    gsl_matrix_set(__matrix.first, 2, 1, rota[2][1].asDouble());
		    gsl_matrix_set(__matrix.first, 2, 2, rota[2][2].asDouble());
			gsl_vector_set(__matrix.second, 0, trans[0].asDouble());
			gsl_vector_set(__matrix.second, 1, trans[1].asDouble());
			gsl_vector_set(__matrix.second, 2, trans[2].asDouble());
		}
		~Matrix() { gsl_matrix_free(__matrix.first); gsl_vector_free(__matrix.second); }
		gsl_matrix* rota() const { return __matrix.first; }
		gsl_vector* trans() const { return __matrix.second; }
		void set_row(const int &row, const matrix_tuple &t) {
			//~ cout << "set_row"<<endl;
			//~ exit(1);
		    gsl_matrix_set(__matrix.first, row, 0, get<0>(t));
		    gsl_matrix_set(__matrix.first, row, 1, get<1>(t));
		    gsl_matrix_set(__matrix.first, row, 2, get<2>(t));
			gsl_vector_set(__matrix.second, row,   get<3>(t));
		}
		friend ostream& operator<< (ostream& stream, const Matrix& m) {
			gsl_matrix *r = m.rota();
			gsl_vector *v = m.trans();
			stream << "[" << gsl_matrix_get(r, 0, 0) << "," << gsl_matrix_get(r, 0, 1) << "," << gsl_matrix_get(r, 0, 2) << "," << gsl_vector_get(v, 0) << "]" << endl;
			stream << "[" << gsl_matrix_get(r, 1, 0) << "," << gsl_matrix_get(r, 1, 1) << "," << gsl_matrix_get(r, 1, 2) << "," << gsl_vector_get(v, 1) << "]" << endl;
			stream << "[" << gsl_matrix_get(r, 2, 0) << "," << gsl_matrix_get(r, 2, 1) << "," << gsl_matrix_get(r, 2, 2) << "," << gsl_vector_get(v, 2) << "]" << endl;
			return stream;
		}
	};
};
#endif

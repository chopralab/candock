#ifndef MATRIX_H
#define MATRIX_H
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector_double.h>
#include <math.h>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include "candock/helper/error.hpp"

namespace Json {
class Value;
}

namespace candock {

namespace geometry {
class Matrix {
    typedef std::pair<gsl_matrix*, gsl_vector*> matrix_pair;
    matrix_pair __matrix;

   public:
    typedef std::tuple<double, double, double, double> matrix_tuple;
    Matrix()
        : __matrix(
              std::make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3))) {}
    Matrix(const double d[9])
        : __matrix(
              std::make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3))) {
        for (int i = 0; i < 9; i++) {
            gsl_matrix_set(__matrix.first, i / 3, i % 3, d[i]);
        }
        gsl_vector_set_all(__matrix.second, 0);
    }
    Matrix(gsl_matrix* U, gsl_vector* t)
        : __matrix(
              std::make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3))) {
        gsl_matrix_memcpy(__matrix.first, U);
        gsl_vector_memcpy(__matrix.second, t);
    }
    Matrix(const Matrix& other) {
        __matrix = std::make_pair(gsl_matrix_alloc(3, 3), gsl_vector_alloc(3));
        gsl_matrix_memcpy(__matrix.first, other.rota());
        gsl_vector_memcpy(__matrix.second, other.trans());
    }  // copy constructor
    Matrix& operator=(const Matrix& other) {
        gsl_matrix_memcpy(__matrix.first, other.rota());
        gsl_vector_memcpy(__matrix.second, other.trans());
        return *this;
    }  // copy constructor
    Matrix(const Json::Value& rota, const Json::Value& trans);
    ~Matrix() {
        gsl_matrix_free(__matrix.first);
        gsl_vector_free(__matrix.second);
    }
    gsl_matrix* rota() const { return __matrix.first; }
    gsl_vector* trans() const { return __matrix.second; }
    void set_row(const int& row, const matrix_tuple& t) {
        //~ cout << "set_row"<<endl;
        //~ exit(1);
        gsl_matrix_set(__matrix.first, row, 0, std::get<0>(t));
        gsl_matrix_set(__matrix.first, row, 1, std::get<1>(t));
        gsl_matrix_set(__matrix.first, row, 2, std::get<2>(t));
        gsl_vector_set(__matrix.second, row, std::get<3>(t));
    }
    friend std::ostream& operator<<(std::ostream& stream, const Matrix& m) {
        gsl_matrix* r = m.rota();
        gsl_vector* v = m.trans();
        stream << "[" << gsl_matrix_get(r, 0, 0) << ","
               << gsl_matrix_get(r, 0, 1) << "," << gsl_matrix_get(r, 0, 2)
               << "," << gsl_vector_get(v, 0) << "]" << std::endl;
        stream << "[" << gsl_matrix_get(r, 1, 0) << ","
               << gsl_matrix_get(r, 1, 1) << "," << gsl_matrix_get(r, 1, 2)
               << "," << gsl_vector_get(v, 1) << "]" << std::endl;
        stream << "[" << gsl_matrix_get(r, 2, 0) << ","
               << gsl_matrix_get(r, 2, 1) << "," << gsl_matrix_get(r, 2, 2)
               << "," << gsl_vector_get(v, 2) << "]" << std::endl;
        return stream;
    }
};
}
}

#endif

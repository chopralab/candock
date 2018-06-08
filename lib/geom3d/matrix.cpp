#include "candock/geom3d/matrix.hpp"
#include <json/json.h>

namespace Geom3D {

    Matrix::Matrix(const Json::Value& rota, const Json::Value& trans) { // c++11 :-)
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

}
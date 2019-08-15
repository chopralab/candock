/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "matrix.hpp"
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
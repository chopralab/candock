/* This is pca.hpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

#ifndef PCA_H
#define PCA_H
#include <tuple>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace candock {

namespace geometry {
typedef std::tuple<gsl_matrix*, gsl_matrix*, gsl_vector*> PrincipalComponents;
PrincipalComponents pca(const gsl_matrix* data, unsigned int L);
}
}

#endif

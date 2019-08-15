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

#include "geom3d/geom3d.hpp"
#include "helper/debug.hpp"
#include "pca.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

using namespace std;

namespace Geom3D {
#ifndef NDEBUG
	string print_matrix(const string &msg, const gsl_matrix *data) {
		stringstream ss;
		ss << msg << endl;
		for (size_t i = 0; i < data->size1; ++i) {
			for (size_t j = 0; j < data->size2; ++j) {
				ss << gsl_matrix_get(data, i, j) << " ";
			}
			ss << endl;
		}
		return ss.str();
	}
	string print_vector(const string &msg, const gsl_vector *vec) {
		stringstream ss;
		ss << msg << endl;
		for (size_t i = 0; i < vec->size; ++i) {
			ss << gsl_vector_get(vec, i) << " ";
			ss << endl;
		}
		return ss.str();
	}
	string print_line(gsl_matrix *eigenmatrix, const gsl_matrix *projection, 
		const gsl_vector *mean) {

		const Geom3D::Vector3 unit_vector = gsl_matrix_column(eigenmatrix, 0).vector;
		dbgmsg(unit_vector);
		vector<double> projected_points;
		for (size_t i = 0; i < projection->size2; ++i) 
			projected_points.push_back(gsl_matrix_get(projection, 0, i)); 
		const Geom3D::Point geom_center = *mean;
		dbgmsg(geom_center);
		auto ret = minmax_element(projected_points.begin(), projected_points.end());
		const double min_value = *ret.first, max_value = *ret.second;
		const double max_dist = max_value - min_value;
		const double interval = max_dist / 11;

		stringstream ss;
		for (int i = 1; i <= 10 ; ++i) {
			const double position = min_value + i * interval;
			ss << "ATOM      1  X   DUM X   1    " 
				<< Geom3D::line_evaluate(geom_center, unit_vector, position).pdb()
				<< "  1.00  0.00           N  " << endl;
		}
		return ss.str();
	}
#endif
	
	PrincipalComponents pca(const gsl_matrix* data, unsigned int L) {
	    /*
	    @param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
	    @param L - dimension reduction
	    */
	    assert(data != nullptr);
	    //~ assert(L > 0 && L < data->size2);
	    assert(L > 0 && L < data->size1);
	    //~ unsigned int i;
	    unsigned int rows = data->size1;
	    unsigned int cols = data->size2;
	    gsl_vector* mean = gsl_vector_alloc(rows);
	 
	    for(unsigned int i = 0; i < rows; i++) {
	        gsl_vector_set(mean, i, gsl_stats_mean(data->data + i * cols, 1, cols));
	    }
		dbgmsg(print_vector("mean", mean));
	    // Get mean-substracted data into matrix mean_substracted_data.
	    gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
	    gsl_matrix_memcpy(mean_substracted_data, data);
	    for(unsigned int i = 0; i < cols; i++) {
	        gsl_vector_view mean_substracted_point_view = gsl_matrix_column(mean_substracted_data, i);
	        gsl_vector_sub(&mean_substracted_point_view.vector, mean);
	    }
	 
	    // Compute Covariance matrix
	    gsl_matrix* covariance_matrix = gsl_matrix_alloc(rows, rows);
	    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
	    //~ gsl_matrix_free(mean_substracted_data);
	 
	    // Get eigenvectors, sort by eigenvalue.
	    gsl_vector* eigenvalues = gsl_vector_alloc(rows);
	    gsl_matrix* eigenvectors = gsl_matrix_alloc(rows, rows);
	    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
	    gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
	    gsl_eigen_symmv_free(workspace);
	    gsl_matrix_free(covariance_matrix);
	 
	    // Sort the eigenvectors
	    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
	    gsl_vector_free(eigenvalues);
	 
	    // Project the original dataset
	    gsl_matrix* projection = gsl_matrix_alloc(L, cols);
	    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
	    gsl_matrix *eigenmatrix = gsl_matrix_alloc(rows, L);
	    //~ gsl_matrix_copy(eigenmatrix, &L_eigenvectors.matrix);
	    gsl_matrix_memcpy(eigenmatrix, &L_eigenvectors.matrix);

		dbgmsg(print_matrix("vectors", eigenmatrix));

	    //~ gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, data, 0.0, projection);
	    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, mean_substracted_data, 0.0, projection);
	 
		dbgmsg(print_matrix("pca matrix", projection));
		dbgmsg(print_line(eigenmatrix, projection, mean));

	    gsl_matrix_free(mean_substracted_data);
	    gsl_matrix_free(eigenvectors);
	
	    // Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L
	    return make_tuple(projection, eigenmatrix, mean);
	}	
};

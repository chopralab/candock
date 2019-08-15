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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include "helper/benchmark.hpp"
#include "helper/error.hpp"
#include "geom3d/geom3d.hpp"
#include <iostream>
#include <exception>
#include <typeinfo>
#include <map>
#include <set>
#include <cmath>
#include "helper/array1d.hpp"
using namespace std;

namespace Molib {

	class Interpolation {
	public:
		static vector<double> derivative(const vector<double> &y, const double step);
		static vector<double> interpolate(const vector<double> &dataX, const vector<double> &dataY, const double step);
		static vector<double> interpolate_bspline(const vector<double> &dataX, const vector<double> &dataY, const double step, const size_t k=4, const size_t ncoeff=24);
	};
};
#endif

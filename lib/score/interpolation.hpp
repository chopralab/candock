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

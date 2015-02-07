#ifndef COORDINATE_H
#define COORDINATE_H
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <tuple>
#include <memory>
#include "helper/error.hpp"
using namespace std;

namespace Geom3D {
	class Matrix;
	class Coordinate {
		double __x, __y, __z;
	public:
		Coordinate() : __x(0), __y(0), __z(0) {}
		Coordinate(double x, double y, double z) : __x(x), __y(y), __z(z) {}
		Coordinate(const double *array) : __x(array[0]), __y(array[1]), __z(array[2]) {}
		Coordinate(const Coordinate& other) : __x(other.__x), __y(other.__y), __z(other.__z) {} // copy constructor
		Coordinate(const gsl_vector& other) : __x(gsl_vector_get(&other, 0)), __y(gsl_vector_get(&other, 1)), __z(gsl_vector_get(&other, 2)) {} // copy constructor
		void set_x(const double &x) { this->__x = x; }
		void set_y(const double &y) { this->__y = y; }
		void set_z(const double &z) { this->__z = z; }
		double x() const { return __x; }
		double y() const { return __y; }
		double z() const { return __z; }
		int i() const { return (int) ::floor(__x); }
		int j() const { return (int) ::floor(__y); }
		int k() const { return (int) ::floor(__z); }
		const Coordinate& crd() const { return *this; }
		Coordinate floor() { return Coordinate(::floor(__x), ::floor(__y), ::floor(__z)); } 
		bool operator==(const Coordinate& right) const { return (fabs(this->x() - right.x()) < numeric_limits<float>::epsilon() && fabs(this->y() - right.y()) < numeric_limits<float>::epsilon() && fabs(this->z() - right.z()) < numeric_limits<float>::epsilon()); }
		bool operator<(const Coordinate& right) const { return (!(*this == right)) && std::tie(this->__x, this->__y, this->__z) < std::tie(right.__x, right.__y, right.__z); }
		static Coordinate cross(const Coordinate& l, const Coordinate& r) { return Coordinate(l.y()*r.z() - l.z()*r.y(), -l.x()*r.z() + l.z()*r.x(), l.x()*r.y() - l.y()*r.x()); }
		static double scalar(const Coordinate& l, const Coordinate& r) { return l.x()*r.x() + l.y()*r.y() + l.z()*r.z(); }
		Coordinate operator+(const double& right) const { return Coordinate(__x + right, __y + right, __z + right); }
		Coordinate operator+(const Coordinate& right) const { return Coordinate(__x + right.x(), __y + right.y(), __z + right.z()); }
		Coordinate operator-(const double& right) const { return Coordinate(__x - right, __y - right, __z - right); }
		Coordinate operator-(const Coordinate& right) const { return Coordinate(__x - right.x(), __y - right.y(), __z - right.z()); }
		Coordinate operator/(const double& right) const { if (right == 0) throw Error("Coordinate::operator/  division by zero\n"); return Coordinate(__x / right, __y / right, __z / right); }
		Coordinate operator*(const double& right) const { return Coordinate(__x * right, __y * right, __z * right); }
		Coordinate operator-() const { return Coordinate(-x(), -y(), -z()); }
		double distance(const Coordinate &c) const { return sqrt(pow(__x-c.x(), 2) + pow(__y-c.y(), 2) + pow(__z-c.z(), 2)); }
		double distance_sq(const Coordinate &c) const { return pow(__x-c.x(), 2) + pow(__y-c.y(), 2) + pow(__z-c.z(), 2); }
		void normalize() { const double length = this->distance(Coordinate(0.0, 0.0, 0.0)); if (length == 0) throw Error("Coordinate::operator/  division by zero\n"); this->__x /= length; this->__y /= length; this->__z /= length; }
		Coordinate norm() { Coordinate c = *this; c.normalize(); return c; }
		string pdb() const { stringstream outs; outs<<fixed<<right<<setprecision(3)<<setw(8)<<__x<<fixed<<right<<setprecision(3)<<setw(8)<<__y<<fixed<<right<<setprecision(3)<<setw(8)<<__z; return outs.str(); }								
		string simple() const { stringstream outs; outs<<fixed<<setprecision(3)<<__x<<" "<<fixed<<setprecision(3)<<__y<<" "<<fixed<<setprecision(3)<<__z; return outs.str(); }								
		string with_underscores() const { stringstream outs; outs<<setprecision(3)<<__x<<"_"<<setprecision(3)<<__y<<"_"<<setprecision(3)<<__z; return outs.str(); }								
		void rotate_inline(const Matrix&);
		unique_ptr<Coordinate> rotate(const Matrix&);
		void inverse_rotate_inline(const Matrix&);
		unique_ptr<Coordinate> inverse_rotate(const Matrix&);
		friend ostream& operator<< (ostream& stream, const Coordinate& c) {
			stream << setprecision(8) << "[" << c.x() << "," << c.y() << "," << c.z() << "]";
			return stream;
		}
	};
};
#endif

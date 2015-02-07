#ifndef GEOM3D_H
#define GEOM3D_H
#include "coordinate.hpp"
#include <vector>
namespace Geom3D {
	typedef Coordinate Point;
	typedef Coordinate Vector3;
	typedef vector<Point> PointVec;
	double degrees(double);
	double radians(double);
	double angle(const Vector3 &, const Vector3 &);
	double angle(const Point &, const Point &, const Point &);
	double dihedral(const Point &, const Point &, const Point &, const Point &);
	Point line_evaluate(const Point &, const Vector3 &, const double);
}
ostream& operator<<(ostream& os, const Geom3D::PointVec &points);

#endif

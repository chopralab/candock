#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <map>
#include <vector>
#include "candock/geometry/coordinate.hpp"

namespace candock {

namespace geometry {
typedef Coordinate Point;
typedef Coordinate Vector3;
typedef std::map<int, geometry::Point::Vec> GridPoints;

double degrees(double);
double radians(double);

double angle(const Vector3&, const Vector3&);
double angle(const Point&, const Point&, const Point&);
double dihedral(const Point&, const Point&, const Point&, const Point&);
Point line_evaluate(const Point&, const Vector3&, const double);

double compute_rmsd_sq(const Point::Vec& crds1, const Point::Vec& crds2);
double compute_rmsd(const Point::Vec& crds1, const Point::Vec& crds2);

Point compute_geometric_center(const geometry::Point::Vec& crds);
Point::Vec uniform_sphere(const int n);

std::ostream& operator<<(std::ostream& os, const geometry::Point::Vec& points);
std::ostream& operator<<(std::ostream& os,
                         const geometry::Point::ConstSet& points);
}
}

#endif

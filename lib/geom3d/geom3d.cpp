#include "geom3d.hpp"

namespace Geom3D {
	double degrees(double radians) { return radians * 57.29577951308232286465; }
	double radians(double degrees) { return degrees / 57.29577951308232286465; }
	double angle(const Vector3 &v0, const Vector3 &v1) { // angle in radians
		double acc = Coordinate::scalar(v0, v1);
		double d0 = v0.distance(Coordinate(0,0,0));
		double d1 = v1.distance(Coordinate(0,0,0));
		if (d0 <= 0 || d1 <= 0)
			return 0;
		acc /= (d0 * d1);
		if (acc > 1)
			acc = 1;
		else if (acc < -1)
			acc = -1;
		return acos(acc);
	}
	double angle(const Point &p0, const Point &p1, const Point &p2) {
		return angle(p0 - p1, p2 - p1);
	}
	double dihedral(const Point &p0, const Point &p1, const Point &p2, const Point &p3) { // in radians
		Vector3 v10 = p1 - p0;
		Vector3 v12 = p1 - p2;
		Vector3 v23 = p2 - p3;
		Vector3 t = Coordinate::cross(v10, v12);
		Vector3 u = Coordinate::cross(v23, v12);
		Vector3 v = Coordinate::cross(u, t);
		double w = Coordinate::scalar(v, v12);
		double acc = angle(u, t);
		if (w < 0)
			acc = -acc;
		return acc;
	}
	Point line_evaluate(const Point &origin, const Vector3 &unit_vector, const double position) {
		return origin + unit_vector * position;
	}
};

ostream& operator<<(ostream& os, const Geom3D::PointVec &points)	{
	for (auto &point : points) {
		os << "ATOM      1   U  DIK     1    " << point.pdb() << endl;
	}
	return os;
}	
ostream& operator<<(ostream& os, const map<int, Geom3D::PointVec> &points)	{
	for (auto &kv : points) {
		const int bsite_id = kv.first;
		os << "MODEL" << setw(9) << right << bsite_id << endl;
		os << kv.second;
		os << "ENDMDL" << endl;
	}
	return os;
}	


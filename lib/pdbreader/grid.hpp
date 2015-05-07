#ifndef GRID_H
#define GRID_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include "geom3d/coordinate.hpp"
#include "helper/debug.hpp"
#include "helper/help.hpp"
using namespace std;

template<class T> // T needs crd(), distance(), element(), and distance(double) interface
class Grid {
public:
	typedef vector<T*> Points;
private:
	Points ***storage;
	int szi, szj, szk;
	Geom3D::Coordinate __min_crd;
	
	Geom3D::Coordinate __correct(const T &point, const double &dist) const {
		Geom3D::Coordinate crd = point.crd() - __min_crd + dist;
		if (dist < 0) {
			if (crd.x() < 0) crd.set_x(0);
			if (crd.y() < 0) crd.set_y(0);
			if (crd.z() < 0) crd.set_z(0);
		} else {
			if (crd.x() >= szi) crd.set_x(szi - 1);
			if (crd.y() >= szj) crd.set_y(szj - 1);
			if (crd.z() >= szk) crd.set_z(szk - 1);
		}
		return crd;
	}
	void __allocate(const int szii, const int szjj, const int szkk) {
		storage = new Points** [szii]; // allocate memory
		for (int i = 0; i < szii; ++i) {
			storage[i] = new Points* [szjj];
			for (int j = 0; j < szjj; ++j) {
				storage[i][j] = new Points [szkk];
			}
		}
	}
	void __deallocate() {
		for (int i = 0; i < szi; ++i) {
			for (int j = 0; j < szj; ++j) {
				for (int k = 0; k < szk; ++k)
					storage[i][j][k].clear();
				delete [] storage[i][j];
			}
			delete [] storage[i];
		}
		delete [] storage;
	}
public:
	Grid& operator=(const Grid &rhs) {
		if (this != &rhs) {
			// Deallocate, allocate new space, copy values...
			this->__deallocate();
			this->__allocate(rhs.szi, rhs.szj, rhs.szk);
			this->__min_crd = rhs.__min_crd;
			this->szi = rhs.szi;
			this->szj = rhs.szj;
			this->szk = rhs.szk;
			for (int i = 0; i < szi; ++i) {
				for (int j = 0; j < szj; ++j) {
					for (int k = 0; k < szk; ++k)
						storage[i][j][k].assign(rhs.storage[i][j][k].begin(), 
							rhs.storage[i][j][k].end());
				}
			}

		}
		return *this;
	}
	Grid() : storage(nullptr), szi(0), szj(0), szk(0) {}
	template<typename P>
	Grid(const P &points) : storage(nullptr), szi(0), szj(0), szk(0) {
		dbgmsg("size of points in grid " << points.size());
		dbgmsg("points is empty " << boolalpha << points.empty());
		if (!points.empty()) {
			// calculate min and max crd
			Geom3D::Coordinate min_crd(HUGE_VAL, HUGE_VAL, HUGE_VAL),
				max_crd(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
			for (auto &point : points) {
				if (point->crd().x() < min_crd.x()) min_crd.set_x(point->crd().x());
				if (point->crd().y() < min_crd.y()) min_crd.set_y(point->crd().y());
				if (point->crd().z() < min_crd.z()) min_crd.set_z(point->crd().z());
				if (point->crd().x() > max_crd.x()) max_crd.set_x(point->crd().x());
				if (point->crd().y() > max_crd.y()) max_crd.set_y(point->crd().y());
				if (point->crd().z() > max_crd.z()) max_crd.set_z(point->crd().z());
			}
			__min_crd = min_crd;
			Geom3D::Coordinate corr = max_crd - min_crd; // calculate size of grid
			szi = corr.i() + 1;
			szj = corr.j() + 1;
			szk = corr.k() + 1;
			this->__allocate(szi, szj, szk);
			for (auto &point : points) { // set data points to grid cells
				Geom3D::Coordinate crd = point->crd() - min_crd;
				storage[crd.i()][crd.j()][crd.k()].push_back(point);
			}
		}
		dbgmsg("points is empty2 " << boolalpha << points.empty());
	}
	~Grid() {
		this->__deallocate();
	}

	Points get_neighbors(const T &point, const double &dist) {
		Geom3D::Coordinate cmin = __correct(point, -dist);
		Geom3D::Coordinate cmax = __correct(point, dist);
		Points points;
		const double dist_sq = pow(dist, 2);
		for (int i = cmin.i(); i <= cmax.i(); ++i)
			for (int j = cmin.j(); j <= cmax.j(); ++j)
				for (int k = cmin.k(); k <= cmax.k(); ++k) {
					for (auto neighbor : storage[i][j][k]) {
						const double d_sq = point.crd().distance_sq(neighbor->crd());
						if (d_sq < dist_sq) {
							if (neighbor != &point) {
								neighbor->distance(d_sq);
								points.push_back(neighbor);
							}
						}
					}
				}
		return points;
	}

	bool has_neighbor_within(const T &point, const double &dist) {
		Geom3D::Coordinate cmin = __correct(point, -dist);
		Geom3D::Coordinate cmax = __correct(point, dist);
		const double dist_sq = pow(dist, 2);
		for (int i = cmin.i(); i <= cmax.i(); ++i)
			for (int j = cmin.j(); j <= cmax.j(); ++j)
				for (int k = cmin.k(); k <= cmax.k(); ++k) {
					for (auto neighbor : storage[i][j][k]) {
						const double d_sq = point.crd().distance_sq(neighbor->crd());
						if ( d_sq < dist_sq) {
							if (neighbor != &point) {
								return true;
							}
						}
					}
				}
		return false;
	}

	bool clashes(const T &point) {
		const double dist = 3.0;
		const double vdw1 = help::vdw_radius[point.idatm_type()];
		Geom3D::Coordinate cmin = __correct(point, -dist);
		Geom3D::Coordinate cmax = __correct(point, dist);
		for (int i = cmin.i(); i <= cmax.i(); ++i)
			for (int j = cmin.j(); j <= cmax.j(); ++j)
				for (int k = cmin.k(); k <= cmax.k(); ++k) {
					for (auto neighbor : storage[i][j][k]) {
						const double d_sq = point.crd().distance_sq(neighbor->crd());
						const double vdw2 = help::vdw_radius[neighbor->idatm_type()];
						dbgmsg("vdw1 = " << vdw1 << " vdw2 = " << vdw2 
								<< " pow = " << pow(0.75 * (vdw1 + vdw2), 2)
								<< " dst_sq = " << d_sq << " dist = " 
								<< sqrt(d_sq));
						if (d_sq < pow(0.75 * (vdw1 + vdw2), 2)) {
							if (neighbor != &point) {
								return true;
							}
						}
					}
				}
		return false;
	}

	Points get_sorted_neighbors(const T &point, const double &dist) {
		Points points = get_neighbors(point, dist);
		sort(points.begin(), points.end(), [](T* i,T* j){ return (i->distance()<j->distance());}); // sort in ascending order
		return points;
	}
};
#endif

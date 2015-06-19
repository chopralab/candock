#ifndef GPOINTS_H
#define GPOINTS_H
#include "geom3d/geom3d.hpp"
#include "helper/array3d.hpp"
#include "helper/error.hpp"
#include "centro/centroids.hpp"

namespace Molib {
	class Score;
};

namespace Docker {

	struct IJK {
		int i, j, k;
		IJK operator+(const IJK& right) const { return IJK{i + right.i, j + right.j, k + right.k}; }
		IJK operator-(const IJK& right) const { return IJK{i - right.i, j - right.j, k - right.k}; }
		friend ostream& operator<<(ostream& os, const IJK &ijk) {
			os << ijk.i << " " << ijk.j << " " << ijk.k << endl;
			return os;
		}
	};
	
	struct Gpoint {
		Geom3D::Point __crd;
		IJK __ijk;
		map<int, double> __energy; // ligand idatm type to energy
		Geom3D::Point& crd() { return __crd; }
		double energy(const int l) { return __energy.at(l); }
		const Geom3D::Point& crd() const { return __crd; }
		IJK& ijk() { return __ijk; }
	};

	typedef vector<Gpoint> GpointVec;

	class Gpoints {
		map<int, GpointVec> __gridpoints;
		//~ Array3d<bool> __gmap;
		Array3d<Gpoint*> __gmap;
		
		void __identify_gridpoints(const double &grid_spacing, const double &radial_check);
		void __identify_gridpoints(const Molib::Score &score, const set<int> &ligand_idatm_types, 
			const Centro::Centroids &centroids, Molib::MolGrid &grid, 
			const double &grid_spacing, const int &dist_cutoff, const double &excluded_radius, 
			const double &max_interatomic_distance);
	
	public:
		Gpoints(const double &grid_spacing, const double &radial_check);
		Gpoints(const Molib::Score &score, const set<int> &ligand_idatm_types, 
			const Centro::Centroids &centroids, Molib::MolGrid &grid, const double &grid_spacing, 
			const int &dist_cutoff, const double &excluded_radius, const double &max_interatomic_distance);
		GpointVec& get_gridpoints0() { if (!__gridpoints.count(0)) throw Error("die : no gridpoints0 found ?"); return __gridpoints[0]; }
		map<int, GpointVec>& get_gridpoints() { return __gridpoints; }
		GpointVec get_gridpoints_as_vec();
		Gpoint& get_center_point();
		//~ Array3d<bool>& get_gmap() { return __gmap; }
		Array3d<Gpoint*>& get_gmap() { return __gmap; }
		
		friend ostream& operator<<(ostream& os, const Gpoints &gpoints);
	};

	ostream& operator<<(ostream& os, const Docker::GpointVec &points);
};


#endif

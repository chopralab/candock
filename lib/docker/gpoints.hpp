#ifndef GPOINTS_H
#define GPOINTS_H
#include "geom3d/geom3d.hpp"
#include "helper/array1d.hpp"
#include "helper/array3d.hpp"
#include "helper/error.hpp"
#include "centro/centroids.hpp"

namespace Molib {
	class Score;
};

namespace Docker {

	class Gpoints {
	public:

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
			Array1d<double> __energy; // ligand idatm type to energy
			Geom3D::Point& crd() { return __crd; }
			double energy(const int l) const { return __energy.data[l]; }
			const Geom3D::Point& crd() const { return __crd; }
			IJK& ijk() { return __ijk; }
		};

		typedef vector<Gpoint> GpointVec;
		typedef vector<Gpoint*> PGpointVec;

	private:

		map<int, GpointVec> __gridpoints;
		Array3d<Gpoint*> __gmap;
		const Molib::Score *__score;
		const set<int> *__ligand_idatm_types;
		
		
		void __identify_gridpoints(const double &grid_spacing, const double &radial_check);
		void __identify_gridpoints(const Centro::Centroids &centroids, Molib::Atom::Grid &grid, 
			const double &grid_spacing, const int &dist_cutoff, const double &excluded_radius, 
			const double &max_interatomic_distance);
	
	public:
		Gpoints(const double &grid_spacing, const double &radial_check);
		Gpoints(const Molib::Score &score, const set<int> &ligand_idatm_types, 
			const Centro::Centroids &centroids, Molib::Atom::Grid &grid, 
			const double &grid_spacing, const int &dist_cutoff, 
			const double &excluded_radius, const double &max_interatomic_distance);
		Gpoints(const Centro::Centroids &centroids, Molib::Atom::Grid &grid, 
			const double &grid_spacing, const int &dist_cutoff, 
			const double &excluded_radius, const double &max_interatomic_distance);
		GpointVec& get_gridpoints0() { try { return __gridpoints.at(0); } catch (const std::out_of_range& oor) { throw Error("die : no gridpoints0 ?"); } }
		const GpointVec& get_gridpoints0() const { try { return __gridpoints.at(0); } catch (const std::out_of_range& oor) { throw Error("die : no gridpoints0 ?"); } }
		//~ const GpointVec& get_gridpoints0() const { return const_cast<const GpointVec&>(get_gridpoints0()); }
		map<int, GpointVec>& get_gridpoints() { return __gridpoints; }
		const Gpoint& get_center_point() const;
		Array3d<Gpoint*>& get_gmap() { return __gmap; }
		
		friend ostream& operator<<(ostream& os, const Gpoints &gpoints);
	};

	ostream& operator<<(ostream& os, const Docker::Gpoints::GpointVec &points);
};


#endif

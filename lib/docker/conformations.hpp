#ifndef CONFORMATIONS_H
#define CONFORMATIONS_H

#include "helper/array2d.hpp"

namespace Docker {

	class Conformations {
	public:
	
		typedef map<int, map<int, map<int, vector<int>>>> ConfMap;
		typedef vector<vector<pair<int, int>>> AtomMatches;
		
		class Conf {
			const AtomMatches &__atom_matches;
			const Molib::AtomVec &__ordered_atoms;
			Gpoints::PGpointVec __points;
		public:
			Conf(const AtomMatches &atom_matches, const Molib::AtomVec &ordered_atoms, Gpoints::PGpointVec &points) : 
				__atom_matches(atom_matches), __ordered_atoms(ordered_atoms), __points(points) {}
			Gpoints::PGpointVec &get_points() { return __points; }
			const Gpoints::PGpointVec &get_points() const { return __points; }
			Gpoints::Gpoint &get_point(const int i) { return *__points[i]; }
			Gpoints::Gpoint &get_point(const int i) const { return *__points.at(i); }
			Molib::Atom &get_atom(const int i) { return *__ordered_atoms[i]; }
			Molib::Atom &get_atom(const int i) const { return *__ordered_atoms.at(i); }
			const AtomMatches &get_atom_matches() const { return __atom_matches; }
			double compute_rmsd(const Conf &other);
			friend ostream& operator<<(ostream& os, const Conf &conf);
		};

	private:

		vector<Conf> __conf_vec;
		Molib::AtomVec __ordered_atoms;
		AtomMatches __atom_matches;
		ConfMap __conf_map;
		Array2d<double> __rmsd;
		
		void __init_conformations(const Molib::Molecule &seed, Gpoints &gpoints, 
			const double &grid_spacing);

	public:
		Conformations(const Molib::Molecule &seed, Gpoints &gpoints, const double &grid_spacing);
		vector<Conf> &get_conformations() { return __conf_vec; }
		ConfMap &get_confmap() { return __conf_map; }
		vector<int> &get_confs_at(const Gpoints::IJK &ijk) { return __conf_map[ijk.i][ijk.j][ijk.k]; }
		double get_rmsd(const int i, const int j) { return __rmsd.data[i][j]; }

		friend ostream& operator<<(ostream& os, const Conformations &conformations);
	};
};

#endif

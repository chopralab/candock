#ifndef CONFORMATIONS_H
#define CONFORMATIONS_H

#include "helper/array2d.hpp"

namespace Docker {

	class Conformations {
	public:
	
		typedef map<int, map<int, map<int, vector<int>>>> ConfMap;
		
	private:

		const Molib::Molecule &__seed;
		
		vector<Gpoints::PGpointVec> __conf_vec;
		ConfMap __conf_map;
		Array2d<double> __rmsd;
		
	public:
		Conformations(const Molib::Molecule &seed, const Gpoints &gpoints, const double &conf_spin, const int num_univec);
			
		vector<Gpoints::PGpointVec> &get_conformations() { return __conf_vec; }
		const vector<Gpoints::PGpointVec> &get_conformations() const { return __conf_vec; }
		ConfMap &get_confmap() { return __conf_map; }
		vector<int> &get_confs_at(const Gpoints::IJK &ijk) { return __conf_map[ijk.i][ijk.j][ijk.k]; }
		double get_rmsd(const int i, const int j) { return __rmsd.data[i][j]; }

		friend ostream& operator<<(ostream& os, const Conformations &conformations);
	};
};

#endif

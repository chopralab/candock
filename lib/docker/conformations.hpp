/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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
		Array2d<double> __rmsd_sq;
		
	public:
		Conformations(const Molib::Molecule &seed, const Gpoints &gpoints, const double &conf_spin, const int num_univec);
			
		vector<Gpoints::PGpointVec> &get_conformations() { return __conf_vec; }
		const vector<Gpoints::PGpointVec> &get_conformations() const { return __conf_vec; }
		ConfMap &get_confmap() { return __conf_map; }
		vector<int> &get_confs_at(const Gpoints::IJK &ijk) { return __conf_map[ijk.i][ijk.j][ijk.k]; }
		double get_rmsd_sq(const int i, const int j) { return __rmsd_sq.data[i][j]; }

		friend ostream& operator<<(ostream& os, const Conformations &conformations);
	};
};

#endif

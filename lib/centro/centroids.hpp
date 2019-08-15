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

#ifndef CENTROIDS_H
#define CENTROIDS_H


#include "geom3d/coordinate.hpp"
#include "ligands/genlig.hpp"

namespace Centro {

	class Centroid {
		Geom3D::Coordinate __centroid;
		double __radial_check;
	public:
		Centroid() {}
		Centroid(const Geom3D::Coordinate centroid, const double radial_check) :
			__centroid(centroid), __radial_check(radial_check) {}
		Geom3D::Coordinate get_centroid() const { return __centroid; }
		double get_radial_check() const { return __radial_check; }
	};
	
	typedef map<int, vector<Centroid>> Centroids;

	Centroids set_centroids(const genlig::BindingSiteClusters &binding_site_clusters, const double centro_clus_rad);
	Centroids set_centroids(const string &centroid_file, const int num_bsites);

	// According to the C++ standard, operator overloads should be done in the same namespace as the
	// object they are overloading
	// See http://clang.llvm.org/compatibility.html#dep_lookup for details
	ostream& operator<<(ostream& os, const Centro::Centroids& centroids);
};

#endif

#ifndef CENTROIDS_H
#define CENTROIDS_H


#include "candock/geometry/coordinate.hpp"
#include "candock/ligands/genlig.hpp"

namespace Centro {

	class Centroid {
		geometry::Coordinate __centroid;
		double __radial_check;
	public:
		Centroid() {}
		Centroid(const geometry::Coordinate centroid, const double radial_check) :
			__centroid(centroid), __radial_check(radial_check) {}
		geometry::Coordinate get_centroid() const { return __centroid; }
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

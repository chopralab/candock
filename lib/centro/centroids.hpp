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

	Centroids set_centroids(const genlig::BindingSiteClusters &binding_site_clusters);
	Centroids set_centroids(const string &centroid_file);

};

ostream& operator<<(ostream& os, const Centro::Centroids& centroids);


#endif

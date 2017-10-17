#include "centroids.hpp"
#include "helper/inout.hpp"
#include "molib/grid.hpp"
#include "molib/molecules.hpp"
#include "helper/benchmark.hpp"
#include "geom3d/geom3d.hpp"
#include "parser/fileparser.hpp"
#include "cluster/greedy.hpp"

namespace Centro {

	/* Centroid stuff
	 * 
	 */
	vector<Centroid> split_binding_site(const Molib::Molecules &binding_site_ligands, const double centro_clus_rad) {

		vector<Centroid> centroids;
		Geom3D::Point::Vec crds = binding_site_ligands.get_crds();
		Geom3D::Point::Vec clustered = Cluster::Cluster::greedy(crds, centro_clus_rad);
		for (auto &point : clustered) {
			// find closest point and calculate distance to this point
			Geom3D::Point closest;
			double min_dist = HUGE_VAL;
			for (auto &point2 : clustered) {
				const double d = point.distance(point2);
				if (d > 0 && d < min_dist) {
					min_dist = d;
					closest = point2;
				}
			}
			centroids.push_back(Centroid(point, point.distance(closest) + 1.0)); // add 1.0 A tolerance
			dbgmsg("ATOM      1  X   GEO X   1    " << point.pdb() << "  1.00  0.00           N  ");
		}
		return centroids;
	}

	Centroids set_centroids(const genlig::BindingSiteClusters &binding_site_clusters, const double centro_clus_rad) {
		Centroids centroids;
		if (binding_site_clusters.empty()) 
			throw Error("die : no binding sites could be predicted for this protein - define its binding site(s) using the centroid option");
		for (auto &kv : binding_site_clusters) {
			const int bsite_id = kv.first;
			const Molib::Molecules &binding_site_ligands = kv.second;
			centroids[bsite_id] = split_binding_site(binding_site_ligands, centro_clus_rad);
			dbgmsg("to better capture the shape of the binding site "
				<< " it has been split into " << centroids[bsite_id].size() << " centroids");
		}
		dbgmsg("setting centroids from ProBiS predicted binding sites : " 
			<< endl << centroids);
		return centroids;
	}
	Centroids set_centroids(const string &centroid_file, const int num_bsites) {
		Centroids centroids;
		vector<string> data;
		Inout::read_file(centroid_file, data);
		for (string &line : data) {
			stringstream ss(line);
			int bsite_id;
			double x, y, z, rc;
			ss >> bsite_id >> x >> y >> z >> rc;
			if (bsite_id <= num_bsites) {
				centroids[bsite_id].push_back(Centroid(Geom3D::Coordinate(x, y, z), rc));
			}
		}
		if (centroids.empty()) 
			throw Error("die: could not find centroid in centroid file " 
				+ centroid_file + "\n");
		dbgmsg("setting centroids from file : " << endl << centroids);
		return centroids;
	}

	ostream& operator<<(ostream& os, const Centro::Centroids& centroids) {
		for (auto &kv : centroids) {
			const int bsite_id = kv.first;
			for (auto &centroid : kv.second) {
				os << bsite_id << " " << centroid.get_centroid().simple() << " " << fixed
					<< setprecision(3) << centroid.get_radial_check() << endl;
			}
		}
		return os;
	}

};

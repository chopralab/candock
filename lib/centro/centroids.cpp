#include "helper/inout.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "geom3d/matrix.hpp"
#include "geom3d/pca.hpp"
#include "geom3d/geom3d.hpp"
#include "kabsch/kabsch.hpp"
#include "score/score.hpp"
#include "pdbreader/pdbreader.hpp"
#include "centroids.hpp"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <iostream>
#include <exception>

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

namespace Centro {

	/* Centroid stuff
	 * 
	 */
	vector<Centroid> split_binding_site(const Molib::Molecules &binding_site_ligands) {

		vector<Centroid> centroids;
		Molib::Atom::Vec atoms = binding_site_ligands.get_atoms();

		gsl_matrix *data = gsl_matrix_alloc(3, atoms.size());
		for (int i = 0; i < atoms.size(); ++i) {
			gsl_matrix_set(data, 0, i, atoms[i]->crd().x());
			gsl_matrix_set(data, 1, i, atoms[i]->crd().y());
			gsl_matrix_set(data, 2, i, atoms[i]->crd().z());
		}

		gsl_matrix *projection, *eigenmatrix;
		gsl_vector *mean;

		// get the first principal component only
		tie(projection, eigenmatrix, mean) = Geom3D::pca(data, 1);

		const Geom3D::Point geom_center = *mean;
		const Geom3D::Vector3 unit_vector = gsl_matrix_column(eigenmatrix, 0).vector;
		vector<double> projected_points;
		const unsigned int cols = projection->size2;

		for (int i = 0; i < cols; ++i) 
			projected_points.push_back(gsl_matrix_get(projection, 0, i)); 

		auto ret = minmax_element(projected_points.begin(), projected_points.end());
		const double min_value = *ret.first, max_value = *ret.second;
		const double max_dist = max_value - min_value;
		const int num_centroids = ceil(max_dist / 5.0); // default radial check is hardcoded to 5.0 A
		dbgmsg("num_centroids = " << num_centroids << " min_value = " 
			<< min_value << " max_dist = " << max_dist);
		const double interval = max_dist / (num_centroids + 1);

		Molib::Atom::Set visited;

		for (int i = 1; i <= num_centroids; ++i) {
			const double position = min_value + i * interval;
			const Geom3D::Point center = Geom3D::line_evaluate(geom_center, 
				unit_vector, position);
			// find the most distant atom left of center (angle is => 90 degrees)
			// exception is the last center where angle is not checked
			double max_l_dist = 0.0;
			for (auto &atom : atoms) { 
				if (i == num_centroids || Geom3D::degrees(Geom3D::angle(atom->crd() - center, unit_vector)) >= 90) {
					if (!visited.count(atom)) {
						visited.insert(atom);
						const double dist = atom->crd().distance(center);
						if (dist > max_l_dist) max_l_dist = dist;
					}
				}
			}
			centroids.push_back(Centroid(center, max_l_dist + 2.0)); // add 2.0 A tolerance
			dbgmsg("original centroid geom_center = " << geom_center
				<< " unit_vector = " << unit_vector << " position = "
				<< position << " : adding splitted centroid[" << i 
				<< "] at center = " << centroids.back().get_centroid() 
				<< " radial_check = " << centroids.back().get_radial_check());
			dbgmsg("ATOM      1  X   GEO X   1    " << center.pdb()	
				<< "  1.00  0.00           N  ");
		}
		gsl_matrix_free(data);
		gsl_matrix_free(projection);
		gsl_matrix_free(eigenmatrix);
		gsl_vector_free(mean);
		return centroids;
	}

	Centroids set_centroids(const genlig::BindingSiteClusters &binding_site_clusters) {
		Centroids centroids;
		if (binding_site_clusters.empty()) 
			throw Error("die : no binding sites could be predicted for this protein - define its binding site(s) using the centroid option");
		for (auto &kv : binding_site_clusters) {
			const int bsite_id = kv.first;
			const Molib::Molecules &binding_site_ligands = kv.second;
			centroids[bsite_id] = split_binding_site(binding_site_ligands);
			dbgmsg("to better capture the shape of the binding site "
				<< " it has been split into " << centroids[bsite_id].size() << " centroids");
		}
		dbgmsg("setting centroids from ProBiS predicted binding sites : " 
			<< endl << centroids);
		return centroids;
	}
	Centroids set_centroids(const string &centroid_file) {
		Centroids centroids;
		vector<string> data;
		inout::Inout::read_file(centroid_file, data);
		for (string &line : data) {
			stringstream ss(line);
			int bsite_id;
			double x, y, z, rc;
			ss >> bsite_id >> x >> y >> z >> rc;
			centroids[bsite_id].push_back(Centroid(Geom3D::Coordinate(x, y, z), rc));
		}
		if (centroids.empty()) 
			throw Error("die: could not find centroid in centroid file " 
				+ centroid_file + "\n");
		dbgmsg("setting centroids from file : " << endl << centroids);
		return centroids;
	}

};

#ifndef GENCLUS_HPP
#define GENCLUS_HPP

#include <string>

namespace genclus {
	void generate_clusters_of_ligands(const std::string &json_file, const std::string &json_with_ligs_file, 
		const std::string &bio_dir, const std::string &names_dir, const bool neighb, 
		const double hetero_clus_rad, const int hetero_min_pts, const double min_z_score, 
		const bool for_gclus=false);
}

#endif

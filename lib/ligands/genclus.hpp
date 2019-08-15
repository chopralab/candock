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

namespace genclus {
	void generate_clusters_of_ligands(const string &json_file, const string &json_with_ligs_file, 
		const string &bio_dir, const string &names_dir, const bool neighb, 
		const double hetero_clus_rad, const int hetero_min_pts, const double min_z_score, 
		const bool for_gclus=false);
}

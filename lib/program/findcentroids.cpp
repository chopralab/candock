#include "findcentroids.hpp"

#include <boost/filesystem.hpp>

#include "probis/probis.hpp"
#include "helper/inout.hpp"
#include "ligands/genclus.hpp"
#include "pdbreader/molecules.hpp"

#include "helper/path.hpp"
#include "options.hpp"

namespace Program {

        FindCentroids::FindCentroids(const Molib::Molecule& receptor) :
                     __receptor( receptor ) {
                if (cmdl.get_string_option("centroid").empty()) {
                        __centroid_file = Path::join(__receptor.name(), "site.cen");
                } else {
                        __centroid_file = cmdl.get_string_option("centroid");
                }
        }

        bool FindCentroids::__can_read_from_files( ) {
                return Inout::file_size( __centroid_file ) > 0;
        }

        void FindCentroids::__read_from_files( ) {
                cout << "Reading " << cmdl.get_int_option("num_bsites") << 
                    " binding sites from " << __centroid_file << endl;
                __result = Centro::set_centroids( __centroid_file, cmdl.get_int_option("num_bsites"));
        }

	void FindCentroids::__continue_from_prev( ) {

		cout << "Running PROBIS for receptor in file: " << __receptor.name() + ".pdb" << endl;

		// Creates an empty nosql file
		Inout::output_file("", Path::join(__receptor.name(), cmdl.get_string_option("nosql"))); // probis local structural alignments

		probis::compare_against_bslib(__receptor.name() + ".pdb",
			Path::join(__receptor.name(), cmdl.get_string_option("srf_file")),
			__receptor.get_chain_ids(Molib::Residue::protein),
			cmdl.get_string_option("bslib"), cmdl.ncpu(),
			Path::join(__receptor.name(), cmdl.get_string_option("nosql")),
			Path::join(__receptor.name(), cmdl.get_string_option("json")));

		genclus::generate_clusters_of_ligands(
			Path::join(__receptor.name(), cmdl.get_string_option("json")),
			Path::join(__receptor.name(), cmdl.get_string_option("jsonwl")),
			cmdl.get_string_option("bio"),
			cmdl.get_string_option("names"),
			cmdl.get_bool_option("neighb"),
			cmdl.get_double_option("probis_clus_rad"),
			cmdl.get_int_option("probis_min_pts"),
			cmdl.get_double_option("probis_min_z_score"));

		auto binding_sites = 
			genlig::generate_binding_site_prediction(
				Path::join(__receptor.name(), cmdl.get_string_option("jsonwl")), 
				cmdl.get_string_option("bio"),
				cmdl.get_int_option("num_bsites"));

		Inout::output_file(binding_sites.first,  Path::join(__receptor.name(), cmdl.get_string_option("lig_clus_file")));
		Inout::output_file(binding_sites.second, Path::join(__receptor.name(), cmdl.get_string_option("z_scores_file")));

		__result = Centro::set_centroids(binding_sites.first, cmdl.get_double_option("centro_clus_rad"));
		Inout::output_file(__result, __centroid_file); // probis local structural alignments
	}
}

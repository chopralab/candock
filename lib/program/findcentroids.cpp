#include "findcentroids.hpp"

#include <boost/filesystem.hpp>

#include "probis/probis.hpp"
#include "helper/inout.hpp"
#include "cmdlnopts.hpp"
#include "ligands/genclus.hpp"
#include "pdbreader/molecules.hpp"

#include "helper/path.hpp"

namespace Program {
	
	bool FindCentroids::__can_read_from_files( const CmdLnOpts& cmdl ) {
		return inout::Inout::file_size( Path::join(__receptor.name(), cmdl.centroid_file()) ) > 0;
	}

	void FindCentroids::__read_from_files( const CmdLnOpts& cmdl ) {
		cout << "Reading " << cmdl.num_bsites() << " binding sites from " << __receptor.name() + "/" + cmdl.centroid_file() << endl;
		__result = Centro::set_centroids( Path::join(__receptor.name(), cmdl.centroid_file()), cmdl.num_bsites());
	}

	void FindCentroids::__continue_from_prev( const CmdLnOpts& cmdl) {

		cout << "Running PROBIS for receptor in file: " << __receptor.name() + ".pdb" << endl;

		// Creates an empty nosql file
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments

		probis::compare_against_bslib(0, 0, __receptor.name() + ".pdb", 
			__receptor.get_chain_ids(Molib::Residue::protein),
			cmdl.bslib_file(),
			cmdl.ncpu(),
			Path::join(__receptor.name(), cmdl.nosql_file()),
			Path::join(__receptor.name(),cmdl.json_file()));

		genclus::generate_clusters_of_ligands(
			Path::join(__receptor.name(), cmdl.json_file()),
			Path::join(__receptor.name(), cmdl.json_with_ligs_file()),
			cmdl.bio_dir(),
			cmdl.names_dir(),
			cmdl.neighb(),
			cmdl.probis_clus_rad(),
			cmdl.probis_min_pts(),
			cmdl.probis_min_z_score());

		auto binding_sites = 
			genlig::generate_binding_site_prediction(
				Path::join(__receptor.name(), cmdl.json_with_ligs_file()), 
				cmdl.bio_dir(),
				cmdl.num_bsites());

		inout::output_file(binding_sites.first,  Path::join(__receptor.name(), cmdl.lig_clus_file()));
		inout::output_file(binding_sites.second, Path::join(__receptor.name(), cmdl.z_scores_file()));

		__result = Centro::set_centroids(binding_sites.first, cmdl.centro_clus_rad());
		inout::output_file(__result, Path::join(__receptor.name(), cmdl.centroid_file())); // probis local structural alignments
	}
}

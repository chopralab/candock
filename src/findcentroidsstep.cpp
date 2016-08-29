#include "findcentroidsstep.hpp"

#include <boost/filesystem.hpp>

#include "probis/probis.hpp"
#include "helper/inout.hpp"
#include "opts_candock.hpp"
#include "ligands/genclus.hpp"
#include "pdbreader/molecules.hpp"

namespace Program {
	
	bool FindCentroidsStep::__can_read_from_files( const CmdLnOpts& cmdl ) {
		return __get_file_size( cmdl.centroid_in_file() );
	}

	void FindCentroidsStep::__read_from_files( const CmdLnOpts& cmdl ) {
		__result = Centro::set_centroids(cmdl.centroid_in_file(), cmdl.num_bsites());
	}

	void FindCentroidsStep::__continue_from_prev( const CmdLnOpts& cmdl, const ProgramStep* prev ) {

		// Creates an empty nosql file
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments

		probis::compare_against_bslib(0, 0, cmdl.receptor_file(), 
			__receptor.get_chain_ids(Molib::Residue::protein), cmdl.bslib_file(), cmdl.ncpu(),
			cmdl.nosql_file(), cmdl.json_file());

		genclus::generate_clusters_of_ligands(cmdl.json_file(), cmdl.json_with_ligs_file(),
			cmdl.bio_dir(), cmdl.names_dir(), cmdl.neighb(), cmdl.probis_clus_rad(),
			cmdl.probis_min_pts(), cmdl.probis_min_z_score());
		auto binding_sites = 
			genlig::generate_binding_site_prediction(cmdl.json_with_ligs_file(), 
			cmdl.bio_dir(), cmdl.num_bsites());

		inout::output_file(binding_sites.first, cmdl.lig_clus_file());
		inout::output_file(binding_sites.second, cmdl.z_scores_file());

		__result = Centro::set_centroids(binding_sites.first, cmdl.centro_clus_rad());	
		inout::output_file(__result, cmdl.centroid_out_file()); // probis local structural alignments
	}
}

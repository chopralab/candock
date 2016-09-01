#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "opts_candock.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "common.hpp"
#include "helper/error.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/pdbreader.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "docker/dock.hpp"
#include "centro/centroids.hpp"
using namespace std;

CmdLnOpts cmdl;

////////////////// TEST BINDING SITE DETECTION ///////////////////////////

int main(int argc, char* argv[]) {
	try {

		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		/* Create empty output files
		 * 
		 */
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments

		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		Centro::Centroids centroids;
		probis::compare_against_bslib(argc, argv, cmdl.receptor_file(), 
			receptors[0].get_chain_ids(Molib::Residue::protein), 
			cmdl.bslib_file(), cmdl.ncpu(),
			cmdl.nosql_file(), cmdl.json_file());
		genclus::generate_clusters_of_ligands(cmdl.json_file(), cmdl.json_with_ligs_file(),
			cmdl.bio_dir(), cmdl.names_dir(), cmdl.neighb(), cmdl.probis_clus_rad(),
			cmdl.probis_min_pts(), cmdl.probis_min_z_score());
		auto binding_sites = 
			genlig::generate_binding_site_prediction(cmdl.json_with_ligs_file(), 
			cmdl.bio_dir(), cmdl.num_bsites());

		//inout::output_file(binding_sites.first, cmdl.lig_clus_file());
		inout::output_file(binding_sites.second, cmdl.z_scores_file());

		centroids = Centro::set_centroids(binding_sites.first, cmdl.centro_clus_rad());	
		inout::output_file(centroids, cmdl.centroid_file()); // probis local structural alignments

		/* Compute atom types for receptor (gaff types not needed since 
		 * they are read from the forcefield xml file)
		 * 
		 */
		receptors.compute_idatm_type();
		
		/* Create receptor grid
		 * 
		 */
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());

		/* Create gridpoints for ALL centroids representing one or more binding sites
		 * 
		 */
		Docker::Gpoints gpoints(centroids, gridrec, 
			cmdl.grid_spacing(), cmdl.dist_cutoff(), cmdl.excluded_radius(), 
			cmdl.max_interatomic_distance());
		inout::output_file(gpoints, cmdl.gridpdb_hcp_file());

		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}

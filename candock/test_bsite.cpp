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
#include "openmm/forcefield.hpp"
#include "openmm/moleculeinfo.hpp"
#include "score/score.hpp"
#include "openmm/omm.hpp"
#include "linker/linker.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"
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
		//~ inout::output_file("", cmdl.gridpdb_hcp_file()); // gridpoints for all binding sites
		//~ inout::output_file("", cmdl.aa_file()); // gridpoints for all binding sites
		//~ inout::output_file("", cmdl.centroid_out_file()); // gridpoints for all binding sites
		//~ inout::output_file("", cmdl.egrid_file()); // output energy grid
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments
		
		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		map<int, vector<common::Centroid>> centroids;
		probis::compare_against_bslib(argc, argv, cmdl.receptor_file(), 
			cmdl.receptor_chain_id(), cmdl.bslib_file(), cmdl.ncpu(),
			cmdl.nosql_file(), cmdl.json_file());
		genclus::generate_clusters_of_ligands(cmdl.json_file(), cmdl.json_with_ligs_file(),
			cmdl.geo_dir(), cmdl.names_dir(), cmdl.neighb(), cmdl.probis_clus_rad(),
			cmdl.probis_min_pts(), cmdl.probis_min_z_score());
		auto binding_sites = 
			genlig::generate_binding_site_prediction(cmdl.json_with_ligs_file(), 
			cmdl.bio_dir(), cmdl.num_bsites());
		inout::output_file(binding_sites.first, cmdl.lig_clus_file());
		inout::output_file(binding_sites.second, cmdl.z_scores_file());
		
		centroids = common::set_centroids(binding_sites.first);	

		inout::output_file(centroids, cmdl.centroid_out_file()); // probis local structural alignments

		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model|Molib::PDBreader::skip_hetatm);
		Molib::Molecules receptors = rpdb.parse_molecule();
		
		receptors[0].filter(Molib::Residue::protein, cmdl.receptor_chain_id());

		/* Compute atom types for receptor (gaff types not needed since 
		 * they are read from the forcefield xml file)
		 * 
		 */
		receptors.compute_idatm_type();
		
		/* Create receptor grid
		 * 
		 */
		Molib::MolGrid gridrec(receptors[0].get_atoms());


		/* Create gridpoints for each binding site represented by a centroid
		 * 
		 */
		map<int, Geom3D::PointVec> gridpoints = common::identify_gridpoints(centroids, 
			gridrec, cmdl.grid_spacing(), cmdl.dist_cutoff(), cmdl.excluded_radius(), 
			cmdl.max_interatomic_distance());
		inout::output_file(gridpoints, cmdl.gridpdb_hcp_file());

		//~ /* Identify amino acids that line each binding site
		 //~ * 
		 //~ */
		//~ map<int, set<Molib::Residue*>> amino_acids = common::identify_amino_acids(centroids, gridrec);
		//~ inout::output_file(amino_acids, cmdl.aa_file(), ios_base::app);

		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}

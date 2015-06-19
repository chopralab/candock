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
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "docker/dock.hpp"
#include "centro/centroids.hpp"
using namespace std;

CmdLnOpts cmdl;

int main(int argc, char* argv[]) {
	try {
		cmdl.init(argc, argv);
		cmdl.display_time("started");
		cout << cmdl << endl;
		/* Create empty output files
		 * 
		 */
		inout::output_file("", cmdl.docked_seeds_file()); // output docked & filtered fragment poses

		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		Centro::Centroids centroids;
		if (cmdl.centroid_in_file().empty()) {
			throw Error("For testing use --centroid option to provide a centroid file");
		} else { // ... or else set binding sites from file
			centroids = Centro::set_centroids(cmdl.centroid_in_file());
		}

		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model|Molib::PDBreader::skip_hetatm);
		Molib::Molecules receptors = rpdb.parse_molecule();
		
		receptors[0].filter(Molib::Residue::protein, cmdl.receptor_chain_id());

		Molib::PDBreader lpdb(cmdl.ligand_file(), Molib::PDBreader::all_models, 
			cmdl.max_num_ligands());

		/* Compute atom types for receptor (gaff types not needed since 
		 * they are read from the forcefield xml file)
		 * 
		 */
		receptors.compute_idatm_type();
		
		/* Create receptor grid
		 * 
		 */
		Molib::MolGrid gridrec(receptors[0].get_atoms());

		/* Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_forcefield_file(cmdl.amber_xml_file());
		receptors[0].prepare_for_mm(ffield, gridrec);

		/* Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */
		Molib::Molecules seeds;
		set<string> added;
		set<int> ligand_idatm_types;

		Molib::Molecules ligands;
		while(lpdb.parse_molecule(ligands)) {
			ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);
			common::create_mols_from_seeds(added, seeds, ligands);
			ligands.clear();
		}
		dbgmsg(seeds);

		/* Erase atom properties since they make the graph matching incorrect
		 *
		 */
		seeds.erase_properties();

		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			gridrec, cmdl.ref_state(),cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond());

		/* Create gridpoints for ALL centroids representing one or more binding sites
		 * 
		 */
		Docker::Gpoints gpoints(score, ligand_idatm_types, centroids, gridrec, 
			cmdl.grid_spacing(), cmdl.dist_cutoff(), cmdl.excluded_radius(), 
			cmdl.max_interatomic_distance());
		inout::output_file(gpoints, cmdl.gridpdb_hcp_file());

		/* Create a zero centered centroid with 10 A radius (max fragment 
		 * radius) for getting all conformations of each seed
		 * 
		 */
		Docker::Gpoints gpoints0(cmdl.grid_spacing(), 10.0);
		//~ inout::output_file(gpoints0, cmdl.gridpdb_hcp_file());

		/* Create template grids using ProBiS-ligands algorithm
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 */

		vector<thread> threads;

		threads.clear();
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(
				thread([&seeds, &gridrec, &score, &gpoints0, &gpoints, i] () {
					// iterate over seeds and dock unique seeds
					for (int j = i; j < seeds.size(); j+= cmdl.ncpu()) {
						try {
							dbgmsg(seeds[j]);
							/* Compute all conformations of this seed with the center
							 * atom fixed on coordinate origin using maximum clique algorithm
							 * 
							 */
							Docker::Conformations conf(seeds[j], gpoints0, cmdl.grid_spacing());

							inout::output_file(conf, "conformations.pdb"); 
							
							/* Dock this seed's conformations to the entire grid
							 * by moving them over all gridpoints and probe where
							 * they clash with the receptor
							 *
							 */
							//~ Docker::Dock dock(gpoints, conf, seeds[j], score);
							Docker::Dock dock(gpoints, conf, seeds[j]);
							Docker::AcceptedConformations accepted = dock.run();

							//~ inout::output_file(accepted, cmdl.docked_seeds_file(), 
								//~ ios_base::app); // output docked & filtered fragment poses
//~ 
							// cluster accepted conformations based on rmsd and 
							// find best-scored cluster representatives
							Molib::Molecules docked_seeds = dock.cluster(accepted);

							inout::output_file(docked_seeds, "tmp/" + seeds[j].name() + "/" 
								+ cmdl.top_seeds_file()); // output docked & clustered seeds
							
							throw Error("exit after dock");
							//~ // cluster non clashing seeds based on rmsd and 
							//~ // find best-scored cluster representatives
							//~ auto clusters_reps_pair = common::cluster_molecules(
								//~ non_clashing_seeds, score, cmdl.clus_rad(), cmdl.min_pts(), 
								//~ cmdl.max_num_clus(), cmdl.max_seeds_to_cluster());
							//~ 
							//~ Molib::Molecules tseeds;
							//~ common::convert_clusters_to_mols(tseeds, clusters_reps_pair.second);
//~ #ifndef NDEBUG
							//~ cluster::MapD<Molib::Molecule> scores = 
								//~ score.many_ligands_score(tseeds);
							//~ dbgmsg(scores);
//~ #endif
							//~ inout::output_file(tseeds, "tmp/" + seeds[j].name() + "/" 
								//~ + cmdl.top_seeds_file()); // output docked & filtered fragment poses
	
						}
						catch (Error& e) {
							cerr << "skipping seed due to : " << e.what() << endl; 
						} 
	
					}
			}));
		}
		for(auto& thread : threads) {
			thread.join();
		}

		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}

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
#include "helper/path.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/nrset.hpp"
#include "pdbreader/pdbreader.hpp"
#include "modeler/forcefield.hpp"
#include "score/score.hpp"
#include "linker/linker.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "docker/dock.hpp"
#include "centro/centroids.hpp"
#include "cluster/greedy.hpp"
#include "modeler/modeler.hpp"

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
		inout::output_file("", cmdl.gridpdb_hcp_file()); // gridpoints for all binding sites
		inout::output_file("", cmdl.prep_file()); // output prepared ligands
		//inout::output_file("", cmdl.docked_file()); // output docked molecule conformations
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
		if (cmdl.centroid_file().empty()) {
			probis::compare_against_bslib(argc, argv, cmdl.receptor_file(), 
				receptors[0].get_chain_ids(Molib::Residue::protein), cmdl.bslib_file(), cmdl.ncpu(),
				cmdl.nosql_file(), cmdl.json_file());
			genclus::generate_clusters_of_ligands(cmdl.json_file(), cmdl.json_with_ligs_file(),
				cmdl.bio_dir(), cmdl.names_dir(), cmdl.neighb(), cmdl.probis_clus_rad(),
				cmdl.probis_min_pts(), cmdl.probis_min_z_score());
			auto binding_sites = 
				genlig::generate_binding_site_prediction(cmdl.json_with_ligs_file(), 
				cmdl.bio_dir(), cmdl.num_bsites());

			inout::output_file(binding_sites.first, cmdl.lig_clus_file());
			inout::output_file(binding_sites.second, cmdl.z_scores_file());
			
			centroids = Centro::set_centroids(binding_sites.first, cmdl.centro_clus_rad());	
			inout::output_file(centroids, cmdl.centroid_file()); // probis local structural alignments

		} else { // ... or else set binding sites from file
			centroids = Centro::set_centroids(cmdl.centroid_file(), cmdl.num_bsites());
		}

		Molib::PDBreader lpdb(cmdl.ligand_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
			cmdl.max_num_ligands());

		/* Compute atom types for receptor and cofactor(s): gaff types for protein, 
		 * Mg ions, and water are read from the forcefield xml file later on while 
		 * gaff types for cofactors (ADP, POB, etc.) are calculated de-novo here
		 * 
		 */
		receptors.compute_idatm_type()
			.compute_hydrogen()
			.compute_bond_order()
			.compute_bond_gaff_type()
			.refine_idatm_type()
			.erase_hydrogen()  // needed because refine changes connectivities
			.compute_hydrogen()   // needed because refine changes connectivities
			.compute_ring_type()
			.compute_gaff_type()
			.compute_rotatable_bonds() // relies on hydrogens being assigned
			.erase_hydrogen();

		/* Create receptor grid
		 * 
		 */
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());

		/* Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */

		vector<thread> threads;
		mutex mtx;
		mutex mtx_read_lock;

		Molib::Molecules seeds;
		set<int> added;
		set<int> ligand_idatm_types;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb, &seeds, &added, &ligand_idatm_types, &mtx, &mtx_read_lock] () {
				Molib::Molecules ligands;
				bool thread_is_not_done = lpdb.parse_molecule(ligands);
				while(thread_is_not_done) {
					// Compute properties, such as idatm atom types, fragments, seeds,
					// rotatable bonds etc.
					ligands.compute_idatm_type()
						.compute_hydrogen()
						.compute_bond_order()
						.compute_bond_gaff_type()
						.refine_idatm_type()
						.erase_hydrogen()  // needed because refine changes connectivities
						.compute_hydrogen()   // needed because refine changes connectivities
						.compute_ring_type()
						.compute_gaff_type()
						.compute_rotatable_bonds() // relies on hydrogens being assigned
						.erase_hydrogen();
					{
						lock_guard<std::mutex> guard(mtx);
						ligands.compute_overlapping_rigid_segments(cmdl.seeds_file());

						ligand_idatm_types = ligands.get_idatm_types(ligand_idatm_types);
						common::create_mols_from_seeds(added, seeds, ligands);
					}
					inout::output_file(ligands, cmdl.prep_file(), ios_base::app);
					ligands.clear();

					lock_guard<std::mutex> guard(mtx_read_lock);
					bool thread_is_done = lpdb.parse_molecule(ligands);
				}
			}));
		}
		for(auto& thread : threads) {
			thread.join();
		}
		dbgmsg(seeds);

		/* Erase atom properties since they make the graph matching incorrect
		 *
		 */
		seeds.erase_properties();

		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		score.define_composition(receptors.get_idatm_types(), ligand_idatm_types)
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		/* Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file())
			.parse_forcefield_file(cmdl.water_xml_file());

		/* Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		receptors[0].prepare_for_mm(ffield, gridrec);

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
		Docker::Gpoints gpoints0(cmdl.grid_spacing(), cmdl.max_frag_radius());

		/* Create template grids using ProBiS-ligands algorithm
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 */

		/* Main loop A with threading support : dock seeds with maximum 
		 * weight clique algorithm to the binding site and filter docked
		 * seeds that clash with the receptor (in the future we could 
		 * give more/less weight to more/less frequent seeds)
		 * 
		 */
		threads.clear();
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(
#ifndef NDEBUG
				thread([&seeds, &gridrec, &score, &gpoints0, &gpoints, i] () {
#else
				thread([&seeds, &gpoints0, &gpoints, i] () {
#endif
					// iterate over docked seeds and dock unique seeds
					for (int j = i; j < seeds.size(); j+= cmdl.ncpu()) {
						try {
							dbgmsg(seeds[j]);
							/* Compute all conformations of this seed with the center
							 * atom fixed on coordinate origin using maximum clique algorithm
							 * 
							 */
							Docker::Conformations conf(seeds[j], gpoints0, 
								cmdl.conf_spin(), // degrees
								cmdl.num_univec() // number of unit vectors
								);

#ifndef NDEBUG
							inout::output_file(conf, "conf_" + seeds[j].name() + ".pdb"); 
#endif			
							/* Dock this seed's conformations to the entire grid by moving them 
							 * over all gridpoints and probe where they clash with the receptor: 
							 * cluster docked conformations based on rmsd and energy and output 
							 * only best-scored cluster representatives
							 *
							 */
#ifndef NDEBUG
							Docker::Dock dock(gpoints, conf, seeds[j], score, gridrec, cmdl.clus_rad());
#else
							Docker::Dock dock(gpoints, conf, seeds[j], cmdl.clus_rad());
#endif

							dock.run();

							//~ inout::output_file(dock.get_docked(), cmdl.top_seeds_dir() + "/" + dock.get_docked().name() + "/" 
								//~ + cmdl.top_seeds_file()); // output docked & clustered seeds
							inout::output_file(dock.get_docked(), Path::join(Path::join(cmdl.top_seeds_dir(), dock.get_docked().name()),
								cmdl.top_seeds_file())); // output docked & clustered seeds
	
						}
						catch (exception& e) {
							cerr << "skipping seed due to : " << e.what() << endl;
						}
					}
			}));
		}
		for(auto& thread : threads) {
			thread.join();
		}

		seeds.clear();

		/** 
		 * Main loop C with threading support : connection of seeds and
		 * minimization
		 * 
		 */
		threads.clear();
		mutex mtx2;
		int ligand_cnt = 0;
		
		Molib::PDBreader lpdb2(cmdl.prep_file(), Molib::PDBreader::all_models, 1);

		/**
		 * Insert topology for cofactors, but not for standard residues
		 * that are already known to forcefield (e.g., amino acid residues)
		 *
		 */
		ffield.insert_topology(receptors[0]); // fixes issue #115

		OMMIface::SystemTopology::loadPlugins();
	
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb2, &receptors, &gridrec, &score, &ffield, &ligand_cnt, &mtx2] () {

				Molib::Molecules ligands;

				while (lpdb2.parse_molecule(ligands)) {
					
					Molib::Molecule &ligand = ligands.first();
					dbgmsg("LINKING LIGAND : " << endl << ligand);
					
					/**
					 * Ligand's resn MUST BE UNIQUE for ffield
					 */
					common::change_residue_name(ligand, mtx2, ligand_cnt); 

					// if docking of one ligand fails, docking of others shall continue...
					try { 
					
						ffield.insert_topology(ligand);
	
						/** 
						 * Read top seeds for this ligand
						 */	 
						Molib::NRset top_seeds = common::read_top_seeds_files(ligand,
							cmdl.top_seeds_dir(), cmdl.top_seeds_file(), cmdl.top_percent());
	
						ligand.erase_properties(); // required for graph matching
						top_seeds.erase_properties(); // required for graph matching
	
						/** 
						 * Jiggle the coordinates by one-thousand'th of an Angstrom to avoid minimization failures
						 * with initial bonded relaxation failed errors
						 */
						top_seeds.jiggle();
	
						/* Init minization options and constants, including ligand and receptor topology
						 *
						 */
						OMMIface::Modeler modeler(ffield, cmdl.fftype(), cmdl.dist_cutoff(),
							cmdl.tolerance(), cmdl.max_iterations(), cmdl.update_freq(), 
							cmdl.position_tolerance(), false, 2.0);
							
						/**
						 * Connect seeds with rotatable linkers, account for symmetry, optimize 
						 * seeds with appendices, minimize partial conformations between linking.
						 * 
						 */
						Linker::Linker linker(modeler, receptors[0], ligand, top_seeds, gridrec, score, 
							cmdl.iterative(), cmdl.dist_cutoff(), cmdl.spin_degrees(), 
							cmdl.tol_seed_dist(), cmdl.lower_tol_seed_dist(), 
							cmdl.upper_tol_seed_dist(), cmdl.max_possible_conf(),
							cmdl.link_iter(), cmdl.clash_coeff(), cmdl.docked_clus_rad(), 
							cmdl.max_allow_energy(), cmdl.max_num_possibles(),
							cmdl.max_clique_size(),	cmdl.max_iterations_final());
							
						Linker::DockedConformation::Vec docks = linker.link();
	
						for (auto &docked : docks) {
							common::change_residue_name(docked.get_ligand(), "CAN"); 
							inout::output_file(Molib::Molecule::print_complex(docked.get_ligand(), docked.get_receptor(), docked.get_energy()), 
								Path::join(cmdl.docked_dir(), ligand.name() + ".pdb" ), ios_base::app); // output docked molecule conformations
						}
					} catch (exception& e) { 
						cerr << "Error: skipping ligand " << ligand.name() << " due to : " << e.what() << endl; 
					} 
					ffield.erase_topology(ligand); // he he
					ligands.clear();
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

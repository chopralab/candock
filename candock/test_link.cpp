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
#include "pdbreader/nrset.hpp"
#include "pdbreader/pdbreader.hpp"
#include "modeler/forcefield.hpp"
#include "score/score.hpp"
#include "linker/linker.hpp"
#include "probis/probis.hpp"
#include "ligands/genclus.hpp"
#include "ligands/genlig.hpp"
#include "cluster/optics.hpp"
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
		inout::output_file("", cmdl.docked_file()); // output docked molecule conformations
		
		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();

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
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());

		/** 
		 * Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */
		set<int> ligand_idatm_types;
		Molib::Molecules ligands;
		while(lpdb.parse_molecule(ligands)) {
			ligand_idatm_types = ligands.get_idatm_types(ligand_idatm_types);
			ligands.clear();
		}

		/**
		 * Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), cmdl.dist_cutoff(), 
			cmdl.step_non_bond());

		score.define_composition(receptors.get_idatm_types(), ligand_idatm_types)
			.process_distributions_file(cmdl.distributions_file())
			.compile_scoring_function()
			.parse_objective_function(cmdl.obj_dir(), cmdl.scale_non_bond());

		/**
		 * Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file())
			.parse_forcefield_file(cmdl.water_xml_file());

		/** 
		 * Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		receptors[0].prepare_for_mm(ffield, gridrec);

		/** 
		 * Main loop C with threading support : connection of seeds and
		 * minimization
		 * 
		 */
		vector<thread> threads;
		threads.clear();
		mutex mtx;
		int ligand_cnt = 0;
		
		Molib::PDBreader lpdb2(cmdl.ligand_file(), Molib::PDBreader::all_models, 1);

		OMMIface::SystemTopology::loadPlugins();
	
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb2, &receptors, &gridrec, &score, &ffield, &ligand_cnt, &mtx] () {

				Molib::Molecules ligands;

				while (lpdb2.parse_molecule(ligands)) {
					
					Molib::Molecule &ligand = ligands.first();
					dbgmsg("LINKING LIGAND : " << endl << ligand);
					
					/**
					 * Ligand's resn MUST BE UNIQUE for ffield
					 */
					common::change_residue_name(ligand, mtx, ligand_cnt); 
					ffield.insert_topology(ligand);
	
					// if docking of one ligand fails, docking of others shall continue...
					try { 
					
						/** 
						 * Read top seeds for this ligand
						 */	 
						Molib::NRset top_seeds = common::read_top_seeds_files(ligand,
							cmdl.top_seeds_dir(), cmdl.top_seeds_file(), cmdl.max_top_seeds());
	
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
							cmdl.top_percent(), cmdl.max_clique_size(),
							cmdl.max_iterations_final());
							
						Linker::DockedConformation::Vec docks = linker.link();
	
						for (auto &docked : docks) {
							common::change_residue_name(docked.get_ligand(), "CAN"); 
							inout::output_file(Molib::Molecule::print_complex(docked.get_ligand(), docked.get_receptor(), docked.get_energy()), 
								cmdl.docked_file(), ios_base::app); // output docked molecule conformations
						}
					} catch (Error& e) { 
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

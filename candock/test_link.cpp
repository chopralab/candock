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
		inout::output_file("", cmdl.mini_file()); // output docked & minimized ligands conformations
		
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
		Molib::Atom::Grid gridrec(receptors[0].get_atoms());

		/* Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */
		set<int> ligand_idatm_types;
		Molib::Molecules ligands;
		while(lpdb.parse_molecule(ligands)) {
			ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);
			ligands.clear();
		}

		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond(),
			cmdl.scale_non_bond());

		/* Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond())
			.parse_forcefield_file(cmdl.amber_xml_file());

		/* Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		receptors[0].prepare_for_mm(ffield, gridrec);

		/* Main loop C with threading support : connection of seeds and
		 * minimization
		 * 
		 */
		vector<thread> threads;
		threads.clear();

		Molib::PDBreader lpdb2(cmdl.ligand_file(), Molib::PDBreader::all_models, 1);

		OMMIface::SystemTopology::loadPlugins();
		
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb2, &receptors, &gridrec, &score, &ffield] () {

				OMMIface::ForceField ffield_copy = ffield; // make a local copy of the forcefield
				Molib::Molecules ligands;

				while (lpdb2.parse_molecule(ligands)) {
					Molib::Molecule &ligand = ligands.first();

					ffield_copy.insert_topology(ligand);
					dbgmsg("LINKING LIGAND : " << endl << ligand);
					try { // if ligand fails docking of others continues...

						// read top seeds for this ligand
						Molib::NRset top_seeds = common::read_top_seeds_files(ligand,
							cmdl.top_seeds_file());
						
						ligand.erase_properties(); // required for graph matching
						top_seeds.erase_properties(); // required for graph matching
						
						/**
						 * Init minization options and constants, including ligand 
						 * and receptor topology
						 */
						 						 
						OMMIface::Modeler modeler(ffield_copy, cmdl.fftype(), cmdl.dist_cutoff(),
							cmdl.tolerance(), cmdl.max_iterations(), cmdl.update_freq(), 
							cmdl.position_tolerance(), false, 2.0);
						
						modeler.add_topology(receptors[0].get_atoms());
						modeler.add_topology(ligand.get_atoms());
						
						modeler.init_openmm();

						modeler.add_crds(ligand.get_atoms(), ligand.get_crds());

						/**
						 * Connect seeds with rotatable linkers, account for symmetry, optimize 
						 * seeds with appendices, minimize partial conformations between linking. 
						 * A graph of segments is constructed in which each rigid segment is 
						 * a vertex & segments that are part of seeds are identified.
						 * 
						 */					
						
						Molib::Internal ic(ligand.get_atoms());
						
						Linker::Linker linker(modeler, receptors[0], ligand, top_seeds, gridrec, score, ic, 
							cmdl.dist_cutoff(), cmdl.spin_degrees(), cmdl.tol_seed_dist(), 
							cmdl.max_possible_conf(), cmdl.link_iter(), cmdl.clash_coeff(),
							cmdl.docked_clus_rad(), cmdl.max_allow_energy(), cmdl.iterative(),
							cmdl.max_num_possibles());

						Linker::Partial::Vec partial_conformations = linker.init_conformations();

						/**
						 * Final minimization of each docked ligand conformation
						 * with full ligand and receptor flexibility
						 * 
						 */
						for (auto &partial : partial_conformations) {
						
							try {
								
								Linker::DockedConformation docked = linker.compute_conformation(partial);

								docked.get_receptor().undo_mm_specific();
								
								modeler.add_crds(receptors[0].get_atoms(), docked.get_receptor().get_crds());
								modeler.add_crds(ligand.get_atoms(), docked.get_ligand().get_crds());
								
								modeler.init_openmm_positions();
								
								modeler.unmask(receptors[0].get_atoms());
								modeler.unmask(ligand.get_atoms());
				
								modeler.set_max_iterations(cmdl.max_iterations_final()); // until converged
								modeler.minimize_state(ligand, receptors[0], score);
	
								// init with minimized coordinates
								Molib::Molecule minimized_receptor(receptors[0], modeler.get_state(receptors[0].get_atoms()));
								Molib::Molecule minimized_ligand(ligand, modeler.get_state(ligand.get_atoms()));
				
								minimized_receptor.undo_mm_specific();
								
								Molib::Atom::Grid gridrec(minimized_receptor.get_atoms());

								const double energy = score.non_bonded_energy(gridrec, minimized_ligand);

								inout::output_file(Molib::Molecule::print_complex(docked.get_ligand(), docked.get_receptor(), docked.get_energy()), 
									cmdl.docked_file(), ios_base::app); // output docked molecule conformations
	
								inout::output_file(Molib::Molecule::print_complex(minimized_ligand, minimized_receptor, energy), 
									cmdl.mini_file(), ios_base::app); // output minimized docked conformations

							} catch(Linker::Linker::ConnectionError &e) {
								cerr << "skipping ligand " << ligand.name() << " due to : " << e.what() << endl;
							}
							catch(OMMIface::Modeler::MinimizationError &e) {
								cerr << "skipping ligand " << ligand.name() << " due to : " << e.what() << endl;								
							}
						}
					}
					catch (Error& e) { 
						cerr << "skipping ligand " << ligand.name() << " due to : " << e.what() << endl; 
					} 
					catch (out_of_range& e) { 
						cerr << "skipping ligand due to : " << e.what() << endl; 
						cerr << "This is the ligand that failed: " << endl << ligands << endl;
					}	
					
					ffield_copy.erase_topology(ligand); // he he
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

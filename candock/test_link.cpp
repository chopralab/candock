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
		inout::output_file("", cmdl.egrid_file()); // output energy grid
		inout::output_file("", cmdl.docked_seeds_file()); // output docked & filtered fragment poses
		inout::output_file("", cmdl.docked_ligands_file()); // output docked molecule conformations
		inout::output_file("", cmdl.mini_ligands_file()); // output docked & minimized ligands conformations
		inout::output_file(OMMIface::print_energies_title(), 
			cmdl.energy_file()); // output energy components of docked & minimized ligands conformations
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments
		
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
			gridrec, cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond());

		/* Forcefield stuff : create forcefield for small molecules (and KB 
		 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
		 * forcefield object
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
			.add_kb_forcefield(score, cmdl.step_non_bond(), cmdl.scale_non_bond())
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

		OMMIface::OMM::loadPlugins();
		
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
						

						/* Connect seeds with rotatable linkers, symmetry, optimize 
						 * seeds with appendices. A graph of segments is constructed in 
						 * which each rigid segment is a vertex & segments that are 
						 * part of seeds are identified.
						 * 
						 */
						
						Molib::Internal ic(ligand.get_atoms());

						// compute max linker lengths for pre-computed linkers
						// generate conneciton graph with fragments
						// use maximum clique algorithm to find all cliques in connection graph
						// score found solutions
						// link them using suitable pre-computed linkers
						
						Molib::Linker linker(ligand, top_seeds, gridrec, score, ic, 
							cmdl.dist_cutoff(), cmdl.spin_degrees(), cmdl.tol_dist(),
							cmdl.tol_max_coeff(), cmdl.tol_min_coeff(), 
							cmdl.max_possible_conf(), cmdl.link_iter());

						Molib::Molecules docked = linker.connect();

						if (!docked.empty()) {

							inout::output_file(docked, cmdl.docked_ligands_file(), 
								ios_base::app); // output docked molecule conformations

							//~ /* Cluster docked conformations and take only cluster
							 //~ * representatives for further minimization
							 //~ * 
							 //~ */
							//~ 
							//~ auto clusters_reps_pair = common::cluster_molecules(docked, score,
								//~ cmdl.docked_clus_rad(), cmdl.docked_min_pts(), cmdl.docked_max_num_clus());
							//~ Molib::Molecules docked_representatives; 
							//~ common::convert_clusters_to_mols(docked_representatives, clusters_reps_pair.second);
//~ 
							//~ /* Minimize each representative docked ligand conformation
							//~ * with full flexibility of both ligand and receptor
							 //~ * 
							 //~ */
							//~ Molib::Molecules mini;
							//~ OMMIface::Energies energies;
	//~ 
							//~ for (auto &docked_ligand : docked_representatives) {
								//~ OMMIface::OMM omm(receptors[0], docked_ligand, ffield_copy, 
									//~ cmdl.fftype(), cmdl.dist_cutoff());
								//~ omm.minimize(cmdl.tolerance(), cmdl.max_iterations()); // minimize
								//~ auto ret = omm.get_state(receptors[0], docked_ligand);
								//~ Molib::Molecule &minimized_receptor = mini.add(new Molib::Molecule(ret.first));
								//~ Molib::Molecule &minimized_ligand = mini.add(new Molib::Molecule(ret.second));
								//~ energies[&minimized_ligand] = omm.get_energy_components(
									//~ minimized_receptor, minimized_ligand, cmdl.dist_cutoff());
								//~ minimized_receptor.undo_mm_specific();
							//~ }
							//~ inout::output_file(mini, cmdl.mini_ligands_file(), 
								//~ ios_base::app); // output docked & minimized ligands
							//~ inout::output_file(energies, cmdl.energy_file(), 
								//~ ios_base::app); // output energies of docked & minimized ligands
						}
					}
					catch (Error& e) { 
						cerr << "skipping ligand due to : " << e.what() << endl; 
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

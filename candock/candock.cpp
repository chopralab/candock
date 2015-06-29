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
		inout::output_file("", cmdl.gridpdb_hcp_file()); // gridpoints for all binding sites
		inout::output_file("", cmdl.prep_file()); // output prepared ligands
		inout::output_file("", cmdl.egrid_file()); // output energy grid
		inout::output_file("", cmdl.docked_seeds_file()); // output docked & filtered fragment poses
		inout::output_file("", cmdl.docked_ligands_file()); // output docked molecule conformations
		inout::output_file("", cmdl.mini_ligands_file()); // output docked & minimized ligands conformations
		inout::output_file(OMMIface::print_energies_title(), 
			cmdl.energy_file()); // output energy components of docked & minimized ligands conformations
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments
		
		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		Centro::Centroids centroids;
		if (cmdl.centroid_in_file().empty()) {
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

			centroids = Centro:set_centroids(binding_sites.first);	
			inout::output_file(centroids, cmdl.centroid_out_file()); // probis local structural alignments

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

		Molib::PDBreader lpdb(cmdl.ligand_file(), 
			Molib::PDBreader::all_models|Molib::PDBreader::hydrogens, 
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

		/* Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */

		vector<thread> threads;
		mutex mtx;

		Molib::Molecules seeds;
		set<string> added;
		set<int> ligand_idatm_types;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(thread([&lpdb, &seeds, &added, &ligand_idatm_types, &mtx] () {
				Molib::Molecules ligands;
				while(lpdb.parse_molecule(ligands)) {
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
						.erase_hydrogen()
						.compute_overlapping_rigid_segments();
					{
						lock_guard<std::mutex> guard(mtx);
						ligands.compute_seeds(cmdl.seeds_file());
						ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);
						common::create_mols_from_seeds(added, seeds, ligands);
					}
					inout::output_file(ligands, cmdl.prep_file(), ios_base::app);
					ligands.clear();
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
		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			gridrec, cmdl.ref_state(), cmdl.comp(), cmdl.rad_or_raw(), 
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
				thread([&seeds, &gridrec, &score, &gpoints0, &gpoints, i] () {
					// iterate over docked seeds and dock unique seeds
					for (int j = i; j < seeds.size(); j+= cmdl.ncpu()) {
						try {
							dbgmsg(seeds[j]);
							/* Compute all conformations of this seed with the center
							 * atom fixed on coordinate origin using maximum clique algorithm
							 * 
							 */
							Docker::Conformations conf(seeds[j], gpoints0, cmdl.grid_spacing());

							inout::output_file(conf, "conformations.pdb"); 
							
							/* Dock this seed's conformations to the entire grid by moving them 
							 * over all gridpoints and probe where they clash with the receptor: 
							 * cluster docked conformations based on rmsd and energy and output 
							 * only best-scored cluster representatives
							 *
							 */
							Docker::Dock dock(gpoints, conf, seeds[j], 2.0);
							Molib::Molecules docked_seeds = dock.run();

							inout::output_file(docked_seeds, "tmp/" + seeds[j].name() + "/" 
								+ cmdl.top_seeds_file()); // output docked & clustered seeds
	
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

		seeds.clear();

		/* Main loop C with threading support : connection of seeds and
		 * minimization
		 * 
		 */
		Molib::PDBreader lpdb2(cmdl.prep_file(), Molib::PDBreader::all_models, 1);
		OMMIface::OMM::loadPlugins();

		threads.clear();
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
						Molib::Linker linker(ligand, top_seeds, gridrec, score, ic, 
							cmdl.dist_cutoff(), cmdl.spin_degrees(), cmdl.tol_dist(),
							cmdl.tol_max_coeff(), cmdl.tol_min_coeff(), 
							cmdl.max_possible_conf(), cmdl.link_iter());

						Molib::Molecules docked = linker.connect();

						if (!docked.empty()) {

							inout::output_file(docked, cmdl.docked_ligands_file(), 
								ios_base::app); // output docked molecule conformations
	
							/* Cluster docked conformations and take only cluster
							 * representatives for further minimization
							 * 
							 */
							auto clusters_reps_pair = common::cluster_molecules(docked, score,
								cmdl.docked_clus_rad(), cmdl.docked_min_pts(), cmdl.docked_max_num_clus());
							Molib::Molecules docked_representatives; 
							common::convert_clusters_to_mols(docked_representatives, clusters_reps_pair.second);
							
							/* Minimize each representative docked ligand conformation
							 * with full flexibility of both ligand and receptor
							 * 
							 */
							Molib::Molecules mini;
							OMMIface::Energies energies;
	
							for (auto &docked_ligand : docked_representatives) {
								OMMIface::OMM omm(receptors[0], docked_ligand, ffield_copy, 
									cmdl.fftype(), cmdl.dist_cutoff());
								omm.minimize(cmdl.tolerance(), cmdl.max_iterations()); // minimize
								auto ret = omm.get_state(receptors[0], docked_ligand);
								Molib::Molecule &minimized_receptor = mini.add(new Molib::Molecule(ret.first));
								Molib::Molecule &minimized_ligand = mini.add(new Molib::Molecule(ret.second));
								energies[&minimized_ligand] = omm.get_energy_components(
									minimized_receptor, minimized_ligand, cmdl.dist_cutoff());
								minimized_receptor.undo_mm_specific();
							}
							inout::output_file(mini, cmdl.mini_ligands_file(), 
								ios_base::app); // output docked & minimized ligands
							inout::output_file(energies, cmdl.energy_file(), 
								ios_base::app); // output energies of docked & minimized ligands
						}
					}
					catch (Error& e) { cerr << "skipping ligand due to : " << e.what() << endl; } 
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

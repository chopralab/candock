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
//~ #include "pdbreader/idatm.hpp"
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
		inout::output_file("", cmdl.top_seeds_file()); // output docked & filtered fragment poses
		inout::output_file("", cmdl.docked_ligands_file()); // output docked molecule conformations
		inout::output_file("", cmdl.mini_ligands_file()); // output docked & minimized ligands conformations
		inout::output_file(OMMIface::print_energies_title(), 
			cmdl.energy_file()); // output energy components of docked & minimized ligands conformations
		inout::output_file("", cmdl.nosql_file()); // probis local structural alignments
		
		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		vector<common::Centroid> centroids;
		if (cmdl.centroid_file().empty()) {
			probis::compare_against_bslib(argc, argv, cmdl.receptor_file(), 
				cmdl.receptor_chain_id(), cmdl.bslib_file(), cmdl.ncpu(),
				cmdl.nosql_file(), cmdl.json_file());
			genclus::generate_clusters_of_ligands(cmdl.json_file(), cmdl.json_with_ligs_file(),
				cmdl.geo_dir(), cmdl.names_dir(), cmdl.neighb(), cmdl.probis_clus_rad(),
				cmdl.probis_min_pts());
			const genlig::BindingSiteClusters binding_site_clusters = 
				genlig::generate_binding_site_prediction(cmdl.json_with_ligs_file(), 
				cmdl.bio_dir(), cmdl.num_bsites());
			inout::output_file(binding_site_clusters, cmdl.lig_clus_file());
			centroids = common::set_centroids(binding_site_clusters, cmdl.def_radial_check());	
		} else { // ... or else set binding sites from file
			centroids = common::set_centroids(cmdl.centroid_file(), 
				cmdl.def_radial_check(), cmdl.num_bsites());
		}
		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), Molib::PDBreader::first_model);
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
		Molib::MolGrid gridrec(receptors[0].get_atoms(cmdl.receptor_chain_id(), 
			Molib::Residue::protein));

		/* Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_forcefield_file(cmdl.amber_xml_file());
		receptors[0].prepare_for_mm(ffield, gridrec);

		/* Create gridpoints for each binding site represented by a centroid
		 * 
		 */
		vector<Geom3D::PointVec> gridpoints;
		for (auto &centroid : centroids) {
			gridpoints.push_back(common::identify_gridpoints(receptors[0], 
				centroid.get_centroid(), gridrec, centroid.get_radial_check(), 
				cmdl.grid_spacing(), cmdl.dist_cutoff(), cmdl.excluded_radius(), 
				cmdl.max_interatomic_distance()));
			inout::output_file(gridpoints.back(), cmdl.gridpdb_hcp_file(), ios_base::app);
		}

		/* Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */
		Molib::Molecules seeds;
		set<string> added;
		set<int> ligand_idatm_types;
		while(1 != 0) {
			Molib::Molecules ligands = lpdb.parse_molecule();
			if (ligands.empty()) break;

			// Compute properties, such as idatm atom types, fragments, seeds,
			// rotatable bonds
			ligands.compute_idatm_type().compute_rotatable_bonds()
				.compute_seeds(cmdl.seeds_file());
			ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);
			common::create_mols_from_seeds(added, seeds, ligands);
			inout::output_file(ligands, cmdl.prep_file(), ios_base::app);
		}
		dbgmsg(seeds);
		
		/* Read distributions file and initialize scores
		 * 
		 */
		Molib::Score score(Molib::get_idatm_types(receptors), ligand_idatm_types, 
			gridrec, cmdl.ref_state(),cmdl.comp(), cmdl.rad_or_raw(), 
			cmdl.dist_cutoff(), cmdl.distributions_file(), cmdl.step_non_bond());

		Molib::NRset non_clashing_seeds;
		map<int, Molib::Molecules*> m_non_clashing_seeds;
		vector<thread> threads;
		std::mutex mutex;

		/* Go over the predicted (or manually set) binding sites
		 * 
		 */
		for (auto &gpoints : gridpoints) {
			
			/* Create energy grids for all atom types that appear in small molecules
			 * 
			 */
			Molib::AtomTypeToEnergyPoint attep = score.compute_energy_grid(
				ligand_idatm_types, gpoints);
			common::HCPoints hcp = common::filter_scores(attep, cmdl.top_percent());
			
			inout::output_file(attep, cmdl.egrid_file(), ios_base::app); // output energy grid

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
					thread([&non_clashing_seeds, &m_non_clashing_seeds, &seeds, &gridrec, &score, &hcp, &mutex, i] () {
						// iterate over docked seeds and dock unique seeds
						for (int j = i; j < seeds.size(); j+= cmdl.ncpu()) {
							// make product graph between seed and hcp grid, cutoffs top_percent ?
							common::ProductGraph graph = 
								common::product_graph(hcp, seeds[j], cmdl.grid_spacing()); 
							common::ProductGraph::Cliques maxclq = 
								graph.max_weight_clique(cmdl.num_iter());

							/* Superimpose fragment coordinates onto clique
							 * and filter out those that clash with the receptor
							 *
							 */
							Molib::Molecules rot_seeds = 
								common::filter_clashes(common::superimpose(maxclq, seeds[j]), 
										gridrec);
							{
								// Guard non_clashing_seeds as it is the only variable that can 
								// get changed concurrently by different threads
								std::lock_guard<std::mutex> guard(mutex);

								// Add non-clashing docked seeds to seeds that
								// may have been found previously in other centroids
								if (!m_non_clashing_seeds.count(j)) {
									m_non_clashing_seeds[j] = &non_clashing_seeds.add(
										new Molib::Molecules());
								}
								for (auto &seed : rot_seeds) {
									m_non_clashing_seeds[j]->add(new Molib::Molecule(seed));
								}
							}
						}
				}));
			}
			for(auto& thread : threads) {
				thread.join();
			}
		}
		inout::output_file(non_clashing_seeds, cmdl.docked_seeds_file(), ios_base::app); // output docked & filtered fragment poses

		/* Main loop B with threading support : cluster non clashing seeds 
		 * based on rmsd and find best-scored cluster representatives
		 * 
		 */
		Molib::NRset top_seeds;
		threads.clear();
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(
				thread([&non_clashing_seeds, &top_seeds, &score, &mutex, i] () {
					// iterate over non clashing seeds and cluster
					for (int j = i; j < non_clashing_seeds.size(); j+= cmdl.ncpu()) {
						cluster::Clusters<Molib::Molecule> representatives = 
							common::cluster_molecules(non_clashing_seeds[j], score,
							cmdl.clus_rad(), cmdl.min_pts(), cmdl.max_num_clus(), 
							cmdl.max_seeds_to_cluster());
						common::get_representatives(top_seeds.add(new Molib::Molecules()),
							representatives);
#ifndef NDEBUG
						cluster::MapD<Molib::Molecule> scores = 
							score.many_ligands_score(top_seeds.last());
						dbgmsg(scores);
#endif
					}
				}));
		}
		for(auto& thread : threads) {
			thread.join();
		}
		inout::output_file(top_seeds, cmdl.top_seeds_file(), ios_base::app); // output docked & filtered fragment poses
		
		non_clashing_seeds.clear();	
		seeds.clear();

		/* Main loop C with threading support : connection of seeds and
		 * minimization
		 * 
		 */
		lpdb.rewind();
		lpdb.set_flags(Molib::PDBreader::all_models|Molib::PDBreader::hydrogens);
		while (1 != 0) {
			Molib::Molecules ligands = lpdb.parse_molecule();
			if (ligands.empty()) break;
			
			// again, compute idatm types etc... 
			// CORRECTED (correct in compute_fragment: seeds is not needed here)
			ligands.compute_idatm_type()
				.compute_hydrogen()
				.compute_bond_order()
				.compute_ring_type()
				.compute_bond_gaff_type()
				.compute_gaff_type()
				.erase_hydrogen()
				.compute_rotatable_bonds()
				.compute_overlapping_rigid_segments();

			/* Forcefield stuff : create forcefield for small molecules (and KB 
			 * non-bonded with receptor) and read receptor's forcefield xml file(s) into 
			 * forcefield object
			 * 
			 */
			OMMIface::ForceField ffield;
			ffield.parse_gaff_dat_file(cmdl.gaff_dat_file())
				.add_residue_topology(ligands)
				.add_kb_forcefield(score, cmdl.step_non_bond(), cmdl.scale_non_bond())
				.parse_forcefield_file(cmdl.amber_xml_file());

			threads.clear();
			for(int i = 0; i < cmdl.ncpu(); ++i) {
				threads.push_back(
					thread([&ligands, &receptors, &top_seeds, &gridrec, &score, &ffield, &mutex, i] () {
						// iterate over docked seeds and dock unique seeds
						for (int j = i; j < ligands.size(); j+= cmdl.ncpu()) {
							// if one ligand fails docking of others shall continue
							try {
								auto &ligand = ligands[j];
								cerr << "Quack quack ligand: " << endl << ligand << endl;
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
								if (docked.empty())	continue;
								
								inout::output_file(docked, cmdl.docked_ligands_file(), 
									ios_base::app); // output docked molecule conformations
		
								/* Cluster docked conformations and take only cluster
								 * representatives for further minimization
								 * 
								 */
								cluster::Clusters<Molib::Molecule> representatives = 
									common::cluster_molecules(docked, score,
									cmdl.docked_clus_rad(), cmdl.docked_min_pts(), cmdl.docked_max_num_clus());
								Molib::Molecules docked_representatives; 
								common::get_representatives(docked_representatives, representatives);
								
								/* Minimize each representative docked ligand conformation
								 * with full flexibility of both ligand and receptor
								 * 
								 */
								Molib::Molecules mini;
								OMMIface::Energies energies;
		
								for (int k = 0; k < docked_representatives.size(); ++k) {
									auto &docked_ligand = docked_representatives[k];
									OMMIface::OMM omm(receptors[0], docked_ligand, ffield, 
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
							catch (Error& e) { cerr << "skipping ligand due to : " << e.what() << endl; } 
							catch (out_of_range& e) { 
								cerr << "skipping ligand due to : " << e.what() << endl; 
								cerr << "This is the ligand that failed: " << endl << ligands[j] << endl;
							} 
						}
				}));
			}
			for(auto& thread : threads) {
				thread.join();
			}
		}
		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}

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
		
		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		vector<common::Centroid> centroids;
		if (cmdl.centroid_file().empty()) {
			throw Error("For testing use --centroid option to provide a centroid file");
		} else { // ... or else set binding sites from file
			centroids = common::set_centroids(cmdl.centroid_file(), 
				cmdl.def_radial_check(), cmdl.num_bsites());
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

			ligand_idatm_types = Molib::get_idatm_types(ligands, ligand_idatm_types);
			common::create_mols_from_seeds(added, seeds, ligands);
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

		/* Go over the predicted (or manually set) binding sites and filter
		 * top scores
		 * 
		 */
		vector<common::HCPoints> hcp_vec;
		for (auto &gpoints : gridpoints) {
			
			/* Create energy grids for all atom types that appear in small molecules
			 * 
			 */
			Molib::AtomTypeToEnergyPoint attep = score.compute_energy_grid(
				ligand_idatm_types, gpoints);
			hcp_vec.push_back(common::filter_scores(attep, cmdl.top_percent()));

			inout::output_file(attep, cmdl.egrid_file(), ios_base::app); // output energy grid
			dbgmsg("after output energy grid");

			/* Create template grids using ProBiS-ligands algorithm
			 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
			 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
			 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
			 */
		 }

		vector<thread> threads;
		//~ std::mutex my_mutex;

		//~ Molib::NRset top_seeds;

		threads.clear();
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(
				thread([&seeds, &gridrec, &score, &hcp_vec, i] () {
					// iterate over docked seeds and dock unique seeds
					for (int j = i; j < seeds.size(); j+= cmdl.ncpu()) {
						dbgmsg(seeds[j]);
						// iterate over top scores for each centroid
						Molib::Molecules non_clashing_seeds;
						for (auto &hcp : hcp_vec) {
							
							// make product graph between seed and hcp grid, cutoffs top_percent ?
							common::ProductGraph graph = 
								common::product_graph(hcp, seeds[j], cmdl.grid_spacing()); 
							common::ProductGraph::Cliques maxclq = 
								graph.max_weight_clique(cmdl.num_iter());

							// Superimpose fragment coordinates onto clique and 
							// filter out those that clash with the receptor
							for (auto &molecule : common::filter_clashes(
								common::superimpose(maxclq, seeds[j]), gridrec)) {
								non_clashing_seeds.add(new Molib::Molecule(molecule));
							}
						}
						inout::output_file(non_clashing_seeds, cmdl.docked_seeds_file(), ios_base::app); // output docked & filtered fragment poses
						// cluster non clashing seeds based on rmsd and 
						// find best-scored cluster representatives
						auto clusters_reps_pair = common::cluster_molecules(
							non_clashing_seeds, score, cmdl.clus_rad(), cmdl.min_pts(), 
							cmdl.max_num_clus(), cmdl.max_seeds_to_cluster());
						
						Molib::Molecules tseeds;
						common::convert_clusters_to_mols(tseeds, clusters_reps_pair.second);
#ifndef NDEBUG
						cluster::MapD<Molib::Molecule> scores = 
							score.many_ligands_score(tseeds);
						dbgmsg(scores);
#endif
						inout::output_file(tseeds, "tmp/" + seeds[j].name() + "/" 
							+ cmdl.top_seeds_file()); // output docked & filtered fragment poses
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

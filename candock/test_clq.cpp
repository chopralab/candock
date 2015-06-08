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
		inout::output_file("", cmdl.docked_seeds_file()); // output docked & filtered fragment poses

		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		common::Centroids centroids;
		if (cmdl.centroid_in_file().empty()) {
			throw Error("For testing use --centroid option to provide a centroid file");
		} else { // ... or else set binding sites from file
			centroids = common::set_centroids(cmdl.centroid_in_file());
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

		/* Create gridpoints for ALL centroids representing one or more binding sites
		 * 
		 */
		Geom3D::GridPoints gridpoints = common::identify_gridpoints(centroids, 
			gridrec, cmdl.grid_spacing(), cmdl.dist_cutoff(), 
			cmdl.excluded_radius(), cmdl.max_interatomic_distance());
		inout::output_file(gridpoints, cmdl.gridpdb_hcp_file());

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

		/* Create energy grids for all atom types that appear in small molecules
		 * for the combined binding site
		 */
		Molib::AtomTypeToEnergyPoint attep = score.compute_energy_grid(
			ligand_idatm_types, gridpoints);
		common::HCPoints hcp = common::filter_scores(attep, cmdl.top_percent());

		inout::output_file(attep, cmdl.egrid_file()); // output energy grid
		dbgmsg("after output energy grid");

		/* Create template grids using ProBiS-ligands algorithm
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS 
		 */

		vector<thread> threads;

		threads.clear();
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(
				//~ thread([&seeds, &gridrec, &score, &hcp, i] () {
				thread([&seeds, &gridrec, &score, &gridpoints, i] () {
					// iterate over seeds and dock unique seeds
					for (int j = i; j < seeds.size(); j+= cmdl.ncpu()) {
						dbgmsg(seeds[j]);
						
						/* Dock seed to the entire grid
						 *
						 */ 
						Molib::Molecules non_clashing_seeds = 
							//~ common::dock_seeds(hcp, seeds[j], cmdl.grid_spacing());
							common::dock_seeds(gridpoints, seeds[j], cmdl.grid_spacing());
						
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

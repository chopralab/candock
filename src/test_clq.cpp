#include <iostream>
#include <exception>
#include <typeinfo>
#include <thread>
#include <mutex>
#include "program/cmdlnopts.hpp"
#include "helper/benchmark.hpp"
#include "helper/inout.hpp"
#include "helper/path.hpp"
#include "helper/error.hpp"
#include "molib/grid.hpp"
#include "molib/molecules.hpp"
#include "molib/pdbreader.hpp"
#include "modeler/forcefield.hpp"
#include "score/score.hpp"
#include "docker/gpoints.hpp"
#include "docker/conformations.hpp"
#include "docker/dock.hpp"
#include "centro/centroids.hpp"
#include "program/common.hpp"

using namespace std;
using namespace Program;

int main(int argc, char* argv[]) {
	try {
		CmdLnOpts cmdl;
		cmdl.init(argc, argv, CmdLnOpts::FRAG_DOCKING | CmdLnOpts::FORCE_FIELD | CmdLnOpts::SCORING);
		cmdl.display_time("started");
		cout << cmdl << endl;
		
		/* Identify potential binding sites using ProBiS algorithm
		 * or alternatively set binding sites from file
		 * 
		 */
		Centro::Centroids centroids;
		if ( inout::Inout::file_size(cmdl.centroid_file()) <= 0 ) {
			throw Error("Provided centroid file not found!");
		} else { // ... or else set binding sites from file
			centroids = Centro::set_centroids(cmdl.centroid_file(), cmdl.num_bsites());
		}

		/* Initialize parsers for receptor (and ligands) and read
		 * the receptor molecule(s)
		 * 
		 */
		Molib::PDBreader rpdb(cmdl.receptor_file(), 
			Molib::PDBreader::first_model);
		Molib::Molecules receptors = rpdb.parse_molecule();
		

		Molib::PDBreader lpdb(cmdl.prep_file(), Molib::PDBreader::all_models, 
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

		/* Prepare receptor for molecular mechanics: histidines, N-[C-]terminals,
		 * bonds, disulfide bonds, main chain bonds
		 * 
		 */
		OMMIface::ForceField ffield;
		ffield.parse_forcefield_file(cmdl.amber_xml_file())
			.parse_forcefield_file(cmdl.water_xml_file());

		receptors[0].prepare_for_mm(ffield, gridrec);

		/* Read ligands from the ligands file - this file may contain millions
		 * of ligands, and we read only a few at one time, to save memory
		 * 
		 */
		Molib::Molecules seeds;
		set<int> added;
		set<int> ligand_idatm_types;

		Molib::Molecules ligands;
		while(lpdb.parse_molecule(ligands)) {
			ligand_idatm_types = ligands.get_idatm_types(ligand_idatm_types);
			common::create_mols_from_seeds(added, seeds, ligands);
			ligands.clear();
		}
#ifndef NDEBUG
		for (auto &seed : seeds) {
			inout::output_file(seed, "seed_" + seed.name() + ".pdb"); 
		}
#endif
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

		dbgmsg("START SCORE" << endl << score << "END SCORE");
		
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

		vector<thread> threads;

		threads.clear();
		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back(
#ifndef NDEBUG
				thread([&seeds, &gridrec, &score, &gpoints0, &gpoints, &cmdl, i] () {
#else
				thread([&seeds, &gpoints0, &gpoints, i, &cmdl] () {
#endif
					// iterate over seeds and dock unique seeds
					for (size_t j = i; j < seeds.size(); j+= cmdl.ncpu()) {
						try {
							dbgmsg(seeds[j]);

#ifndef NDEBUG
							Molib::Molecule orig_seed = seeds[j];
#endif
							
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

							inout::output_file(dock.get_docked(), Path::join(Path::join(cmdl.top_seeds_dir(), dock.get_docked().name()), 
								cmdl.top_seeds_file())); // output docked & clustered seeds

#ifndef NDEBUG
							int nn = 0;
							for (auto &d : dock.get_docked()) {
								cout << "SEED" << endl << orig_seed << endl << "DOCKED" << endl << d << endl;
								cout << "seed=" << dock.get_docked().name() << " conf=" << ++nn << " rmsd=" << d.compute_rmsd_ord(orig_seed) << " energy=" << d.name() << endl;
							}
#endif
							
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

		cmdl.display_time("finished");
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}

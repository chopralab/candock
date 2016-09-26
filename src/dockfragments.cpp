#include "dockfragments.hpp"

#include "helper/path.hpp"

namespace Program {

	bool DockFragments::__can_read_from_files ( const CmdLnOpts& cmdl )
	{
		bool all_seeds_are_present = true;

		const Molib::Molecules& all_seeds = __fragmented_ligands.seeds();

		// No early return so that we have the ability to redock missing seeds
		for (const auto& seed : all_seeds) {
			all_seeds_are_present &= inout::Inout::file_size( __name + "_" + Path::join( Path::join(cmdl.top_seeds_dir(), seed.name()),
			                                                  cmdl.top_seeds_file() ) ) > 0;
		}
		
		return all_seeds_are_present;
	}

	void DockFragments::__read_from_files ( const CmdLnOpts& cmdl )
	{
		cout << "All seeds are present in " << cmdl.top_seeds_dir() << " for " << __name << ". Docking of fragments skipped." << endl;
	}

	void DockFragments::__dock_fragment ( int start, const Docker::Gpoints& gpoints, const Docker::Gpoints& gpoints0, const CmdLnOpts& cmdl) {
		// iterate over docked seeds and dock unique seeds
		for (int j = start; j < __fragmented_ligands.seeds().size(); j+= cmdl.ncpu()) {
			try {
				
				int file_exist_check = inout::Inout::file_size( __name + "_" + Path::join( Path::join(cmdl.top_seeds_dir(), __fragmented_ligands.seeds()[j].name()),
			                                                  cmdl.top_seeds_file() ) ) > 0;
				
				if ( file_exist_check > 0 ) {
					cout << "Skipping docking of seed: " << __fragmented_ligands.seeds()[j].name() << " because it is already docked!" << endl;
					continue;
				}
				dbgmsg(__fragmented_ligands.seeds()[j]);
				/* Compute all conformations of this seed with the center
				 * atom fixed on coordinate origin using maximum clique algorithm
				 * 
				 */
				Docker::Conformations conf(__fragmented_ligands.seeds()[j], gpoints0, 
						                   cmdl.conf_spin(), // degrees
						                   cmdl.num_univec() // number of unit vectors
									  );

#ifndef NDEBUG
				inout::output_file(conf, "conf_" + __fragmented_ligands.seeds()[j].name() + ".pdb"); 
#endif			
				/* Dock this seed's conformations to the entire grid by moving them 
				 * over all gridpoints and probe where they clash with the receptor: 
				 * cluster docked conformations based on rmsd and energy and output 
				 * only best-scored cluster representatives
				 *
				 */
#ifndef NDEBUG
				Docker::Dock dock(gpoints, conf, __fragmented_ligands.seeds()[j], __score, __gridrec, cmdl.clus_rad());
#else
				Docker::Dock dock(gpoints, conf, __fragmented_ligands.seeds()[j], cmdl.clus_rad());
#endif

				dock.run();

				inout::output_file(dock.get_docked(), __name + "_" + Path::join(Path::join(cmdl.top_seeds_dir(), dock.get_docked().name()),
					cmdl.top_seeds_file())); // output docked & clustered seeds
	
			}
			catch (exception& e) {
				cerr << "skipping seed due to : " << e.what() << endl;
			}
		}
	}

	void DockFragments::__continue_from_prev ( const CmdLnOpts& cmdl )
	{

		cout << "Docking fragments into: " << __name << "_" << cmdl.top_seeds_dir() << ". Files will be named: " << cmdl.top_seeds_file() << endl;

		/* Create gridpoints for ALL centroids representing one or more binding sites
		 * 
		 */
		Docker::Gpoints gpoints(__score, __fragmented_ligands.ligand_idatm_types(), __found_centroids.centroids(),
			__gridrec, cmdl.grid_spacing(), cmdl.dist_cutoff(), cmdl.excluded_radius(), 
			cmdl.max_interatomic_distance());
		inout::output_file(gpoints, __name + "_" + cmdl.gridpdb_hcp_file());

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

		std::vector<std::thread> threads;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back( std::thread([&,this,i] {__dock_fragment(i, gpoints, gpoints0, cmdl);} ) );
		}
		for(auto& thread : threads) {
			thread.join();
		}
		
		cout << "Done with fragment docking" << std::endl;

	}

	std::vector<std::pair<double, std::string>> DockFragments::get_best_seeds(const CmdLnOpts& cmdl) {
		const Molib::Molecules& all_seeds = __fragmented_ligands.seeds();

		std::vector<std::pair<double, std::string>> seed_score_map;

		for ( const auto &seed : all_seeds ) {
			Molib::PDBreader spdb (__name + "_" + Path::join( Path::join(cmdl.top_seeds_dir(), seed.name()),
			                                                  cmdl.top_seeds_file()), Molib::PDBreader::first_model);

			Molib::Molecules seed_molec = spdb.parse_molecule();
			seed_score_map.push_back( {std::stod( seed_molec.first().name()), seed.name()} );
		}

		std::sort(seed_score_map.begin(), seed_score_map.end() );

		return seed_score_map;
	}

}

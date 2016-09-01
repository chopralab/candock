#include "dockfragmentsstep.hpp"

#include "helper/path.hpp"

namespace Program {

	bool DockFragmentsStep::__can_read_from_files ( const CmdLnOpts& cmdl )
	{
		return true;
	}

	void DockFragmentsStep::__read_from_files ( const CmdLnOpts& cmdl )
	{

	}

	void DockFragmentsStep::__dock_fragment ( int start, const Docker::Gpoints& gpoints, const Docker::Gpoints& gpoints0, const CmdLnOpts& cmdl) {
		std::cout << "Fragment Docking has begun" << std::endl;
		// iterate over docked seeds and dock unique seeds
		for (int j = start; j < __fragmented_ligands.seeds().size(); j+= cmdl.ncpu()) {
			try {
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

				//~ inout::output_file(dock.get_docked(), cmdl.top_seeds_dir() + "/" + dock.get_docked().name() + "/" 
				//~ + cmdl.top_seeds_file()); // output docked & clustered seeds
				inout::output_file(dock.get_docked(), Path::join(Path::join(cmdl.top_seeds_dir(), dock.get_docked().name()),
					cmdl.top_seeds_file())); // output docked & clustered seeds
	
			}
			catch (exception& e) {
				cerr << "skipping seed due to : " << e.what() << endl;
			}
		}
	}

	void DockFragmentsStep::__continue_from_prev ( const CmdLnOpts& cmdl )
	{
		std::cout << "ASFASFASDFSA" << std::endl;
		/* Create gridpoints for ALL centroids representing one or more binding sites
		 * 
		 */
		Docker::Gpoints gpoints(__score, __fragmented_ligands.ligand_idatm_types(), __found_centroids.centroids(),
			__gridrec, cmdl.grid_spacing(), cmdl.dist_cutoff(), cmdl.excluded_radius(), 
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

		std::vector<std::thread> threads;

		for(int i = 0; i < cmdl.ncpu(); ++i) {
			threads.push_back( std::thread([&,this,i] {__dock_fragment(i, gpoints, gpoints0, cmdl);} ) );
		}
		for(auto& thread : threads) {
			thread.join();
		}
		
		cout << "Done with fragment docking" << std::endl;

	}

}
